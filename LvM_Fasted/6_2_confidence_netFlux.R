# This script calculates the confidence interval for the net flux of a pair of exchange fluxes
# for each reaction, different tracers are processed in series, not in parallel.

library(readxl)
library(matrixcalc)
library(Ryacas) # combine like terms
library(tidyverse)
library(parallel)


rm(list = ls())

time.start <- Sys.time()

# set path to the folder of the current active script
path <- rstudioapi::getActiveDocumentContext()$path %>% dirname()
setwd(path); getwd()

# cleaned labeling - infusion data (from the 'data' folder)
load(file = "../data/cleaned_labeling_data.RData")


state_fasted.or.refed <- "fasted"


if (state_fasted.or.refed == "fasted") {
  load("3_list.stoich.EMU_allTracers_fasted.RData") # fast-specific EMU decomposed result 
  d.13C.labeling <- d.13C.fasted
  totalCO2       <- 1800 
}



# Import supplementary functions
load("1_Supplement_functions.RData")

# Import the optimal flux solution from the repeated MFA runs
load("5_2_optimal_solution.RData")




time.Start.CI <- Sys.time()


# tissues included in the EMU decomposition model
TissuesModeled <- target.EMU %>% str_extract("(?<=\\.)[:alpha:]{1,9}(?=\\_)") %>% unique()


# Select labeling data with the same scope of tissues included in the EMU decomposition model
d.13C.labeling  <- d.13C.labeling  %>% filter(tissue %in% TissuesModeled)




# make sure the modeled tissues are also included in the labeling dataset
if ( (sum(TissuesModeled %in% unique(d.13C.labeling $tissue)) != n_distinct(d.13C.labeling $tissue)) ) error("Modeled tissues are not included in the labeling dataset.")

# in the labeling dataset, include only compound-tissue that is included in the target EMU list of the model network
d.13C.labeling  <- d.13C.labeling  %>% filter(Compound.tissue.seq %in% target.EMU)



# observed labeling
tracers <- names(l.stoich.EMU_allTracers) # [1:4]


# reaction index to calculate confidence interval


# specify a pair of fluxes, whose net fluxes is to be determined 

# Lac <-> Pyr
f1 <- reaction_data %>% filter(reactions == "Lac.Blood->Pyr.Lvc") %>% pull(enzyme); f1
f2 <- reaction_data %>% filter(reactions == "Pyr.Lvc->Lac.Blood") %>% pull(enzyme); f2


# Gln <-> aKG 
f5 <- reaction_data %>% filter(reactions == "Gln.Blood->aKG.Lv") %>% pull(enzyme); f5
f6 <- reaction_data %>% filter(reactions == "aKG.Lv->Gln.Blood") %>% pull(enzyme); f6


# Ala <-> Pyr
f9  <- reaction_data %>% filter(reactions == "Ala.Blood->Pyr.Lvc") %>% pull(enzyme); f9
f10 <- reaction_data %>% filter(reactions == "Pyr.Lvc->Ala.Blood") %>% pull(enzyme); f10

# AcAct <-> HB
f11  <- reaction_data %>% filter(reactions == "AcAct.Lv->HB.Blood") %>% pull(enzyme); f11
f12  <- reaction_data %>% filter(reactions == "HB.Blood->AcAct.Lv") %>% pull(enzyme); f12


f.all <- list(c(f1, f2), c(f5, f6), c(f9, f10), c(f11, f12))


# Create a cluster
cl <- makeCluster( 
  # if mac, use 6 cores at max; otherwise if lab computer use 20 cores
  ifelse(.Platform$OS.type == "unix", min(length(f.all), 6), 20), 
  # if mac, print progress in console; otherwise if lab computer windows save progress in txt file
  outfile = ifelse(.Platform$OS.type == "unix", "", "log_CI.txt") )  




# Make libraries accessible to each working node
clusterEvalQ(cl, {
  library(parallel)  
  library(readxl)
  library(Ryacas)
  library(tidyverse)
})



# foreach (f = c(11, 17, 23:26, 28, 32)) %dopar% { # loop through all flux index

func.calculate_confidence_interval <- function(f){
  # for (f in c(10)){  
  
  # check 
  # f=23
  
  
  # The reaction whose C.I. we seek to determine is fixed to a given level, e.g. increase by 100% of the optimal value
  flux.C.I_index1 <- f[1] # reaction index of each of the flux pair
  flux.C.I_index2 <- f[2]
  
  # optimal value of the net flux
  flux.C.I_value.optimal <- optimal_allFlux[flux.C.I_index1] - optimal_allFlux[flux.C.I_index2] 
  
  # determine the adjustment step h
  if ( round(flux.C.I_value.optimal) == 0) {
    h0 <- 50} else {
      h0 <- flux.C.I_value.optimal # change by 100% incremental at the initial setp
    }
  
  h     <- h0
  stage <- "UB" # upper bound of C.I. 
  
  # record the flux value and the associated cost
  d.flux.cost <- tibble(flux.index = paste0(flux.C.I_index1, "-", flux.C.I_index2),
                        flux.value = flux.C.I_value.optimal,
                        cost = optimal_cost)
  
  total.binary.steps <- 8 # make a max of 6 iterations of binary search for each side
  binary.index <- 1
  lock <- F
  
  while (binary.index <= total.binary.steps) { 
    
    message("net flux ", "v", f[1], "-", "v", f[2], "; stage ", stage, "; binary index ", binary.index)
    message(paste("h = ", round(h, 2)))
    
    
    # -------below is the standalone procedure for flux estimation using the classic algorithm  ------
    
    # measured IDV; note that the tracer is in order of 'l.stoich.EMU_allTracers'
    l.IDV.obs  <- list() # measured IDV; note that the tracer is in order of 'l.stoich.EMU_allTracers'; each tracer data is one element
    IDV.obs    <- c()    # the vector form of l.IDV.obs
    IDV_SD.obs <- c()    # the vector form of SD of labeling
    
    l.Compound.tissue.seq      <- list() # compound-tissue to extract repeatedly for difference mice to construct the sensitivity matrix; for extraction elements of list
    l.Compound.tissue.seq.m.id <- list() # same list as ablove, but adding mouse id; for naming purpose
    
    l.13C <- list() # a list of labeling data
    
    
    for (b in seq(tracers)){ # b the index
      # print(b)
      # check
      # b=9
      # b=5
      a <- tracers[b]        # a the tracer name
      # print(a)
      
      # metabolite measured labeling
      d.a <- d.13C.labeling %>% filter(Infusate == a)
      
      # Create a list, with each element corresponding to a tracer. Each element is a vector of compound-tissue|mouse replicate used to extract simulated labeling and derivatives
      # Choose only compound-tissue that is in the targeted EMU list
      Compound.tissue.seq.m.id_a <- d.a %>% 
        select(mouse.when.who, Compound.tissue.seq_m.id) %>% distinct() %>% pull(Compound.tissue.seq_m.id)
      
      l.Compound.tissue.seq.m.id[[b]] <- Compound.tissue.seq.m.id_a                                 # for naming elements,                 with    mouse ID
      
      
      
      
      ### The mouse id is in the format of m1, m2...as in Tony's data, or just letters a, b, c, d...without prefix m 
      # 1) For Tony's data only
      # l.Compound.tissue.seq     [[b]] <- Compound.tissue.seq.m.id_a %>% unique() %>% str_remove("\\|m\\d{1,2}$") # for extracting elements from a list, without mouse ID 
      
      # 2) inclusive format for both Tony data and 13C protein tracing data 
      l.Compound.tissue.seq     [[b]] <- Compound.tissue.seq.m.id_a %>% unique() %>% str_remove("\\|[a-zA-Z0-9]{1,3}$")
      
      
      
      names(l.Compound.tissue.seq.m.id)[b] <- a # for naming elements,                 with    mouse ID
      names(l.Compound.tissue.seq)     [b] <- a # for extracting elements from a list, without mouse ID 
      
      # IDV observed
      IDV.obs.a <- d.a$labeling
      names(IDV.obs.a) <- d.a %>% mutate(name = str_c(Infusate, "-", Compound.tissue.seq, "_#", C_Label, "|", m.id)) %>% pull(name)
      
      # IDV SD observed
      IDV_SD.obs.a <- d.a$labeling.sd 
      if (sum(is.na(IDV_SD.obs.a)) >= 1 ) stop("there is NA in standard deviations...")
      names(IDV_SD.obs.a) <- names(IDV.obs.a)
      
      # compile tracers in the form of a vector
      IDV.obs    <- c(IDV.obs,    IDV.obs.a)     # labeling
      IDV_SD.obs <- c(IDV_SD.obs, IDV_SD.obs.a)  # labeling SD
      
      # compile tracers in the form of a list, each tracer as an element
      l.IDV.obs[[b]] <- IDV.obs.a
      names(l.IDV.obs)[b] <- a
    }                  
    
    
    # IDV.covariance.inverse <- IDV.covariance %>% func.pseudo_inverse_svd() 
    
    names(IDV.obs)
    
    n.opt.steps <- 30 # maximum iterate optimization steps; 
    PRINT <- F # if True, progress at major calculation steps is printed
    
    for (optim.step in 1:n.opt.steps) { 
      
      # cat("\n")
      # message("------------------iteration ", optim.step, " ------------------------------------------------------")
      # check
      # optim.step <- 1
      if (optim.step == 1){
        
        # option (1) initiate random flux, s.t. all fluxes ≥ 0
        u.initial <- optimal_free_flux
        # names(u.initial) <- paste0("v", names(u.initial))
        
        u <- u.initial # free fluxes
        
        # Replace a single flux notation with tangible flux value; 
        # define this appears important in parallel computation; otherwise u can't be found
        func.FillwithFluxValue <- function(x, whichFlux = "v10"){
          str_replace_all(x, pattern = paste0("\\b", whichFlux, "\\b"), replacement = u[whichFlux] %>% as.character() )
        }
      }
      
      IDV.simulated_allTracers <- NULL # collecting the               simulated glucose labeling from all tracer experiments
      dX.target_allTracers     <- NULL # collecting the derivative of simulated glucose labeling from all tracer experiments; dX is same as derivative of IDV
      
      # check
      # tracer.k="13CGlc"
      # tracer.k="13CLac"
      # tracer.k="13CGlycerol"
      for (tracer.k in tracers){
        message("net flux ", "v", f[1], "-v", f[2], "; stage ", stage, "; binary index ", binary.index, "; optim.step ", optim.step,   "; tracer ", tracer.k)
        
        # load constant quantities unique to tracer k, yet globally shared in all optimization iterations for tracer k
        # the quantities below are pre-calculated in "Stoic_EMU_calculator.R"
        l.stoich.EMU.k         <- l.stoich.EMU_allTracers[[tracer.k]]
        l.AXBY                 <- l.stoich.EMU.k  $  l.AXBY
        l.EMU.substrates       <- l.stoich.EMU.k  $  l.EMU.substrates
        l.deriv_EMU.substrates <- l.stoich.EMU.k  $  l.deriv_EMU.substrates
        l.deriv_A              <- l.stoich.EMU.k  $  l.deriv_A
        l.deriv_B              <- l.stoich.EMU.k  $  l.deriv_B
        
        # the sensitive matrix of the target compound
        dX.target <- NULL
        
        # the lists of EMU source and its derivatives will extend to include more EMUs while looping through EMU sizes,
        # and is updated after each optimization iteration
        l.EMU.sources <- l.EMU.substrates
        l.deriv_EMU.sources <- l.deriv_EMU.substrates
        
        for (j in 1:n.EMU_sizes) { # note that j itself is not necessarily the same as the actual EMU size, as EMU size may not be continuous
          # j=1
          # j=4
          
          # Note that A and B are subnetwork (EMU size dependent)
          # For a given EMU size sub-network, A and B are functions of fluxes, and gets updated after each iteration
          A <- l.AXBY[[j]]$A %>% as.matrix()
          B <- l.AXBY[[j]]$B %>% as.matrix()
          
          # replace A and B flux notation with real numbers
          if (nrow(A) == 1 | nrow(B) == 1) {
            for (a in paste0("v", index.freeFluxes) ) { # loop through free fluxes
              A <- A %>% as.matrix() %>% apply(MARGIN = c(1,2), func.FillwithFluxValue, whichFlux = a)
              B <- B %>% as.matrix() %>% apply(MARGIN = c(1,2), func.FillwithFluxValue, whichFlux = a)
            }
          } else { # faster speed
            for (a in paste0("v", index.freeFluxes) ) { # loop through v1, v2, v3...
              A <- A %>% as.matrix() %>% apply(MARGIN = 2, func.FillwithFluxValue, whichFlux = a)
              B <- B %>% as.matrix() %>% apply(MARGIN = 2, func.FillwithFluxValue, whichFlux = a)
            }
          }
          # parse and evaluate the flux into numeric values
          A <- A %>% as.matrix() %>% apply(MARGIN = c(1, 2), function(x){ if (x=="0") {return(0)} else { eval(parse(text = x)) }   } )
          B <- B %>% as.matrix() %>% apply(MARGIN = c(1, 2), function(x){ if (x=="0") {return(0)} else { eval(parse(text = x)) }   } )
          
          A.inverse <- func.pseudo_inverse_svd(A)
          
          
          X.names <- l.AXBY[[j]]$X # each element in X is a single EMU
          EMU.actualSize <- X.names[1] %>% str_extract("(?<=_)\\d+") %>% str_count() # count digits after underscore to get the size of the current EMU network
          Y.names <- l.AXBY[[j]]$Y # an element can be a single EMU or a Cauchy product
          
          # Calculate Y as a Cauchy product
          # In the first (smallest) EMU network, rows corresponding to substrate IDV
          # for Y of larger EMU network, it requires calculation of Cauchy product from the X matrix of prior iterations, or substrate EMU; the source individual EMU sizes can be different
          Y <- matrix(ncol = EMU.actualSize + 1, nrow = length(Y.names)) # empty matrix to be updated in the loop below
          
          for (iii in 1:length(Y.names)){
            
            # for single EMU, directly get the IDV from the EMU knowledge pool (list)
            if (! Y.names[iii] %>% str_detect("\\+") ) {
              Y [iii, ] <- l.EMU.sources      [Y.names[iii]] %>% as_tibble() %>% as.matrix() %>% t()
              
              # optional - checking point: labeling of all isotopologues of the given molecule should sum up to 1
              if (sum(Y [iii, ]) %>% round(3) != 1) stop("Isotopologues should sum up to 100%.")
              
            } else { # for EMU in Y being a Cauchy product
              Y.names.split <- Y.names[iii] %>% str_split_1(pattern = "\\+") # split Cauchy products into each individual EMU names
              x1 <- l.EMU.sources[ Y.names.split[1] ] %>% unlist() %>% unname() # 1st individual EMU
              x2 <- l.EMU.sources[ Y.names.split[2] ] %>% unlist() %>% unname() # 2nd individual EMU
              Cauchy <- convolve(x1, rev(x2), type = "open") %>% unname() # calculate the Cauchy product
              Y[iii, ] <- Cauchy # update the iii-th row of the Y
              
              # optional - checking point: labeling of all isotopologues of the given molecule should sum up to 1
              if (sum(Y [iii, ]) %>% round(3) != 1) stop("Isotopologues should sum up to 100%.")
              
            }
          }
          if (PRINT == T) print("Successfully calculated Y matrix")
          
          
          # Calculate X, which contains the IDV of newly produced EMU
          X <- A.inverse %*% B %*% Y
          rownames(X) <- X.names
          
          # optional - checking point: each row should sum up to 1
          for (xr in seq(nrow(X))){
            if (sum(X[xr, ]) %>% round(2) != 1) stop(paste("in row", xr, ": Isotopologues do not sum up to 100% as it should."))
          }
          
          # collect IDV of newly produced EMU (which is first converted from a matrix to a list to be bound to the compiling list)
          l.EMU.sources <- append(l.EMU.sources, func.turnMatrixToList(X))
          
          
          # Loop through free fluxes to calculate derivatives: dY and then dX
          for (flux.i in paste0("v", index.freeFluxes ) ){
            
            # check
            # flux.i="v1"
            # calculate derivatives of Y
            # In the first (smallest) EMU network, rows corresponding to substrate IDV derivatives
            # for Y of larger EMU network, it requires calculation of Cauchy product from derivatives of the X matrix of prior iterations, or substrate EMU; the source individual EMU sizes can be different
            dY <- matrix(ncol = EMU.actualSize + 1, nrow = length(Y.names)) # empty matrix to be updated in the loop below; this dY should be of the same dimension as Y calculated outside of the free flux loop
            
            for (iii in 1:length(Y.names)){
              
              # for single EMU, directly get the IDV from the EMU knowledge pool (list)
              if (! Y.names[iii] %>% str_detect("\\+") ) {
                dY[iii, ] <- l.deriv_EMU.sources[ paste0(Y.names[iii], "-", flux.i) ] %>% as_tibble() %>% as.matrix() %>% t()
                
              } else { # for EMU in Y being a Cauchy product
                Y.names.split <- Y.names[iii] %>% str_split_1(pattern = "\\+") # split Cauchy products into each individual EMU names
                # EMU
                x1 <- l.EMU.sources[ Y.names.split[1] ] %>% unlist() %>% unname() # 1st individual EMU
                x2 <- l.EMU.sources[ Y.names.split[2] ] %>% unlist() %>% unname() # 2nd individual EMU
                # derivatives
                dx1 <- l.deriv_EMU.sources[ paste0(Y.names.split[1], "-", flux.i) ] %>% unlist() %>% unname()
                dx2 <- l.deriv_EMU.sources[ paste0(Y.names.split[2], "-", flux.i) ] %>% unlist() %>% unname()
                
                # both calculation gives same result, but the first line is the right equation!
                dCauchy <- convolve(dx1, rev(x2), type = "open") + convolve(x1, rev(dx2), type = "open")
                # dCauchy <- convolve(dx1, rev(dx2), type = "open")
                dY[iii, ] <- dCauchy
              }
            }
            if (PRINT == T) print("Successfully calculated Y derivative matrix")
            
            # derivatives of A and B
            dA <- l.deriv_A[[ paste0("dA", j, "_", flux.i) ]]
            dB <- l.deriv_B[[ paste0("dB", j, "_", flux.i) ]]
            
            # The newly produced EMU derivatives from dX do NOT contain substrate information (as the substrate information is contained in matrix Y)
            # The new derivative EMU is named after the free flux postscript; the substrates derivatives (defined in the very beginning) has no subscript, as their derivatives are all zero, independent of fluxes
            # Thus the new dX can be directly appended to the original list of substrate derivatives without duplication
            dX <- A.inverse %*% (dB%*%Y + B%*%dY - dA%*%X)
            rownames(dX) <- paste0(X.names, "-", flux.i)
            
            # collect derivatives of newly produced EMU (which is first converted from a matrix to a list to be bound to the compiling list)
            l.deriv_EMU.sources <- append(l.deriv_EMU.sources, func.turnMatrixToList(dX))
            
            if (PRINT == T) print(paste("Calculated derivatives for flux:", flux.i))
          }
          if (PRINT == T) print(paste("Successfully calculated derivatives for all free fluxes within EMU size of", j))
        }
        
        if (PRINT == T) cat("Completed analysis under the tracer:", tracer.k, "\n\n")
        
        
        
        
        
        # testing: all mouse labeling data to be used in model, its predicted label must be available from 'l.EMU.sources'; 
        # otherwise, add the missing compound-tissue EMU in 'target.EMU' in 'Stoic_EMU_Calculator.R' 
        if (sum(!(l.Compound.tissue.seq[[tracer.k]] %in% names(l.EMU.sources))) > 0) stop ("For the given compound-tissue measured labeling, its predicted EMU is not available. Set it up in the target EMU.")
        
        # create a vector of simulated labeling
        positions <- match(l.Compound.tissue.seq[[tracer.k]], names(l.EMU.sources)) # index of a mouse-compound-tissue EMU in the compiled fitted EMU list
        IDV.simulated.k <- l.EMU.sources[positions] %>% unlist()
        
        # add names
        IDV.simulated.k.names  <- d.13C.labeling %>% filter(Infusate == tracer.k) %>% mutate(x = str_c(Infusate, "-", Compound.tissue.seq, "_#", C_Label, "|", m.id)) %>% pull(x)
        names(IDV.simulated.k) <- IDV.simulated.k.names
        
        # checkpoint: the fitted data length matches the measured data labeling length (optional)
        if (length(IDV.simulated.k) != length(IDV.simulated.k.names))             stop("Length mismatches")
        if (identical(names(l.IDV.obs[[tracer.k]]), names(IDV.simulated.k)) != T) stop("naming mismatches")
        
        # compile simulated IDV from all tracers 
        IDV.simulated_allTracers <- c(IDV.simulated_allTracers, IDV.simulated.k)
        
        if (PRINT == T) print(paste("constructing sensitivity matrix of ", tracer.k))
        
        s <- data.frame()
        for (i in seq( l.Compound.tissue.seq[[tracer.k]] ) ) {
          
          # which compound-tissue (of a mouse replicate) to look for
          w <- l.Compound.tissue.seq[[tracer.k]][i] 
          # the sensitivity matrix of this compound-tissue (of a mouse replicate) with respect to all free fluxes
          x <- l.deriv_EMU.sources [ names(l.deriv_EMU.sources) %>% str_detect( w )] %>% as.data.frame() %>% rename_all(~ str_extract(., "v[0-9]+$")) # blood glucose   M+0 to M+6 derivatives to free fluxes
          # add row names
          rownames(x) <- str_c(tracer.k, "-", l.Compound.tissue.seq.m.id[[tracer.k]][i], " #", seq(nrow(x))-1) # here use # for number of 13C tags
          
          s <- rbind(s, x) # bind all different compound-tissue|mouse derivatives
        }
        
        # checkpoint:
        if (nrow(s) != length(IDV.simulated.k)) stop("Length of simulated labeling vector should be the same as the row number of the sensitivity matrix, for any given tracer.")
        
        # compile IDV derivatives from all tracers
        dX.target_allTracers <- rbind(dX.target_allTracers, s)
        
      }
      if (PRINT == T) cat("Looped through all tracers")
      
      
      
      # optional: print the cost function dynamically
      cost.i          <- sum( ( IDV.simulated_allTracers - IDV.obs)             ^2)  # not weighted
      cost.i.weighted <- sum( ((IDV.simulated_allTracers - IDV.obs)/IDV_SD.obs) ^2)  # weighted by standard deviation
      
      
      # keep a record of the history of the cost
      if (optim.step == 1) cost.history <- c()
      cost.history <- c(cost.history, cost.i)
      
      
      
      
      
      
      # # plot : iteration converging, within each binary search
      # cost shown in log scale
      # if (optim.step == 1) {
      #   plot(x = 1,
      #        y = cost.i %>% log(),
      #        xlim = c(1, 10),
      #        ylim = c(optimal_cost %>% log10() - .1, cost.i %>% log10() + 1),
      #        pch = 19, ylab = "log10(cost)", cex.lab = 1.4, cex.axis = 1.4,
      #        main = str_c("working on: stage - ", stage,  ", binary index - ", binary.index))
      # 
      #   abline(h = optimal_cost %>% log10(), col = 'green', lty = 2, lwd = 2)
      # }
      # if (optim.step %% 1 == 0) { # mark the objective function every x iteration steps
      #   points(x = optim.step, y = cost.i %>% log10(), pch = 19)
      # }
      
      
      
      # Enter the quadratic programming sub-problem -----
      if (PRINT == T) print("Moving on to stack the sensitive matrices, and solve the quadratic programming problem.\n\n")
      
      dX.target_allTracers <- as.matrix(dX.target_allTracers) # convert sensitive tibble to matrix
      
      # not weighted
      J <- t(dX.target_allTracers) %*% (IDV.simulated_allTracers - IDV.obs)
      H <- t(dX.target_allTracers) %*% (dX.target_allTracers)
      
      # # # # weighted by standard deviation
      # J <- t(dX.target_allTracers) %*% diag(1/(IDV_SD.obs^2)) %*% (IDV.simulated_allTracers - IDV.obs)
      # H <- t(dX.target_allTracers) %*% diag(1/(IDV_SD.obs^2)) %*% (dX.target_allTracers)
      
      
      
      # # Ensure H is positive definite
      # # method 1
      eg <- eigen(H)
      eigenValues <- eg$values; eigenValues
      
      tolerance <- max( max(eigenValues) * .Machine$double.eps, .Machine$double.eps) 
      
      if ( any((eigenValues < tolerance) == T) ){ # if any eigenvalues are close to zero, add a small number
        eigenValues[eigenValues < tolerance] <- tolerance
        H <- eg$vectors %*% diag(eigenValues) %*% t(eg$vectors) # reconstruct H
      }
      
      # method 2
      # H <- H + diag(1, nrow = nrow(H)) * .Machine$double.eps * 100
      H <- H + diag(1, nrow = nrow(H)) * (10^-10)
      
      
      
      
      # first set up the constraint of u after update (then convert to constraint of delta u later)
      # first set up the constraint of u after update (then convert to constraint of delta u later)
      Amat0 <- m.freeFlux_to_fullSetFlux
      myMin <- .01 # non-negativity lowest bound; use a small non-zero value, instead of abs. zero, to avoid that isotopologues fail to sum up to 1
      bo <- rep(myMin, nrow(Amat0)) # all fluxes >= 0 # instead of setting zero, set it to a small positive number; otherwise certain cases isotopologue does not sum to 1
      
      if (nrow(Amat0) != length(bo)) stop("Dimension wrong.")
      
      
      
      
      
      # For the given tissue, for reactions involving CO2, add row index for reactions producing CO2, and minus row index for reactions consuming CO2 
      func.add_rows_by_index <- function(myIndex){ # input is vector of reaction (enzyme) index involving CO2 production (positive index) and consumption (negative index)
        
        result <- rep(0, ncol(Amat0))
        
        for (index in myIndex) {
          if (index > 0) {
            # note!!! Use 'Amat0' the original Stoich Null space free-to-full-flux mapping matrix
            # which is before the addition of constraint rows
            result <- result + Amat0[index, ] 
          } else if (index < 0) {
            result <- result - Amat0[-index, ]
          }
        }
        return(result)
      }
      
      
      # CO2 total flux to sink = 1800
      # stacking on TOP!!! of 'Amat0' for meq = 1
      Amat <- rbind(  func.add_rows_by_index(i.CO2_sink), Amat0 ) 
      Amat <- rbind( -func.add_rows_by_index(i.CO2_sink), Amat )
      bo   <- c (totalCO2-1, bo)
      bo   <- c(-totalCO2-1, bo)
      
      
      if (nrow(Amat) != length(bo)) stop("Dimension wrong.")
      
      
      # The reaction whose C.I. we seek to determine is fixed to a given level, e.g. increase by 10% of the optimal value
      # For equality constraints, stack on TOP!!! of 'Amat0' for meq = 1
      Amat <- rbind(  Amat0[flux.C.I_index1, ] - Amat0[flux.C.I_index2, ] , Amat)
      Amat <- rbind(-(Amat0[flux.C.I_index1, ] - Amat0[flux.C.I_index2, ]), Amat)
      bo   <- c(  optimal_allFlux[flux.C.I_index1] - optimal_allFlux[flux.C.I_index2] + h   , bo)
      bo   <- c(-(optimal_allFlux[flux.C.I_index1] - optimal_allFlux[flux.C.I_index2] + h)-1, bo)
      
      
      # when working towards the lower bound, approaching LB -> 0, 
      # if the tested value is smaller than the non-negativity constraint, stop the binary search
      if  ((optimal_allFlux[flux.C.I_index1] - optimal_allFlux[flux.C.I_index2] + h) < myMin) {
        # 'break' to break out the iterating converging
        # and with binary index set to +Inf binary index, stop the binary search
        binary.index <-  +Inf # break out the allowed binary search total steps
        break 
      }
      
      
      # find the index by reaction name
      func.findIndex <- function(myR = "Pyr.Lvm+CO2.Blood->OAA.Lvd"){ # reaction
        ri <- reaction_data %>% filter(reactions == myR) %>% pull(enzyme)
        if (is_empty(ri)) stop("Did not find the index. Check spelling of reaction names")
        return(ri)
      }
      
      
      # Fluxes cannot be infinitely big (e.g., ≤ 10^8)
      # Amat <- rbind( Amat, -Amat0 )
      # bo   <- c(bo, rep(-10^28, nrow(Amat0)))
      
      
      # # PEPCK flux cannot be (e.g., 100) times bigger than citrate synthesis flux (due to ATP constraint)
      # # Lv
      # i.PEPCK.Lv  <- func.findIndex("Pyr.Lvm+CO2.Blood->OAA.Lv")
      # i.CS.Lv     <- func.findIndex("OAA.Lv+AcCoA.Lv->Cit.Lv")
      # Amat <- rbind( Amat, 100 * Amat0[i.CS.Lv, ] - Amat0[i.PEPCK.Lv, ] )
      # bo   <- c(bo, 0)
      # 
      # 
      # # pyruvate-lactate exchange flux: the relative flux order of magnitude follows the TCA fluxes 
      # i.Pyr_Lac.Lv  <- func.findIndex("Pyr.Lvc->Lac.Blood")
      # i.Pyr_Lac.M   <- func.findIndex("Pyr.Mc->Lac.Blood")
      # ## M > Lv
      # Amat <- rbind( Amat, Amat0[i.Pyr_Lac.M, ] - Amat0[i.Pyr_Lac.Lv, ])
      # bo <- c(bo, 0)
      # 
      # 
      # the alanine release flux pyr->ala roughly follows the tissue mass / TCA activity
      # i.Pyr_Ala.M  <- func.findIndex("Pyr.Mc->Ala.Blood")
      # i.Pyr_Ala.Lv <- func.findIndex("Pyr.Lvc->Ala.Blood")
      # ## M > Lv
      # Amat <- rbind( Amat, Amat0[i.Pyr_Ala.M, ] - Amat0[i.Pyr_Ala.Lv, ])
      # bo <- c(bo, 0)
      # 
      # 
      # # most PEP-derived pyruvate leaks out to circulation (Lv and K): Pyr->lac.blood >> PEP->Pyr; i.e., Lac.Blood input flux dominates the influx of pyruvate in tissues
      # # this hypothesis is relaxed in brain, where glucose->pyr->Ox is a major pathway
      # i.PEP_Pyr.M  <- func.findIndex("PEP.M->Pyr.Mc")
      # i.PEP_Pyr.Lv <- func.findIndex("PEP.Lv->Pyr.Lvc")
      # ## M
      # Amat <- rbind( Amat, Amat0[i.Pyr_Lac.M, ] - 1 * Amat0[i.PEP_Pyr.M, ])
      # bo <- c(bo, 10) # Pyr->lac.blood above some lower bound (e.g. above 100) at least if PEP->Pyr is small flux
      # ## Lv
      # Amat <- rbind( Amat, Amat0[i.Pyr_Lac.Lv, ] - 1 * Amat0[i.PEP_Pyr.Lv, ])
      # bo <- c(bo, 10) 
      # 
      # 
      # 
      # # OAA-MA-Suc exchange flux proceeds at reasonably fast flux, though not necessarily at complete equilibrium
      # ## Lv
      # i.OAA.Mal.Lv <- func.findIndex("OAA.Lv->Mal.Lv")
      # i.Mal.Suc.Lv <- func.findIndex("Mal.Lv->Suc.Lv")
      # Amat <- rbind( Amat, Amat0[i.OAA.Mal.Lv, ] - .1 * Amat0[i.CS.Lv, ])
      # bo <- c(bo, 0)
      # Amat <- rbind( Amat, Amat0[i.Mal.Suc.Lv, ] - .1 * Amat0[i.CS.Lv, ])
      # bo <- c(bo, 0)
      # ## M
      # i.OAA.Mal.M <- func.findIndex("OAA.M->Mal.M")
      # i.Mal.Suc.M <- func.findIndex("Mal.M->Suc.M")
      # i.CS.M <- func.findIndex("OAA.M+AcCoA.M->Cit.M")
      # Amat <- rbind( Amat, Amat0[i.OAA.Mal.M, ] - .1 * Amat0[i.CS.M, ])
      # bo <- c(bo, 0)
      # Amat <- rbind( Amat, Amat0[i.Mal.Suc.M, ] - .1 * Amat0[i.CS.M, ])
      # bo <- c(bo, 0)
      # 
      # 
      # # ignore muscle malic enzyme activity
      # # i.ME.M <- func.findIndex("Mal.M->Pyr.M+CO2.Blood")
      # # Amat <- rbind( Amat, -Amat0[i.ME.M, ])
      # # bo <- c(bo, -0.2)
      
      
      
      
      bo <- bo - Amat %*% u  # convert to constraint of delta v
      
      # If the QP has error (H is not positive definite, due to numeric instability), skip to the next for iteration
      qp <- NULL
      qp <- tryCatch({
        # Attempt to execute the problematic line
        quadprog::solve.QP(Dmat = 2*H, dvec = -2*J, Amat = t(Amat), bvec = bo, meq = 0) # CI-interested flux fixed, and total CO2 fixed to 1800 
      }, error = function(e) {
        # Handle the error and jump to the next iteration
        message(e$message)
        
        return(NULL)  # Return NULL or some default value
      })
      
      if ( is.null(qp) ) break
      
      # fumarate scrambling is not counted as the optimization program, as it's a dead fixed value
      u.step.update <- qp$solution # update the rest of the full set of fluxes
      
      u <- u + u.step.update # update the free fluxes # -------
      
      if (PRINT == T) print(paste("finish optimization step:", optim.step))
      
      
      
      
      
      # after a given number of iterations, check convergence:
      # the last several cost values should have a small variation (low relative percentage error)
      if (optim.step >= 5) {
        cost_latest5iterations <- cost.history[(optim.step - 2) : optim.step]
        cost_latest5iterations.error.pct <- sd(cost_latest5iterations) / mean(cost_latest5iterations) * 100
        
        if (cost_latest5iterations.error.pct < 2) {
          # dev.off()
          break 
        }
      }
      # if the optimization stops, keep a record of the current fixed flux value and the cost
      
      
      
      # track the parallel progress
      # update by saving the new current marker file
      # save(IDV.obs, file = paste0("parallelism/progress ", "MC-",  MC, "  iteration-", optim.step, "_R.Data"))
      
      # deletion optional: delete last iteration of the same MC loop
      # list.files("parallelism", pattern = paste0("progress ", "MC-",  MC, "  iteration-", optim.step-1, "_R.Data"), full.names = TRUE) %>% unlink() 
      # if (optim.step == n.opt.steps) list.files("parallelism", pattern = paste0("progress ", "MC-",  MC), full.names = TRUE) %>% unlink() 
    }                                   
    
    
    # F stats
    n.minus.p <- nrow(dX.target_allTracers) - ncol(dX.target_allTracers)
    F.critical.95 <- qf(p = .95,  df1 = 1, df2 = n.minus.p, lower.tail = T)
    cost.threshold.95 <- F.critical.95 / n.minus.p * optimal_cost + optimal_cost
    
    
    
    # dev.off()
    # # optional: dynamic plotting for flux - cost update (bottom half of the viz. device)
    # # y axis shown in linear scale (not log)
    # if (binary.index == 1 & stage == "UB") { # create the base plot showing the initial random flux-generated IDV
    # 
    #   plot(x = flux.C.I_value.optimal,
    #        y = optimal_cost,
    #        xlim = c(0, ifelse( round(flux.C.I_value.optimal)==0, 200, flux.C.I_value.optimal * 4)),
    #        ylim = c(optimal_cost, cost.threshold.95 * 1.3),
    #        col = "green3", pch = 19, ylab = "cost", xlab = "flux values", cex.lab = 1.4, cex.axis = 1.4)
    # 
    #   abline(h = optimal_cost,           col = 'green3', lty = 2, lwd = 2) # horizontal line - the optimal cost value
    #   abline(v = flux.C.I_value.optimal, col = 'green3', lty = 2, lwd = 2) # vertical line   - the optimal flux value
    #   abline(h = cost.threshold.95, col = 'black', lty = 2, lwd = 2)      # threshold
    # } else {
    #   replayPlot(p.CI)
    # }
    # points(x = flux.C.I_value.optimal + h, y = cost.i, col = 'black', pch = 19)
    # 
    # if (stage == "UB") abline(v = flux.C.I_value.optimal + h, col = "tomato")  # upper bound
    # if (stage == "LB") abline(v = flux.C.I_value.optimal + h, col = "skyblue") # lower bound
    # 
    # p.CI <- recordPlot()
    
    
    
    # cost at the current h adjustment 
    d.flux.cost.h <- tibble(flux.index = paste0(flux.C.I_index1, "-", flux.C.I_index2),
                            flux.value = flux.C.I_value.optimal + h,
                            cost = cost.i)
    print(d.flux.cost.h)
    
    # compile the result
    d.flux.cost <- bind_rows(d.flux.cost, d.flux.cost.h)
    
    # Update the h for the next binary search step
    binary.index <- binary.index + 1
    
    if (stage == "UB") {
      
      if (lock == F) {
        # if not higher than threshold, keep moving to the right with big steps
        if (cost.i < cost.threshold.95 & binary.index <= total.binary.steps) {
          h <- h + h0
          
        } else { 
          lock <- T # locking the UB within the defined interval 
          
          R = flux.C.I_value.optimal + h
          L = R - h0
          h <- (L + R)/2 - flux.C.I_value.optimal
          next
        }
      }
      
      if (lock == T) {  
        
        if (cost.i >= cost.threshold.95){
          # search the right side
          R <- flux.C.I_value.optimal + h
          L <- L # being the same keft boundary of last iteration
        } else {
          # search the left side
          L <- flux.C.I_value.optimal + h
          R <- R # being the same right boundary of last iteration
        }
        
        h <- (L + R)/2 - flux.C.I_value.optimal
      }
      
    }
    
    
    # Working on CI lower bound, immediately with lock = T; 
    # if flux optimal value is already zero then do not determine the lower boundary
    if (stage == "UB" & binary.index == total.binary.steps & round(flux.C.I_value.optimal) != 0 ) {
      stage <- "LB"
      binary.index <- 1 # restart the binary counting index
      h <- - h0 / 2
      L <- 0.01 # should be 0; but to avoid isotopologues do not sum up to 100%, set it to a small non-zero number (smaller than QP meq constraint)
      R <- flux.C.I_value.optimal
      next
    }
    
    if (stage == "LB"){
      if (cost.i >= cost.threshold.95) {
        # working on the right half side of the interval
        L <- flux.C.I_value.optimal + h
        R <- R # right boundary being the same as the last iteration
      } else {
        # working on the left half side of the interval
        R <- flux.C.I_value.optimal + h
        L <- L # left boundary being the same as the last iteration
      }
      h <- (L + R)/2 - flux.C.I_value.optimal
    }
    
  }                                     
  
  save(d.flux.cost, 
       n.minus.p,
       F.critical.95,
       cost.threshold.95,
       optimal_allFlux, # optimal from repeated run
       optimal_cost, 
       reaction_data,
       file = paste0("./CI/", flux.C.I_index1, "-", flux.C.I_index2, ".RData"))
  
} 




# Export necessary variables and functions to the cluster

# # Update the working nodes OPTION 1: export directly to working nodes
# time.Export2Nodes <- system.time(
#   clusterExport(cl, varlist = ls())
# )
# time.Export2Nodes



# Update the working nodes OPTION 2: save current ls() as a separate file, and read it directly into each working nodes
# this approach is much FASTER (e.g. x30 faster for full body multiorgan model, but only slightly faster for the simplified Lv-M model)
save(list = ls(), file = "temporary_workspace_for_workers.RData")

time.Export2Nodes <- system.time( 
  clusterEvalQ(cl, {
    load("temporary_workspace_for_workers.RData")
  }) 
)
time.Export2Nodes

file.remove("temporary_workspace_for_workers.RData")




# run parallel computation
parLapply(cl, f.all, func.calculate_confidence_interval)


time.End.CI <- Sys.time()

sink("log_CI.txt")
time.End.CI - time.Start.CI
sink()
beepr::beep(sound = "coin")


message("\n\nCompleted running confidence interval for net fluxes of reversible reactions!\n\n")

