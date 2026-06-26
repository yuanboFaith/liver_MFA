# This scripts performs iteration convergence and quadratic programming to find the optimal fluxes. 

# ---------
library(readxl)
library(matrixcalc)
library(Ryacas) # combine like terms
library(foreach)
library(parallel)
library(doParallel)
library(tidyverse)


# rm(list = ls()) # this line is not active, so that this file can be run and controlled in a loop of another file

# rstudioapi::getActiveDocumentContext()$path %>% dirname() %>% setwd(); getwd() # set path to the folder of the current active script



# # Import 13C tracing labeling data, and infusion rate
# load("../data/13C_tracing_labeling_data.RData")
# d.labeling.summary





# Import Stoic and EMU decomposition result from different tracing experiments
# 'state_fasted.or.refed' is a master control variable. It may be commented out here and specified in other files
state_fasted.or.refed <- "fasted"


if (state_fasted.or.refed == "fasted") {
  load("3_list.stoich.EMU_allTracers_fasted.RData") # refed-specific EMU decomposed result 
  d.13C.labeling <- d.13C.fasted
  totalCO2       <- 1800 
}


# if (state_fasted.or.refed == "refed") {
#   load("3_list.stoich.EMU_allTracers_refed.RData") # refed-specific EMU decomposed result 
#   d.13C.labeling <- d.13C.refed.hpAA  
#   totalCO2       <- 1800 * 1.2
# }


# cleaned labeling - infusion data (from the 'data' folder)
# !!! important to use the new data set (combining the AA tracing data in the refed dataset)
load(file = "../data/cleaned_labeling_data.RData")




# each state, fasted or refed, will generate its own .RData file containing the correction factors for coarse-grained modeling
# however, for both states, regarding their the intermediately outputs,
# e.g., the simulated IDV labeling, simulated EMUs labeling, and sensitivity matrix dX, 
# they will be saved in the same folder tracer_k_parallel

# clear the tracer_k_parallel folder that saves intermediary outputs of the QP iteration
folder_path <- "tracer_k_parallel"
rdata_files <- list.files(path = folder_path, pattern = "\\.RData$", full.names = TRUE)
file.remove(rdata_files) # remove 


# Import supplementary functions
load("1_Supplement_functions.RData")



time.Start.iteration <- Sys.time()



# tissues included in the EMU decomposition model
TissuesModeled <- target.EMU %>% str_extract("(?<=\\.)[:alpha:]{1,9}(?=\\_)") %>% unique()



# Select labeling data with the same scope of tissues included in the EMU decomposition model
d.13C.labeling <- d.13C.labeling %>% filter(tissue %in% TissuesModeled)

d.13C.labeling$tissue %>% unique()
d.13C.labeling$Compound %>% unique()


# make sure the modeled tissues are also included in the labeling dataset
if ( (sum(TissuesModeled %in% unique(d.13C.labeling$tissue)) != n_distinct(d.13C.labeling$tissue)) ) error("Modeled tissues are not included in the labeling dataset.")

# in the labeling dataset, include only compound-tissue that is included in the target EMU list of the model network
d.13C.labeling <- d.13C.labeling %>% filter(Compound.tissue.seq %in% target.EMU)
d.13C.labeling$Compound %>% unique()


# this line is used for both non-parallel and parallel computation;
# for non-parallel, it is compiled during the last iteration; for parallel, it's compiled from exported files after conclusion of iteration
# record the simulated labeling in all tracer infusion at the last step of iteration
# this is used to calculate the simulated total 13CO2 produced from each organ, to further update the coarse grained model
l.EMU.sim <- list()




# Create initial starting values of free fluxes, with constraints that all fluxes are positive
# Here we do not have any equality constraints.
GG <- rbind(m.freeFlux_to_fullSetFlux,  # lower bound (non-negativity constraint)
            -m.freeFlux_to_fullSetFlux) #  upper bound (e.g., <1000)

HH <- c(rep(0, nrow(m.freeFlux_to_fullSetFlux)),   # lower bound (non-negativity constraint)
        rep(-500, nrow(m.freeFlux_to_fullSetFlux))) #  upper bound (e.g., <1000)

# don't select the 1st row of sampled initial start; there are too many zeros in the first sampled row
message("\n\nGenerating initial random fluxes...")
mat.u.initial <- limSolve::xsample(
  G = GG, H = HH, # lower (non-negativity) and upper (non-infinity) bound
  iter = 50)$X
mat.u.initial <- mat.u.initial[-1, ] # remove the 1st all-zero row


tracers <- names(l.stoich.EMU_allTracers)



l.IDV.obs <- list() # measured IDV; note that the tracer is in order of 'l.stoich.EMU_allTracers'; each tracer data is one element
IDV.obs   <- c()    # the vector form of l.IDV.obs
IDV_SD.obs  <- c()    # the vector form of SD of labeling
d.u.iterations <- tibble(.rows = 0)  # record the convergence of fluxes 
cost.iterations <- c() # record the convergence of the cost
l.Compound.tissue.seq      <- list() # compound-tissue to extract repeatedly for difference mice to construct the sensitivity matrix; for extraction elements of list
l.Compound.tissue.seq.m.id <- list() # same list as ablove, but adding mouse id; for naming purpose

l.13C <- list() # a list of labeling data

for (b in seq(tracers)){ # b the index
  # print(b)
  # check
  # b=1
  # b=6
  a <- tracers[b]        # a the tracer name
  # print(a)
  
  # metabolite measured labeling
  d.13C.labeling$Compound %>% unique()
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


n.opt.steps <- 25 # iterate optimization steps
PRINT <- F # if True, progress at major calculation steps is printed




# Define a function to process each tracer
func.process_tracer <- function(tracer.k) {
  
  message("iteration ", optim.step,  ": working on ", tracer.k)
  
  # tracer.k= "13ChpAAb"
  
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
        Y [iii, ] <- l.EMU.sources [Y.names[iii]] %>% as_tibble() %>% as.matrix() %>% t()
        
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
  
  
  # # This line is for non-parallel computation
  # # compile simulated IDV from all tracers 
  # IDV.simulated_allTracers <- c(IDV.simulated_allTracers, IDV.simulated.k)
  # for parallel computation, instead output the result to file
  save(IDV.simulated.k, file = paste0("./tracer_k_parallel/IDV.simulated_", tracer.k, ".RData"))
  
  
  
  
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
  
  
  # # # This line is for non-parallel computation
  # # compile IDV derivatives from all tracers
  # dX.target_allTracers <- rbind(dX.target_allTracers, s)
  # for parallel computation, instead output the result to file
  save(s, file = paste0("./tracer_k_parallel/dX_", tracer.k, ".RData"))
  
  
  # # # # This line is for non-parallel computation
  # # for the last iteration, record the simulated labeling
  # if (optim.step == n.opt.steps) {
  #   l.EMU.sim[[tracer.k]] <- l.EMU.sources
  # }
  
  
  # Also save the simulated EMU labeling 
  # Only the last step iteration will be used for constructing the correctino matrix
  save(l.EMU.sources, file = paste0("./tracer_k_parallel/l.EMU.sources_", tracer.k, ".RData"))
  
}


# # test
# func.process_tracer(tracer.k = "13CGlc")
# func.process_tracer(tracer.k = "13CLac")
# func.process_tracer(tracer.k = "13CGlycerol")
# func.process_tracer(tracer.k = "13CAla")
# func.process_tracer(tracer.k = "13CGln")
# func.process_tracer(tracer.k = "13CPalm")
# func.process_tracer(tracer.k = "13COle")
# func.process_tracer(tracer.k = "13CLino")
# func.process_tracer(tracer.k = "13CHB")
# func.process_tracer(tracer.k = "13ChpAAb")
# func.process_tracer(tracer.k = "13ChpAAc")
# func.process_tracer(tracer.k = "13ChpAAd")
# func.process_tracer(tracer.k = "13ChpAAe") # 
# func.process_tracer(tracer.k = "13ChpAAf") #
# func.process_tracer(tracer.k = "13ChpAAj") #
# func.process_tracer(tracer.k = "13ChpAAk") #
# func.process_tracer(tracer.k = "13ChpAAl") #




# Create a cluster
cl <- makeCluster(5, outfile = "" )  

# Make libraries accessible to each working node
clusterEvalQ(cl, {
  library(parallel)  
  library(readxl)
  library(Ryacas)
  library(tidyverse)
})


# # Export necessary variables and functions to the cluster
# clusterExport(cl, varlist = ls())




for (optim.step in 1:n.opt.steps) {
  
  
  cat("\n")
  message("------------------iteration ", optim.step, " ------------------------------------------------------")
  # check
  # optim.step <- 1
  if (optim.step == 1){
    
    # option (1) initiate random flux, s.t. all fluxes ≥ 0
    
    set.seed(as.integer(Sys.time())) # ensure each run generates random initial fluxes
    
    u.initial <- mat.u.initial[ sample(1:nrow(mat.u.initial), size = 1), ] * sample(x=c(.1, .5, 1, 5, 10, 50, 100, 500, 1000), size = 1)
    u.initial %>% plot(main = "initiation values of free fluxes")
    
    c(m.freeFlux_to_fullSetFlux %*% u.initial) %>% round()
    
    names(u.initial) <- paste0("v", names(u.initial))
    
    # testing
    # u.initial=c(2359032, 4718065, 7231253, 4937636, 0, 7301842, 3013842, 1590596, 6447854, 780148, 6726560)
    # names(u.initial)=paste0("v", index.freeFluxes)
    
    u <- u.initial # free fluxes
    
    # Replace a single flux notation with tangible flux value; 
    # define this appears important in parallel computation; otherwise u can't be found
    func.FillwithFluxValue <- function(x, whichFlux = "v10"){
      str_replace_all(x, pattern = paste0("\\b", whichFlux, "\\b"), replacement = u[whichFlux] %>% as.character() )
    }
    
    message("exporting global environment variables to working nodes...")
    
    
    # # Update the working nodes OPTION 1: export directly to working nodes
    # time.Export2Nodes <- system.time(
    #   clusterExport(cl, varlist = ls())
    # )
    
    # Update the working nodes OPTION 2: save current ls() as a separate file, and read it directly into each working nodes
    # this approach is x30 FASTER
    save(list = ls(), file = "temporary_workspace_for_workers.RData")
    
    time.Export2Nodes <- system.time( 
      clusterEvalQ(cl, {
        load("temporary_workspace_for_workers.RData")
      }) 
    )
    
    message("Finished global objects export to working nodes.")
    print(time.Export2Nodes)
    
    message("\nCalculating the first converging step...")
    
  }
  
  IDV.simulated_allTracers <- NULL # collecting the               simulated glucose labeling from all tracer experiments
  dX.target_allTracers     <- NULL # collecting the derivative of simulated glucose labeling from all tracer experiments; dX is same as derivative of IDV
  
  
  # check
  # tracer.k="13CGlc"
  # tracer.k="13CLac"
  # tracer.k="13CGlycerol"
  
  ### Below is the line looping through tracers
  # for (tracer.k in tracers){}
  # In parallel computation, each loop is defined as function 'func.process_tracer', and put above outside of the iteration loop
  
  
  
  # Execute in parallel
  parLapply(cl, tracers, func.process_tracer)
  
  # func.process_tracer(tracer.k = "13CGlc")
  # func.process_tracer(tracer.k = "13CGlycerol")
  
  
  # # Instead of compiling during non-parallel run using the following code, compile from the exported files
  # # compile simulated IDV from all tracers
  # IDV.simulated_allTracers <- c(IDV.simulated_allTracers, IDV.simulated.k)
  # # compile IDV derivatives from all tracers
  # dX.target_allTracers <- rbind(dX.target_allTracers, s)
  
  
  
  
  # Compile data from the exported files of the simulated IDV
  files.IDV.simulated <- list.files("./tracer_k_parallel",  full.names = TRUE, pattern = "^IDV\\.simulated_")
  # !!! Reorder the files according to 'tracer' ; files are shown by default in alphabetical order
  files.IDV.simulated <- files.IDV.simulated[match(tracers, str_extract(files.IDV.simulated, "13C[:alpha:]{1,1000}"))]
  
  for (file.i in files.IDV.simulated) {
    IDV.simulated.k <- load(file.i) %>% get()
    IDV.simulated_allTracers <- c(IDV.simulated_allTracers, IDV.simulated.k)
  }
  
  
  
  # Compile data from the exported files of the sensitivity matrix
  files.dX <- list.files("./tracer_k_parallel",  full.names = TRUE, pattern = "^dX_")
  # !!! Reorder the files according to 'tracer' ; files are shown by default in alphabetical order
  files.dX <- files.dX[match(tracers, str_extract(files.dX, "13C[:alpha:]{1,1000}"))]
  
  for (file.j in files.dX) {
    dX.k <- load(file.j) %>% get()
    dX.target_allTracers <- rbind(dX.target_allTracers, dX.k)
  }
  
  
  
  if (PRINT == T) cat("Looped through all tracers")
  
  
  # optional: print the cost function dynamically
  cost.i          <- sum( ( IDV.simulated_allTracers - IDV.obs)             ^2)  # not weighted
  cost.i.weighted <- sum( ((IDV.simulated_allTracers - IDV.obs)/IDV_SD.obs) ^2)  # weighted by standard deviation
  cost.iterations <- c(cost.iterations, cost.i)
  
  if (optim.step == 1) { # create the base plot showing the initial random flux-generated IDV
    plot(x = 1, 
         y = log10(cost.i),
         xlim = c(1,  n.opt.steps), 
         ylim = c(-1, log10(cost.i) + 1),
         pch = 19, ylab = "log10(cost)", cex.lab = 1.4, cex.axis = 1.4)
    
    abline(h =    0, col = 'red4', lty = 2, lwd = 1)
    abline(h =  -.2, col = 'red2', lty = 2, lwd = 1)
    abline(h =  -.4, col = 'orange', lty = 2, lwd = 1)
    abline(h =  -.6, col = 'skyblue', lty = 2, lwd = 1)
    abline(h =  -.8, col = 'green4', lty = 2, lwd = 1)
    abline(h =   -1, col = 'green', lty = 2, lwd = 1)
    
  }
  if (optim.step %% 1 == 0) { # mark the objective function every five iteration steps
    points(x = optim.step, y = log10(cost.i), col = 'black', pch = 19)
  }
  
  message("cost is: ", cost.i)
  
  
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

  tolerance <- max( max(eigenValues) * .Machine$double.eps, .Machine$double.eps) * 100

  if ( any((eigenValues < tolerance) == T) ){ # if any eigenvalues are close to zero, add a small number
    eigenValues[eigenValues < tolerance] <- tolerance
    H <- eg$vectors %*% diag(eigenValues) %*% t(eg$vectors) # reconstruct H
  }

  # # method 2
  # # H <- H + diag(1, nrow = nrow(H)) * .Machine$double.eps * 100
  # H <- H + diag(1, nrow = nrow(H)) * (10^-10)
  
  
  
  # first set up the constraint of u after update (then convert to constraint of delta u later)
  Amat0 <- m.freeFlux_to_fullSetFlux
  bo <- rep(.01, nrow(Amat0)) # all fluxes >= 0 # instead of setting zero, set it to a small positive number; otherwise certain cases isotopelogue does not sum to 1
  
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
  
  # tissue TCA constraints
  # Stack on the TOP !! of Amat (for equality constraints)
  # Stack on the BOTTOM !! of Amat  for inequality constraints
  # Alternatively, considering the floating error, use constraint as inequality: constraint - 1 ≤ flux ≤ constraint + 1
  
  
  # CO2 total flux to sink = 1800
  # stacking on TOP!!! of 'Amat0' for meq = 1
  Amat <- rbind(  func.add_rows_by_index(i.CO2_sink), Amat0 ) 
  Amat <- rbind( -func.add_rows_by_index(i.CO2_sink), Amat )
  bo   <- c (totalCO2-1, bo)
  bo   <- c(-totalCO2-1, bo)
  
  if (nrow(Amat) != length(bo)) stop("Dimension wrong.")
  
  
  
  
  
  
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
  # 
  # 
  # # the alanine release flux pyr->ala roughly follows the tissue mass / TCA activity
  # i.Pyr_Ala.M  <- func.findIndex("Pyr.Mc->Ala.Blood")
  # i.Pyr_Ala.Lv <- func.findIndex("Pyr.Lvc->Ala.Blood")
  # ## M > Lv
  # Amat <- rbind( Amat, Amat0[i.Pyr_Ala.M, ] - Amat0[i.Pyr_Ala.Lv, ])
  # bo <- c(bo, 0)
  # 
  # 
  # 
  # 
  # # most PEP-derived pyruvate leaks out to circulation (Lv and K): Pyr->lac.blood >> PEP->Pyr; i.e., Lac.Blood input flux dominates the influx of pyruvate in tissues
  # # this hypothesis is relaxed in brain, where glucose->pyr->Ox is a major pathway
  # ## M
  # i.PEP_Pyr.M  <- func.findIndex("PEP.M->Pyr.Mc")
  # Amat <- rbind( Amat, Amat0[i.Pyr_Lac.M, ] - 1 * Amat0[i.PEP_Pyr.M, ])
  # bo <- c(bo, 10) # Pyr->lac.blood above some lower bound (e.g. above 100) at least if PEP->Pyr is small flux
  # # ## Lv
  # # i.PEP_Pyr.Lv <- func.findIndex("PEP.Lv->Pyr.Lvc")
  # # Amat <- rbind( Amat, Amat0[i.Pyr_Lac.Lv, ] - 1 * Amat0[i.PEP_Pyr.Lv, ])
  # # bo <- c(bo, 10) 
  # 
  # 
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
  
  
  
  
  bo <- bo - Amat %*% u  # convert to constraint of delta v
  qp <- quadprog::solve.QP(Dmat = 2*H, dvec = -2*J, Amat = t(Amat), bvec = bo, meq = 0)
  
  u.step.update <- qp$solution # update the rest of the full set of fluxes
  
  
  # osqp::osqpSettings(eps_abs = .Machine$double.eps, eps_rel = .Machine$double.eps,
  #                    eps_prim_inf = .Machine$double.eps,# primal infeasibility tol
  #                    eps_dual_inf = .Machine$double.eps,# dual infeasibility tol
  #                    max_iter  = 10000000) # allow more iterations
  # qp <- osqp::solve_osqp(P = 2*H, q = 2*J, A = Amat, l = bo)   # relative tolerance)
  # u.step.update <- qp$x # update the rest of the full set of fluxes
  # # 
  
  # # If the QP has error (H is not positive definite, due to numeric instability), skip to the next for iteration
  # qp <- NULL
  # qp <- tryCatch({
  #   # Attempt to execute the problematic line
  #   quadprog::solve.QP(Dmat = 2*H, dvec = -2*J, Amat = t(Amat), bvec = bo, meq = 0) 
  # }, error = function(e) {
  #   # Handle the error and jump to the next iteration
  #   message(e$message)
  #   
  #   return(NULL)  # Return NULL or some default value
  # })
  # 
  # if ( is.null(qp) ) break
  
  
  # fumarate scrambling is not counted as the optimization program, as it's a dead fixed value
  
  
  
  u <- u + u.step.update # update the free fluxes # -------
  
  
  # record the convergence of fluxes
  d.u.iterations.i <- tibble(ite = optim.step, index = names(u), u = u)
  d.u.iterations <- bind_rows(d.u.iterations, d.u.iterations.i)
  
  
  
  
  # after a given number of iterations, check convergence:
  # the last several cost values should have a small variation (low relative percentage error)
  if (optim.step >= 15) {
    cost_latest5iterations <- cost.iterations[(optim.step - 6) : optim.step]
    cost_latest5iterations.error.pct <- sd(cost_latest5iterations) / mean(cost_latest5iterations) * 100
    
    if (cost_latest5iterations.error.pct < .1) {
      # dev.off()
      break 
    }
  }
  
  
  
  
  message("updating u to nodes")
  # clusterExport(cl, varlist = ls())
  clusterExport(cl, varlist = list("u"))
  
  
  # count number of negative free fluxes
  neg.fluxes <- sum(m.freeFlux_to_fullSetFlux %*% u < 0.001)
  if (neg.fluxes >= 1) stop("there is negative flux values calcualted")
  
  
  if (PRINT == T) print(paste("finish optimization step:", optim.step))
  
  
}


stopCluster(cl)

# -----Post simulation analysis ----------------------------------------------------------------
message("Completed MFA iteration!")

dX.target_allTracers %>% dim()
qr(dX.target_allTracers)$rank




# Compare the simulated vs observed labeling
d.obs.simu <-  tibble(obs = IDV.obs, 
                      sim = IDV.simulated_allTracers,
                      IDV = names(IDV.obs)) %>% 
  separate(IDV, into = c("tracer", "metabolite", "fullSeq", "M+#"), sep = "[_-]")  %>% 
  separate(`M+#`, into = c("label", "mouse.id"), sep = "\\|") # %>% 
# filter(tracer == "13CGln") 
# mutate(tracer = factor(tracer, levels = tracers),
#        metabolite = factor(metabolite, levels = c("Glucose", "Lactate", "Alanine", "Glycerol", "CO2")))


# place tracers in order
ordered.tracers = paste0(
  "13C", c("Glc", "Lac", "Ala", "Gln", "Glycerol", "HB", "Palm", "Ole", "Lino", 
           paste0("hpAA", c("b", "c", "d", "e", "f", "j", "k", "l"))))

d.obs.simu <- d.obs.simu %>% mutate(tracer = factor(tracer, levels = ordered.tracers))


# sim vs. obs plot1
plt.obs.sim <- d.obs.simu %>% 
  # mutate(metabolite = str_replace(metabolite, "\\.Blood", "\\.Bld")) %>% 
  ggplot(aes(x = obs, y = sim, color = metabolite)) +
  geom_text(data = d.obs.simu %>%
              group_by(tracer, metabolite, label) %>%
              summarise(sim = mean(sim),
                        obs = mean(obs)),
            aes(color = label, label = label %>% str_remove("#")),
            size = 3, fontface = "bold") +
  scale_x_continuous(transform = "sqrt") +
  scale_y_continuous(transform = "sqrt") +
  facet_grid(metabolite ~ tracer) +
  geom_abline(slope = 1, intercept = 0, color = "grey") +
  coord_cartesian(xlim = c(0, .3), ylim = c(0, .3)) +
  # coord_cartesian(xlim = c(0.8, 1), ylim = c(0.8, 1)) +
  theme_bw(base_size = 13) +
  theme(axis.text.x = element_text(angle = 40, hjust = 1),
        legend.position = "none",
        panel.spacing = unit(1, "pt"),
        panel.grid = element_line(size = .3)) 
  
plt.obs.sim

# ggsave(filename = "./plots/sim vs obs.pdf", width = 5, height = 24)



# sim vs. obs plot2
p <- d.obs.simu %>% 
  mutate(metabolite = str_replace(metabolite, "\\.Blood", "\\.Bld")) %>% 
  ggplot(aes(x = obs, y = sim, color = tracer)) +
  # geom_text(aes(label = label)) +
  geom_point(size = 3, shape = 21, stroke = .3) +
  # geom_text(aes(label = label), size = 3 , fontface = "bold") +
  scale_x_continuous(transform = "sqrt") +
  scale_y_continuous(transform = "sqrt") +
  # facet_grid(metabolite ~ tracer) +
  # facet_wrap(~tracer) +
  geom_abline(slope = 1, intercept = 0, color = "black", linewidth = .3) +
  theme_bw(base_size = 13) +
  theme(axis.text.x = element_text(angle = 40, hjust = 1),
        # legend.position = "none",
        panel.spacing = unit(1, "pt"),
        panel.grid = element_line(size = .3)) +
  guides(color = guide_legend(override.aes = list(stroke = 1, shape = 19)))  +
  
  geom_point(data = d.obs.simu %>%
               group_by(tracer, metabolite, label) %>%
               summarise(sim = mean(sim),
                         obs = mean(obs)),
             aes(fill = tracer), color = "black", fill = "black", size = .5,
             shape = 21)

p02 <- p + coord_equal(xlim = c(0, .2), ylim = c(0, .2));   p02
# ggsave(filename = "./plots/sim vs obs LOW range.pdf", width = 5, height = 5)

p81 <- p + coord_equal(xlim = c(0.9, 1), ylim = c(0.9, 1)); p81
# ggsave(filename = "./plots/sim vs obs HIGH range.pdf", width = 5, height = 5)




# sim vs. obs plot4
p <- d.obs.simu %>% 
  mutate(metabolite = str_replace(metabolite, "\\.Blood", "\\.Bld")) %>% 
  ggplot(aes(x = obs, y = sim, color = metabolite)) +
  scale_x_continuous(transform = "sqrt") +
  scale_y_continuous(transform = "sqrt") +
  facet_wrap(~tracer, nrow = 2) +
  geom_abline(slope = 1, intercept = 0, color = "black", linewidth = .3) +
  theme_bw(base_size = 13) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        legend.position = "bottom",
        panel.spacing = unit(5, "pt"),
        panel.grid = element_line(size = .3)) +
  geom_text(data = d.obs.simu %>%
              group_by(tracer, metabolite, label) %>%
              summarise(sim = mean(sim),
                        obs = mean(obs)),
            aes(color = metabolite, label = label %>% str_remove("#")),
            size = 3) +
  guides(color = guide_legend(override.aes = list(size = 5, fontface = "bold"))) +
  labs(color = NULL)

p02.facet <- p + coord_equal(xlim = c( 0, .2), ylim = c(0, .2));   p02.facet
# ggsave(filename = "./plots/sim vs obs LOW range faceted.pdf", width = 6, height = 6)

p81.facet <- p + coord_equal(xlim = c(.8, 1),  ylim = c(.8, 1));   p81.facet
# ggsave(filename = "./plots/sim vs obs HIGH range faceted.pdf", width = 6, height = 6)





# plot cost converging profile
d.cost.iterations <-  tibble(cost = cost.iterations) %>% mutate(iteration = 1:nrow(.)) 

p.cost <- d.cost.iterations %>% 
  ggplot(aes(x = iteration, y = cost)) + 
  geom_line(color = "firebrick", linewidth = 1) + 
  ggrepel::geom_text_repel(aes(label = round(cost, 4)), min.segment.length = 0, color = "turquoise4") +
  scale_y_log10(breaks = c(seq(.1, 1, .1), 1:10, seq(10, 100, 10)), expand = expansion(add = 0.1)) + 
  annotation_logticks(sides = "l") +
  theme_classic(base_size = 15) + 
  labs(y = "cost") + 
  scale_x_continuous(breaks = 1:30) +
  annotate(geom = "label", x = 9, 
           y = max(d.cost.iterations$cost) * .5, 
           label = paste("cost =", slice_tail(d.cost.iterations, n = 1)$cost %>% round(6)),
           size = 5, fontface = "bold", fill = "snow2")
p.cost



# Calculated fluxes
d.v.final <- m.freeFlux_to_fullSetFlux %*% as.matrix(u) %>% # calculate all fluxes
  as_tibble() %>% 
  rename(flux = "V1") %>% 
  
  # add reaction name
  mutate(Flux.index = rownames(m.freeFlux_to_fullSetFlux) %>% as.integer()) %>% 
  left_join(reaction_data %>% rename(Flux.index = enzyme),  # 'reaction_data' originally loaded from the 'Stoic_EU_calculator.RData' -> 'decomposition_parallel.RData'
            by = "Flux.index") %>% 
  
  # arrange fluxes in order
  arrange(Flux.index) %>% 
  mutate(Flux.index = factor(Flux.index, levels = unique(.$Flux.index))) %>% 
  mutate(flux = round(flux, 1))

# remove fast equillibrium fluxes
d.v.final2 <- d.v.final %>% 
  # filter(! reactions %>% str_detect("Mal")) %>% 
  filter(! reactions %>% str_detect("Suc2")) 


# # plot simulated result
# d.v.final2 %>%
#   ggplot(aes(x = Flux.index, y = flux)) +
#   # geom_point() +
#   geom_text(aes(label = flux %>% round())) +
#   scale_y_continuous(n.breaks = 8, expand = expansion(mult = c(0, .1))) +
#   theme_light() +
#   theme(panel.grid.minor.x = element_blank(),
#         legend.position = "none",
#         axis.text.y = element_text(size = 7),
#         axis.title = element_text(size = 12),
#         plot.margin = margin(l = 180)) +
#   geom_hline(yintercept = 150, linetype = "dashed", size = .2) +
#   geom_text(data = d.v.final2 %>% select(Flux.index, reactions) %>% distinct(),
#             aes(y = -400, label = reactions),
#             hjust = 1, size = 4, color = "red3") +
#   coord_flip(ylim = c(0, 2800), clip = "off") +
#   labs(x = NULL)






# net CO2 output
d.flux.CO2 <- d.CO2 %>% rename(Flux.index  = "enzyme") %>% mutate(Flux.index = factor(Flux.index)) %>% 
  select(Flux.index, tissue, coef) %>% 
  left_join(d.v.final, by = "Flux.index") %>% 
  mutate(CO2.produce.consume = coef * flux) # if > 0, produce CO2; if < 0, consume CO2

d.flux.CO2$CO2.produce.consume %>% sum()

d.flux.CO2.summary <- d.flux.CO2 %>% group_by(tissue) %>% 
  summarise(netCO2 = sum(CO2.produce.consume)) %>% 
  mutate(pct = netCO2 / sum(netCO2))

p.lv.M <- d.flux.CO2.summary %>% 
  mutate(tissue = fct_reorder(tissue, pct)) %>% 
  ggplot(aes(x = 1, y = netCO2, fill = tissue)) +
  geom_col(color = "black", alpha = .7)  +
  geom_text(aes(label = paste(tissue, round(pct * 100, 1), "%")),
            position = position_stack(vjust = .5), angle = 90, size = 5, fontface = "bold") +
  scale_y_continuous(breaks = seq(0, 4000, 100),
                     expand = expansion(mult = c(0, 0.05))) +
  scale_fill_brewer(palette = "Set2") +
  # coord_polar(theta = "y", direction = -1)
  theme_void() +
  theme(axis.text.y = element_text(), 
        legend.position = "none") +
  scale_x_continuous(expand = expansion(add = 0))
p.lv.M



# plot residuals
d.obs.simu %>% mutate(residuals = sim - obs) %>%
  ggplot(aes(x = residuals, fill = metabolite)) +
  geom_vline(xintercept = 0, linewidth = .1) +
  geom_histogram(color = "black", position = "dodge", linewidth = 0) +
  # facet_grid(metabolite~tracer) +
  theme_bw() +
  theme(# legend.position = "none",
    panel.spacing = unit(0, units = "pt"),
    panel.grid = element_blank()) +
  scale_y_continuous(expand = expansion(mult = 0))



# plot convergence of fluxes
d.u.iterations %>% 
  # filter(u < 1000) %>% 
  ggplot(aes(x = ite, y = u, color = index)) +
  geom_point(size = 1) +
  geom_line() +
  geom_text(data = d.u.iterations %>% filter(ite == max(ite)),
            aes(x = ite + 1.5, label = index), size = 3) +
  theme_minimal() +
  theme(legend.position = "none") +
  coord_cartesian(ylim = c(0, 300)) +
  labs(x = "number of interations", y = "free fluxes (nmol/min/g BW)")

# ggsave("./plots/flux convergence.pdf", width = 6, height = 6)



d.u.iterations %>% filter(ite > 15) %>% 
  group_by(index) %>% 
  summarise(u_sd.percent = sd(u)/mean(u) * 100) %>% 
  filter(u_sd.percent > 100) %>% 
  arrange(index)


dX.target_allTracers %>% dim()
dX.target_allTracers %>% qr() %>% .$rank

# l.EMU.sources[names(l.EMU.sources) %>% str_detect("Pyr_2")]
# l.EMU.sources[names(l.EMU.sources) %>% str_detect("Pyr.M_2")]






# calculate all fluxes from the free flux
# Get the final optimized flux to determine confidence interval (in a separate script)
optimal_allFlux   <- m.freeFlux_to_fullSetFlux %*% as.matrix(u)
optimal_free_flux <- u
optimal_cost      <- cost.i


# save.image(file = "optimal_solution.RData")


# if (state_fasted.or.refed == "fasted") {
#   save(optimal_allFlux, optimal_free_flux, optimal_cost, file = "4_optimal_solution_fasted.RData")
# }
# if (state_fasted.or.refed == "refed") {
#   save(optimal_allFlux, optimal_free_flux, optimal_cost, file = "4_optimal_solution_refed.RData")
# }

# The small chunk outputing the best fit is eventually replaced by the best fit from repeated MFA runs in script #7_1 and #7_2



# Sys.sleep(5) # Figured that you need to pause a few seconds to the code below to work...but don't know why!
# message("Calculate the correction factors for 13C and 12C...")
# 
# # combine the simulated EMU labeling from the last step (using the exported files)
# 
# # This line is for non-parallel computation
# # l.EMU.sim[[tracer.k]] <- l.EMU.sources
# 
# Compile data from the exported files of the simulated IDV
files.EMU <- list.files("./tracer_k_parallel",  full.names = TRUE, pattern = "^l.EMU.sources_")
# !!! Reorder the files according to 'tracer' ; files are shown by default in alphabetical order
files.EMU <- files.EMU[match(tracers, str_extract(files.EMU, "13C[:alpha:]{1,1000}"))]

for (file.i in files.EMU) {
  EMU.k <- load(file.i) %>% get()
  tracer.k <- file.i %>% str_extract("13C[:alpha:]{1,1000}")
  l.EMU.sim[[tracer.k]] <- EMU.k
}
# 
# 
# 
# # Calculate the correction factor 
# 
# # reactions involving CO2 production
# d.v.final %>% filter(str_detect(reactants, "Pyr") & str_detect(products, "AcCoA")) # CO2 made by PDH
# d.v.final %>% filter(str_detect(reactants, "Cit") & str_detect(products, "aKG"))   # CO2 made by IDH
# d.v.final %>% filter(str_detect(reactants, "aKG") & str_detect(products, "Suc"))   # CO2 made by OGDH
# d.v.final %>% filter(str_detect(reactants, "OAA") & str_detect(products, "PEP"))   # CO2 made by PEPCK
# d.v.final %>% filter(str_detect(reactants, "Mal") & str_detect(products, "Pyr"))   # Co2 made by ME 
# d.v.final %>% filter(str_detect(reactants, "Pyr") & str_detect(products, "OAA"))   # CO2 used by PC
# 
# l.EMU.sim$'13CGlc'$Pyr.Lv_1    [2] # Pyr C1 lost to CO2 by PDH
# l.EMU.sim$'13CGlc'$Cit.Lv_6    [2] # Cit C6 lost to CO2 by IDH 
# l.EMU.sim$'13CGlc'$aKG.Lv_1    [2] # aKG C1 lost to CO2 by OGDH
# l.EMU.sim$'13CGlc'$OAA.Lv_4    [2] # OAA C4 lost to CO2 by PEPCK
# l.EMU.sim$'13CGlc'$Mal.Lv_4    [2] # OAA C4 lost to CO2 by ME
# l.EMU.sim$'13CGlc'$CO2.Blood_1 [2] # CO2 used by PC 
# 
# 
# # these metabolites each has its specific carbon involved in CO2 mass balance
# compound.carbonPosition        <- c(    1,     6,     1,     4,      4,     1)
# names(compound.carbonPosition) <- c("Pyr", "Cit", "aKG", "OAA",  "Mal", "CO2")
# 
# # generate EMU carbon index that exchange with CO2
# d.EMU_exchCO2 <- d.flux.CO2 %>% select(-c(transitions, contains(".C"))) %>% 
#   
#   # if CO2 is reactant (i.e., for PC reaction), e.g., with reactant being "Pyr.Lv+CO2.Blood", keep only CO2.Blood
#   mutate(reactants = ifelse(str_detect(reactants, "\\+CO2.Blood"), "CO2.Blood", reactants)) %>% 
#   
#   mutate(compound = str_extract(reactants, "[[:alpha:][2]]{1,10}(?=\\.)")) %>%  # extract compound name
#   mutate(carbonIndex = compound.carbonPosition[compound]) %>%  # add index of carbon in the compound that exchanges with CO2
#   mutate(compartment = str_extract(reactants, "[:alpha:]{1,10}$")) %>% # add compartment
#   mutate(EMU_exchCO2 = str_c(compound, ".", compartment, "_", carbonIndex), .keep = "unused")  # create EMU full names that exchanges with CO2
# d.EMU_exchCO2
# # check duplication
# if (sum(duplicated(d.EMU_exchCO2)) >= 1) stop("Duplication found. This error greatly distorts the correction marix!")
# 
# # extract the single carbon EMU labeling from each tracer experiment
# d.singleCarbonEMU <- tibble(.rows = 0)
# for (tracer.k in tracers) {
#   # !! here unique is used, as both CO2 is by PC in both liver and KD
#   l <- sapply( unique(d.EMU_exchCO2$EMU_exchCO2), function(name) l.EMU.sim[[tracer.k]][[name]][2] ) # label of the single carbon EMUs
#   t.k <- tibble(tracer = tracer.k,
#                 EMU_exchCO2 = names(l),
#                 label =  l)
#   d.singleCarbonEMU <- bind_rows(d.singleCarbonEMU, t.k)
# }
# 
# d.singleCarbonEMU
# # check duplication
# if (sum(duplicated(d.singleCarbonEMU)) >= 1) stop("Duplicated rows found. This error greatly distorts the correction marix!")
# 
# 
# 
# # combine with the flux dataset
# d.EMU_exchCO2.all <- d.singleCarbonEMU %>% left_join(d.EMU_exchCO2) %>% 
#   mutate(flux.13C = flux * label * coef,
#          flux.12C = flux * coef) # flux 12C is the signed flux
# d.EMU_exchCO2.all
# 
# 
# 
# 
# 
# # calculate the 13CO2 calculated from only IDH and OGDH
# # total net 13CO2 produced from all pathways
# x1 <- d.EMU_exchCO2.all %>% 
#   group_by(tracer, tissue) %>% 
#   summarise(flux.allPath.13C = sum(flux.13C), # 13CO2 production
#             flux.allPath.12C = sum(flux.12C))       # the endogenous total CO2 production
# 
# # total 13CO2 produced from only IDH and OGDH
# x2 <- d.EMU_exchCO2.all %>% 
#   filter(str_detect(reactants, "aKG") | str_detect(reactants, "Cit")) %>% 
#   group_by(tracer, tissue) %>% 
#   summarise(flux.IDH_OGDH.13C  = sum(flux.13C),
#             flux.IDH_OGDH.12C = sum(flux))
# 
# d.correct.factors <- left_join(x1, x2) %>% 
#   mutate(correct.fctr.13C = flux.allPath.13C  / flux.IDH_OGDH.13C ) %>% 
#   mutate(correct.fctr.12C = flux.allPath.12C / flux.IDH_OGDH.12C ) %>% 
#   arrange(tissue) %>% 
#   ungroup()
# 
# d.correct.factors.simplified.tidy <- d.correct.factors %>% 
#   select(tracer, tissue, contains("correct.fctr")) %>% 
#   pivot_longer(-c(tracer, tissue), names_to = "state", values_to = "factor")
# 
# plt.correction.factor <- d.correct.factors.simplified.tidy %>% 
#   mutate(tracer = factor(tracer, levels = tracers)) %>% 
#   ggplot(aes(x = tracer, y = factor, fill = tissue)) +
#   geom_col(position = "dodge", color = "black", width = .7, linewidth = .3) +
#   scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, .1)),
#                      breaks = seq(0, 2, .1)) +
#   geom_hline(yintercept = 1, linetype = "dashed") +
#   theme_classic() +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1),
#         panel.grid.major.y = element_line()) +
#   facet_wrap(~state)
# 
# plt.correction.factor
# 
# ggsave("./plots/correction factor.pdf", height = 4, width = 8)
# 
# 
# 
# 
# # correction factor for 13CO2
# d.correct.factor.13C <-  d.correct.factors %>% 
#   select(tracer, correct.fctr.13C, tissue) %>% 
#   pivot_wider(names_from = tissue, values_from = correct.fctr.13C)
# 
# m.correct.factor.13C <- d.correct.factor.13C %>% select(-tracer) %>% as.matrix()
# rownames(m.correct.factor.13C) <- d.correct.factor.13C$tracer
# 
# 
# # correction factor for total 12CO2
# d.correct.factor.12C <-  d.correct.factors %>% 
#   select(tracer, correct.fctr.12C, tissue) %>% 
#   pivot_wider(names_from = tissue, values_from = correct.fctr.12C)
# 
# m.correct.factor.12C <- d.correct.factor.12C %>% select(-tracer) %>% as.matrix()
# rownames(m.correct.factor.12C) <- d.correct.factor.12C$tracer
# 
# 
# # fasted
# if (state_fasted.or.refed == "fasted") {
#   m.correct.factor.13C.fasted <- m.correct.factor.13C
#   m.correct.factor.12C.fasted <- m.correct.factor.12C
#   save(
#     m.correct.factor.13C.fasted, 
#     m.correct.factor.12C.fasted, 
#     optimal_cost,
#     file = "4_correct_factor_fasted.RData"
#   )
# }
# 
# 
# # refed
# if (state_fasted.or.refed == "refed") {
#   m.correct.factor.13C.refed <- m.correct.factor.13C
#   m.correct.factor.12C.refed <- m.correct.factor.12C
#   save(
#     m.correct.factor.13C.refed, 
#     m.correct.factor.12C.refed, 
#     optimal_cost,
#     file = "4_correct_factor_refed.RData")
# }



d.v.final



# sketch out liver fluxes

# create metabolites coordinate

# create metabolites coordinate
d.coord <- tribble(
  ~reactants, ~x,  ~y,
  "Glc.Lv",     15, 15,
  "Glc.Blood",  25, 15,
  "Glycogen.Lv",   15, 20,
  "GAP.Lv",     10,   10,
  "DHAP.Lv",    20,   10,
  "Glycerol.Blood", 20, 5,
  
  "Lys.hp",     20, 0,
  "Trp.hp",     20, -4,
  "Tyr.hp",     20, -10,
  "Thr.hp",     0, -10,
  
  
  "Palm.Blood", 25, - 20,
  "Ole.Blood",  22, - 23,
  "Lino.Blood", 19, - 25,
  "srcAcCoA",   15, -27,
  
  "Pyr.Lvc",     0,   0,
  "Lac.Blood",   10, 0,
  "Ala.Blood",   -10, 10,
  "Pyr.Lvm",     0,   -5,
  "AcCoA.Lv",   10, -10,
  "Cit.Lv",     10, -18,
  "aKG.Lv",      5, -30,
  "Gln.Blood",   5, -35,
  "Suc.Lv",    -5, -30,
  "Mal.Lv",    -12, -20,
  "OAA.Lv",    -10, -10,
  "PEP.Lv",      0,  10,
  "AcAct.Lv",    15, -30,
  "HB.Blood",    25, -30,
  "sink",        20, -35
  # "protein",  -30,-35, 
)

d.Lv <- d.v.final %>%
  filter(reactions %>% str_detect("Lv")) %>%
  mutate(reactants = str_remove(reactants, "\\+CO2.Blood")) %>%
  mutate(products = str_remove(products, "\\+CO2.Blood")) %>%
  select(-c(transitions, contains(".C")))

# if there is a '+' in reactants or products, split into two separate rows
for (i in 1:nrow(d.Lv)) {
  # i=9
  d.Lv.i <- d.Lv[i, ]
  if (d.Lv.i$reactants %>% str_detect("\\+")){
    d.Lv <- d.Lv[-i, ] # remove that row
    d.new  <- d.Lv.i %>% separate(reactants, into = c("r1", "r2"), sep = "\\+")
    d.new1 <- d.new %>% select(-r2) %>% rename(reactants = r1)
    d.new2 <- d.new %>% select(-r1) %>% rename(reactants = r2)
    d.Lv <- bind_rows(d.Lv, d.new1, d.new2)
  }
}

d.Lv


p.map <- d.Lv %>%
  left_join(d.coord) %>% rename(reactants.x = x, reactants.y = y) %>%
  left_join(d.coord %>% rename(products = reactants)) %>% rename(product.x = x, product.y = y) %>%
  # filter(complete.cases(.)) %>%
  # calculate flux position
  mutate(f.x = (reactants.x + product.x)/2) %>%
  mutate(f.y = (reactants.y + product.y)/2) %>%
  # slice_head(n = 2) %>%
  
  ggplot(aes(x = reactants.x, reactants.y)) +
  geom_segment(aes(x = reactants.x, y = reactants.y, xend = product.x, yend = product.y)) +
  geom_label(aes(label = reactants), fill = "turquoise", fontface = "bold") +
  geom_label(aes(label = products, x = product.x, y = product.y), fill = "turquoise", fontface = "bold") +
  ggrepel::geom_text_repel(aes(label = round(flux), x = f.x, y = f.y), force = .1,
                           box.padding = unit(0, "pt"), color = "tomato", fontface = "bold") +
  scale_x_continuous(expand = expansion(add = 3), name = NULL) +
  scale_y_continuous(expand = expansion(add = 3), name = NULL) +
  theme_minimal() 
p.map
# theme_void()




# Plot together
# cost + liver flux map
p.top <- cowplot::plot_grid(p.cost, p.lv.M, p.map, rel_widths = c(4, .5, 5), nrow = 1, align = "h")                   

# obs vs sim
p.btm <- p02.facet  + facet_wrap(~tracer, nrow = 2)  + theme(legend.position = "right") +
  guides(color = guide_legend(ncol = 1, override.aes = list(size = 5, fontface = "bold"))) 

cowplot::plot_grid(p.top, p.btm, rel_heights = c(3, 1.8), ncol = 1)



# calcualte time used
time.End.iteration <- Sys.time()
time.used <- time.End.iteration - time.Start.iteration


if (! exists("repeat_index")) {repeat_index <- 0}

ggsave(filename = str_c(
  "./plots/converging_history/result ", repeat_index, "  -  ", round(time.used, 1), " min.pdf"), 
  height = 13, width = 14)





# Log up
sink("./plots/converging_history/best converge.txt", append = T)
aaaaa <- paste0("repeat No. ", repeat_index, ";  ", 
                optim.step, " steps; ",
                round(time.used, 2), " min; ",
                "cost: ", slice_tail(d.cost.iterations, n = 1)$cost %>% round(8))
aaaaa
sink()


 
# mark completion
beepr::beep(2) 






# Log up
message("record this converging in the log text in 'converging_history' folder")
#sink("./plots/converging_history/log converge.txt", append = T)
#sink()

aaaaa <- paste0(
  "iteration steps ", optim.step, "; Time spent to converge ",
  round(time.used, 2), " min; ",
  "cost at terminal step ", slice_tail(d.cost.iterations, n = 1)$cost %>% round(8))

# export explicit to the log text
cat(aaaaa, "\n", file = "./plots/converging_history/log converge.txt", append = TRUE)

