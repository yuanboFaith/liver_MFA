# Based on the carbon transition table (two input sheets, global and tracing design, in the global environment)
# we'll create the stoichiometric matrix, identify a set of free fluxes, perform EMU decomposition, calculate A and B matrices and their derivatives.
# The stoic matrix and free fluxes are identical in all tracing experiments. 
# The matrices A and B and their derivatives are constants for a given tracing experiment and EMU size, independent of optimization iterations. 
# In the end, the script outputs all needed objects for iterative optimization procedure (in a different R file).


# This script does NOT require loading 'supplement_function.RData'. 
# However, it requires the tracer infusion rate, which can be input manually, or imported from '13C_tracing_labeling_data.RData' for the infusion parameters


# The proper function of the program is limited to the following:
# (1) Stoichiometric reactions with coefficient 1
# (2) Metabolites with NO more than 9 atoms
# (3) NO more than two reactants (involved in carbon transition) in the same reaction (limited to Cauchy product of only two EMUs). This limit is usually not a concern for vast majority of reactions involving carbon transition. 

rm(list = ls())

library(readxl)
library(tidyverse)
# library(Ryacas)

# set path to the folder of the current active script
rstudioapi::getActiveDocumentContext()$path %>% dirname() %>% setwd(); getwd()

load("1_Supplement_functions.RData")

filePath <- "LvM_fasted.xlsx"

# sheet of global carbon transitions of the metabolic network (independent of tracing design)
reaction_data0 <- read_excel(filePath, sheet = "LvM") %>% select(1:3)
tail(reaction_data0)

# sheet of carbon transition of tracer input; in 'func.stoicEMU()', a row will be selected based on the specified tracer name, and bind to the bottom of the global transition sheet
infusion_data <- read_excel(filePath, sheet = "tracing") %>% select(1:4)

load("../data/cleaned_labeling_data.RData")


func.stoicEMU <- function(
    
  # input value must be identical to the name used in the global carbon transition table
  # input values serves two roles: (1) select a row in the tracing sheet, and bind to global sheet; (2) appended to the global substrate vector     
  myTracer,             # e.g., myTracer = "13CGlc",  
  tracerInputRate) {    # e.g., tracerInputRate = 20
  
  # check
  # myTracer= "13CGlc" ;   tracerInputRate = 20
  # myTracer= "13CLac" ;   tracerInputRate = 20
  # myTracer= "13CMal.M" ; tracerInputRate = 0
  # myTracer= "13ChpAAb" ; tracerInputRate = 0
  
  # combine 
  reaction_data <- bind_rows(reaction_data0, 
                             filter(infusion_data, tracer == myTracer) %>% select(-tracer) )
  
  n.reactions <- nrow(reaction_data) # total number of reactions, including tracer input and fumarate scrambling flux
  fluxes <- paste0("v", 1:n.reactions)
  
  
  reaction_data <- reaction_data %>% 
    # remove white space
    mutate(reactions = str_remove_all(reactions, " "),
           transitions = str_remove_all(transitions, " ")) %>% 
    # names of reactants and products as separate columns
    separate(reactions,   into = c("reactants",   "products"),   sep = "->", remove = F) %>% 
    separate(transitions, into = c("reactants.C", "products.C"), sep = "->", remove = F)
  
  reaction_data
  
  # check if there are duplicated rows
  e <- reaction_data %>% select(-enzyme)
  if ( sum(duplicated(e)) >=1 ) {
    duplicated_rows <- e[duplicated(e) | duplicated(e, fromLast = TRUE), ]
    print(duplicated_rows)
    stop(str_c("Above rows are duplicated"))
  }
  
  # create the stoichiometric matrix
  # create mass balance for molecules being both a product and a reactant in the network
  print("Preparing the stoichiometric matrix...")
  
  reactants.all <- unlist(sapply(reaction_data$reactants, function(x) strsplit(x, "\\+"))) %>% unique()
  products.all <- unlist(sapply(reaction_data$products, function(x) strsplit(x, "\\+"))) %>% unique()
  intermediary.all <-  intersect(reactants.all, products.all) # the intermediary molecules being both products and reactants in different reactions
  n.intermediary <- length(intermediary.all) # number of intermediary molecules
  
  Stoich <- matrix(0, nrow = n.intermediary, ncol = n.reactions); 
  rownames(Stoich) <- intermediary.all
  colnames(Stoich) <- 1:n.reactions
  Stoich # row names indicate the intermediary the mass balance eq. is based on
  
  for (i in 1:n.intermediary){ # !!! only handles coefficient 1 for maximum of two reactants and two products
    # check
    # i=1
    for (j in 1:nrow(reaction_data)){ # check each reaction if the given intermediary is involved as reactant or product
      # check
      # j=5
      x <- reaction_data[j, ]
      rs <- x$reactants %>% str_split_1("\\+") # all reactant(s) (one or two) involved in reaction j
      ps <- x$products %>% str_split_1("\\+") # all product(s) (one or two) involved in reaction j
      
      # if found the intermediary as a reactant, than mark -1
      if (intermediary.all[i] %in% rs) { Stoich[i, j] <- -1 } 
      # if found the intermediary as a product, than mark 1
      if (intermediary.all[i] %in% ps) { Stoich[i, j] <- 1 } 
    }
  }
  Stoich
  
  # !!! Important: 
  # (1) When CO2 recycling is considered, the product CO2 SHOULD be included the stoichiometric matrix.
  # Stoich <- Stoich[rownames(Stoich) != "CO2", ]
  
  # # (2) Remove the Stoich last row and second-to-the-last column corresponding to fumarate2 exchange flux, as it does not affect the mass balance; however, it should be reflected in the EMU network.  
  # Stoich <- Stoich[-nrow(Stoich), - (ncol(Stoich)- 1) ]
  
  # (3) Remove the Stoich last column corresponding to 13C tracer input; thus the Stoich is identical in different 13C-tracing experiments
  Stoich <- Stoich[, - (ncol(Stoich)) ]
  Stoich
  tail(Stoich)
  
  
  
  
  
  # Identify a set of free fluxes from the Stoich matrix.
  p <- Stoich %>% pracma::rref() %>% qr() %>% .$pivot; print(p)
  n.free <- ncol(Stoich) - pracma::Rank(Stoich) # number of free independent fluxes
  n.dep <- ncol(Stoich) - n.free                # number of dependent fluxes
  
  index.dep        <- p[1:n.dep]                # column index of dependent fluxes (corresponding to the left side identity matrix in the reduced row echelon form of Stoich)
  index.freeFluxes <- p[(n.dep+1) : length(p)]  # column index of free independent fluxes
  index.freeFluxes
  
  rank_Stoich <- qr(Stoich)$rank
  n.freeFluxes <- length(index.freeFluxes)
  
  S_dep  <- Stoich[, index.dep] # S matrix associated with dependent fluxes
  S_free <- Stoich[, index.freeFluxes] # S matrix associated with free fluxes
  
  # m.freeFlux_to_fullSetFlux <- -solve(S2) %*% Su # matrix mapping free fluxes to dependent fluxes
  
  # For the simple TCA flux network (Pyr, AcCoA, OAA, Cit, aKG, Suc), the associated stoic is square, and S2 is not square, not inversible; use pseudo-inverse instead
  m.freeFlux_to_fullSetFlux <- - func.pseudo_inverse_svd(S_dep)  %*% S_free %>% round(8)
  rownames(m.freeFlux_to_fullSetFlux) <- colnames(S_dep)
  
  
  # matrix mapping free fluxes to the full set of fluxes
  m.freeFlux_to_fullSetFlux <- rbind(
    diag(1, n.free), # diag matrix, mapping free to itself
    m.freeFlux_to_fullSetFlux) # dependent fluxes
  
  rownames(m.freeFlux_to_fullSetFlux)[1:n.freeFluxes] <- index.freeFluxes # names for rows mapping free fluxes to itself
  
  # Update the row order to original sequential number
  m.freeFlux_to_fullSetFlux <- m.freeFlux_to_fullSetFlux[as.numeric(rownames(m.freeFlux_to_fullSetFlux)) %>% order(), ]
  
  # check that the free to full flux conversion matrix should be a null space of Stoich
  if (all(Stoich %*% m.freeFlux_to_fullSetFlux == 0) == F) {
    stop("The mapping matrix converting free to full set of fluxes should be a null space of the stoichiometric matrix")
  }
  
  tail(m.freeFlux_to_fullSetFlux)
  # for suc - suc2 scramble reactions, express them as an arbitrarily big numbers, e.g., 100 fold of sum of all free fluxes
  n.scrambleReactions <- reaction_data$reactants %>% str_detect("Suc2") %>% sum() # number of suc - suc2 scramble reactions
  
  # turn last 'n.scrambleReactions' rows to 100 fold of sum of all free fluxes
  # the last row corresponding to tracer mass balance has been removed earlier
  m.freeFlux_to_fullSetFlux[(nrow(m.freeFlux_to_fullSetFlux)-n.scrambleReactions+1):nrow(m.freeFlux_to_fullSetFlux), ] <- 10 
  
  
  
  # When working with free fluxes, express A and B matrices fluxes with free ones (except the fumarate scrambling flux and tracer input flux as fixed values)
  # Note that the last reaction, the fumarate scrambling flux, is not included here at this step.  
  v.expressWithFreeFlux <- c()
  for (i in 1:nrow(m.freeFlux_to_fullSetFlux)){
    x <- paste0(m.freeFlux_to_fullSetFlux[i, ], paste0("v", index.freeFluxes)); x # express fluxes as a linear combination of free fluxes
    
    x <- x[!str_detect(x, "^0v")]     # remove fluxes with zero coefficient
    x <- str_replace(x, "^1v", "v")   # no need to show the coefficient 1
    x <- str_replace(x, "^-1v", "-v") # show -v instead of -1v
    x <- paste(x, collapse = "+") %>% str_replace_all("\\+\\-", "-") # connect all together, and turn '+-' to '-'
    
    v.expressWithFreeFlux[i] <- x # record this flux
    names(v.expressWithFreeFlux)[i] <- rownames(m.freeFlux_to_fullSetFlux)[i] 
  }
  
  v.expressWithFreeFlux %>% unname()
  
  
  #################### Add the fumarate scrambling flux and tracer infusion flux 
  #################### tracerInfusionRate <- as.character(tracerInputRate); # tracer input rate, from the manual input
  
  # add tracer infusion rate
  v.expressWithFreeFlux <- c(v.expressWithFreeFlux, tracerInputRate)
  
  names(v.expressWithFreeFlux)[n.reactions] <- n.reactions # n.reactions is the total biophyical reaction number, including succinate scrambling and tracer infusion
  v.expressWithFreeFlux
  
  
  
  # create the EMU framework
  print("Performing EMU decomposition...")
  
  target.EMU <- c(
    # Blood
    "Glc.Blood_123456", "Lac.Blood_123", "Ala.Blood_123", "Glycerol.Blood_123", "HB.Blood_1234", "Gln.Blood_12345",
    "Palm.Blood_12", "Ole.Blood_12", "Lino.Blood_12",
    # tissue
    "Suc.Lv_1234", "Mal.Lv_1234" # , "aKG.Lv_12345" # "Cit.Lv_123456", 
    # "Suc.M_1234", "Mal.M_1234",
    
    # CO2
    # "CO2.Blood_1"
  )
  
  
  
  new.EMU <- target.EMU # start from target EMU; it'll keep growing to store newly identified EMU
  
  # the substrate that feeds into the system
  substrate <- c("srcAcCoA", "Glycogen.Lv", "Glycogen.M", "TAG.W", "protein") # global
  substrate <- c(substrate, myTracer) # add the tracer input
  
  # add amino acid input from the portal system
  portalAAs <- c(
    "Ala.hp", "Ser.hp", "Thr.hp", "Trp.hp", "Ile.hp", "Leu.hp",
    "Lys.hp", "Phe.hp", "Tyr.hp", "Arg.hp", "Glu.hp", "Gln.hp", "Val.hp", "Asp.hp")
  substrate <- c(substrate, portalAAs)
  
  
  # !!! IMPORTANT: for 13C dietary protein tracing, add the arbitrary tracer (placeholder) to the substrate,
  # so as not to search the carbon transition in the metabolic network; 
  # the tracer has to be a intermediary metabolite; otherwise if being an arbitrary name, include it as substrate (the starting material)
  if (str_detect(myTracer, "13ChpAA")) substrate <- c(substrate, myTracer)
  
  
  decomposition <- c()
  
  # start the search with product EMU
  # !!! not able to handle non-unit coefficient in the while loop below
  
  i <- 1
  # sink("EMU decomposition.txt")
  
  while (i <= length(new.EMU)){
    products.EMU <- new.EMU[i]
    print(paste0("Check the current ", i, "th EMU, ", products.EMU, ", as the product EMU"))
    
    # check
    # products.EMU <- "Cit_12345"
    
    product.name <- products.EMU %>% str_extract("^[^_]+") # e.g., "Fum2" will match "Fum2"; match anything before the underscore
    products.index <- products.EMU %>% str_extract_all("(?<=_)\\d+", simplify = T) %>% str_split_1("") %>% as.numeric() # e.g., match digits after the underscore; e.g., "Fum2_123" matches the vector of 1, 2, 3
    
    # reactions that make the specified product EMU
    # search exactly the specified product; 
    # e.g., Fum would be selected, while Fum2 would not be selected
    # e.g., X.m will match X, not X.m
    d.i <- reaction_data %>% filter(str_detect(products, pattern = paste0("\\b", product.name, "(?!\\.)\\b"), )) 
    
    print(paste0("The following reaction makes ", products.EMU, ":"))
    print(d.i)
    
    # iterate through each reaction to make the specified product EMU
    for (j in 1:nrow(d.i)){
      d.j <- d.i[j, ]
      
      # names of each individual reactants and products as vector elements
      reactants <- str_split_1(d.j$reactants, "\\+")
      products  <- str_split_1(d.j$products, "\\+")
      
      
      # carbon transition of each individual reactants and products
      transition.reactants <- str_split_1(d.j$reactants.C, "\\+") 
      transition.products  <- str_split_1(d.j$products.C, "\\+")
      
      # compounds' names as the name of vector elements of carbon transitions 
      names(transition.reactants) <- reactants
      names(transition.products) <- products
      
      # split product carbon letters into individual vector elements
      products.split <- transition.products[product.name] %>% str_split_1(pattern = "") 
      # selected product carbon position
      products.split <- products.split[products.index] 
      
      # check origin of product carbons from the reactants
      n.reactants <- length(reactants)
      
      # if there is only one reactant
      if (n.reactants == 1) {
        print(paste0("In row ", j, ", only one reactant found to make the product EMU ", products.EMU))
        # which carbons letters in the first reactant are found in the product carbon letters?
        reactant.split <- transition.reactants %>% str_split_1(pattern = "") 
        reactant.index <- which(reactant.split %in% products.split) 
        
        # collect newly produced EMU of the first reactant 
        reactant.EMU <- paste0(reactants, "_", paste(reactant.index, collapse = ""))
        print(paste("Newly generated reactant EMU", reactant.EMU))
        
        # collect only if NOT present in prior iterations, and not being a network substrate
        if ( (!reactant.EMU %in% new.EMU) & !(names(transition.reactants) %in% substrate) ) { 
          new.EMU <- append(new.EMU, reactant.EMU) 
        }
        
        # create the EMU transition to create the given iterated product EMU
        decomposition.i <- paste(reactant.EMU, "->", products.EMU) 
      }
      
      # If there are two reactants 
      # Note: one or both two reactants may be a carbon source
      if (n.reactants == 2){
        print(paste0("In row ", j, ", two reactants found. Check each of them either being a source EMU:"))
        # which carbons letters in the first reactant are found in the product carbon letters?
        reactant.1.split <- transition.reactants[1] %>% str_split_1(pattern = "") 
        reactant.1.index <- which(reactant.1.split %in% products.split) 
        
        # proceed if a product carbon is found in the first reactant
        if(!is_empty(reactant.1.index)){
          # collect newly produced EMU of the first reactant 
          reactant.1.EMU <- paste0(reactants[1], "_", paste(reactant.1.index, collapse = ""))
          print(paste(reactant.1.EMU, "is a newly found source reactants EMU"))
          # collect only if NOT present in prior iterations, and not being a network substrate
          if ( (!reactant.1.EMU %in% new.EMU) & !(names(transition.reactants[1]) %in% substrate) ) { 
            new.EMU <- append(new.EMU, reactant.1.EMU) 
          }
        }
        
        # which carbons letters in the second reactant are found in the product carbon letters?
        reactant.2.split <- transition.reactants[2] %>% str_split_1(pattern = "") 
        reactant.2.index <- which(reactant.2.split %in% products.split) 
        
        # proceed if a product carbon is found in the second reactant
        if (! is_empty(reactant.2.index)){
          # collect newly produced EMU of the second reactant (if not present in prior iterations)
          reactant.2.EMU <- paste0(reactants[2], "_", paste(reactant.2.index, collapse = ""))
          print(paste(reactant.2.EMU, "is a newly found source reactants EMU"))
          if ( (!reactant.2.EMU %in% new.EMU) & !(names(transition.reactants[2]) %in% substrate) ) { 
            new.EMU <- append(new.EMU, reactant.2.EMU) 
          }
        }
        
        # Create the EMU transition to create the given iterated product EMU
        if (
          # if only reactant 1 is a source
          !is_empty(reactant.1.index) & is_empty(reactant.2.index)
        ){
          decomposition.i <- paste(reactant.1.EMU, "->", products.EMU)   
        } else if (
          # if only reactant 2 is a source
          is_empty(reactant.1.index) & !is_empty(reactant.2.index)
        ) {
          decomposition.i <- paste(reactant.2.EMU, "->", products.EMU)   
        } else if (
          # both reactants 1 and 2 are a source
          !is_empty(reactant.1.index) & !is_empty(reactant.2.index)
        ) {
          decomposition.i <- paste(reactant.1.EMU, "+", reactant.2.EMU, "->", products.EMU)   
        }
        
      }
      
      decomposition <- append(
        decomposition, 
        # mark the index of reaction
        paste0("#", d.j$enzyme, ": ", decomposition.i))
    }
    print(paste0("complete the ", i,  "th EMU."))
    cat("\n\n\n")
    
    i <- i + 1
  }
  
  # sink()
  
  d.EMU.reaction <- tibble(EMU.reaction = decomposition) %>% 
    mutate(EMU.reaction = str_remove_all(EMU.reaction, " ")) %>% 
    # enzyme names
    separate(EMU.reaction, into = c("enzyme", "EMU.reaction"), sep = ":") %>% 
    mutate(enzyme = str_remove(enzyme, "#")) %>% 
    # separate columns for EMU of reactants and products
    separate(EMU.reaction, into = c("reactants", "products"), sep = "->", remove = F) %>% 
    mutate(size = str_extract_all(products, "(?<=_)\\d+") %>% as.numeric() %>% str_count("\\d")) %>% # count number of digits (atoms) after the underscore; !!! can't handle molecules with 9 or more atoms
    # important: arrange in ascending order of the EMU size
    # this will be saved in order in the l.AXBY list
    arrange(size) 
  
  d.EMU.reaction # -----
  
  
  
  
  
  # print the needed substrate EMU names
  x <- c()
  for (a in substrate) { x <- append( x,  filter(d.EMU.reaction, reactants %>% str_detect(a))$reactants %>% unique() ) }
  x
  
  ##### obsolete. only needed when CO2 is not an intermediate (i.e., not recycled) !!! need manual check; for bicarbonate required reaction, keep only the bicarbonate
  # x[x %>% str_detect("CO2.Blood_1")] <- "CO2.Blood_1"
  # x
  
  
  # substrate EMU, needed for derivative computation; this list will extend by adding newly produced EMUs (the X matrix)
  l.EMU.substrates <- list()
  
  # Create the IDV for the substrate EMU; no labeling for entire set, except the tracers
  for (name in x) {
    digits <- str_extract(name, pattern = "(?<=_)\\d+") # Extract the digits after the underscore 
    vector <- c(1, rep(0, nchar(digits))) # Create the vector with initial 1 and a length determined by the number of digits
    l.EMU.substrates[[name]] <- vector # Add this vector to the list with the name of the element
  }
  
  # for fully labeled tracers (names started with "13C), turn IDV's last digit to 1, else to zero
  for (i in 1:length(l.EMU.substrates)){
    if (names(l.EMU.substrates[i]) %>% str_detect("^13C")) {
      l.EMU.substrates[[i]] <- c( rep(0, length(l.EMU.substrates[[i]])-1 ), 1 )
    }
  }
  
  
  # if the tracer is portal amino acids, update their IDV first and last digit to the measured value
  if (myTracer %>% str_detect("13ChpAA")) {
    for (i in 1:length(l.EMU.substrates)){
      
      i_EMU <- names(l.EMU.substrates[i]) # get the EMU
      
      if (i_EMU %>% str_detect("\\.hp")) {
        
        # if (i_EMU == "Ser.hp_1") stop()
        # stop()  
        print(i_EMU)
        
        # extract this specific amino acid's labeling in portal blood
        AA.i <- i_EMU %>% str_extract("[a-zA-Z]{1,8}(?=\\.hp)") # extact compound name without EMU notation
        
        # get the M+full isotopologue labeling
        labeling.i <- d.13C.refed.hpAA %>% 
          filter(Compound == AA.i) %>% 
          filter(tissue == "hp") %>% 
          filter(Infusate == myTracer) %>%  # each tracer corresponds to a mouse
          filter(C_Label == max(C_Label)) %>%  # select the full label
          pull(labeling)
        
        print(labeling.i)
        
        # l.EMU.substrates[[i]] <- c(
        #   rep(0, length(l.EMU.substrates[[i]])-1 ), 
        #   labeling.i) # this labeling is extracted from the measured portal AA labeling specific to each AA and mouse ID 
        
        l.EMU.substrates[[i]] [length(l.EMU.substrates[[i]])] <- labeling.i # update the last digit
        l.EMU.substrates[[i]] [1] <- 1-labeling.i # update the first digit
      }
    }
  }
  
  
  
  
  # derivative matrix of the EMU; for substrates, the derivative relative to all fluxes are 0
  print("Calculating derivatives of substrate EMUs...")
  
  l.deriv_EMU.substrates <- sapply(
    l.EMU.substrates, function(x)str_replace_all(x, pattern = ".", replacement = "0") %>% as.numeric)
  
  # Updated with respect to each free flux (all zeros, as they are independent of fluxes)
  # this format is compatible with the EMU derivative formats, and makes it easier to program in the optimization loops
  x <- list()
  i <- 1
  for ( a in 1:length(l.deriv_EMU.substrates) ) {
    for (f in index.freeFluxes) {
      x[i] <- l.deriv_EMU.substrates[a] # keep an entry for each combination of EMU substrates and free fluxes
      names(x)[i] <- paste0( names(l.deriv_EMU.substrates[a]), "-v", f  ) # rename using names of EMU substrate and free flux
      i <- i+1
    }
  }
  
  l.deriv_EMU.substrates <- x
  
  
  
  
  # create EMU mass balance equations 
  print("Calculating matrices A and B, and express fluxes as a function of free fluxes...")
  print("yacas combine like terms element-wise in matrices A and B will take a few moments.")
  
  l.AXBY <- list()
  n.EMU_sizes <- n_distinct(d.EMU.reaction$size)
  
  for (i in 1:n.EMU_sizes) { # loop through different EMU sizes
    
    print(paste("Working on EMU size of", i)); cat("\n")
    # check
    # i=5
    size.i <- d.EMU.reaction %>% filter(size == unique(d.EMU.reaction$size)[i])
    
    # create a single full set of matrix multiplication of Ax = 0, without separating the unknown and known fluxes
    
    EMU.i <- c(size.i$reactants, size.i$products) %>% unique()
    n.row <- n_distinct(size.i$products) # matrix row number equal to number of unique products in the sub-network
    n.col <- length(EMU.i) # column number equal to total number of unique reactants and products in the sub-network
    
    
    # write out the balance equations 
    # This part is not used in later code. It's here to help troubleshooting
    # this section also helped me clear my mind of the later algorithm to create the stoichiometric matrix
    size.i.balanced <- size.i %>% 
      arrange(products) %>% 
      
      # for each product, multiply the reaction rate with the reactant EMU, and sum the reactions up
      mutate(incomeFlux = paste0("v", enzyme, "⋅", reactants)) %>% 
      group_by(products) %>% 
      mutate(incomeFlux.all = str_flatten(incomeFlux, collapse = "+")) %>% 
      
      # for each product, the outgoing flux equals the total income flux
      mutate(outcomeFlux.all = str_flatten(paste0("v", enzyme), last = "+")) %>% ungroup() %>% # intermediary step
      mutate(outcomeFlux.all = paste0("(", outcomeFlux.all, ")⋅", products)) %>% 
      
      # complete balanced equation, with income flux = outcome flux
      # note: by this step, for rows of the same product, the eq value is the same (the final balanced equation)
      mutate(eq = paste0(incomeFlux.all, "-", outcomeFlux.all, "= 0"))
    
    
    # assign a position EMU has the position in the matrix
    position <- 1:length(EMU.i)
    names(position) <- EMU.i
    
    # full set of stoichiometric matrix, with S ⋅ EMU = 0 
    S <- matrix(0, nrow = n.row, ncol = n.col)
    
    # iterate through each product (no. unique products = no. reactions = row number)
    for (k in 1:n.row){ 
      # k=1
      product.k <- unique(size.i$products)[k] 
      print(paste("check product #", k, ":", product.k))
      
      size.i.k <- size.i %>% filter(products == product.k)
      
      # product flux (outgoing flux)
      # at kth row, at the associated column position of the product, sum up fluxes as NEGATIVE values.
      
      ### opt (1): this line is used when working with FULL set of fluxes
      ### S[k, position[product.k]] <- paste0("-v", size.i.k$enzyme, collapse = "")
      
      ### opt (2): this line is used when working with FREE fluxes, expressing fluxes with free fluxes
      x <- paste0("+", v.expressWithFreeFlux[size.i.k$enzyme], collapse = "") # first connect with plus sign, and then swap the signs
      x <- str_replace_all(x, "\\+\\-", "\\-") # replace +- with - sign
      
      x <- x %>% str_replace_all("\\-", "$") %>% # mark the original minus sign
        str_replace_all("\\+", "\\-") %>% # plus turned to minus
        str_replace_all("\\$", "\\+") # original minus turned to plus
      
      
      S[k, position[product.k]] <- x
      
      print( paste("outgoing flux (expressed as free fluxes):", S[k, position[product.k]]) )
      
      # reactant flux (incoming flux)
      # loop through the reactants
      for (l in 1:nrow(size.i.k)){
        # for each reactant, at the kth row, at the associated column position of the reactant, mark the flux as POSITIVE value
        size.i.k.l <- size.i.k[l, ]
        print(paste("check incoming flux from reactant,", size.i.k.l$reactants,  ", #", l, "/ total", nrow(size.i.k), "reactions" ))
        
        ### opt (1): this line is used when working with FULL set of fluxes
        # S[k, position [size.i.k.l$reactants]] <- paste0("v", size.i.k.l$enzyme)
        
        ### opt (2): this line is used when working with FREE fluxes, expressing fluxes with free fluxes
        S[k, position [size.i.k.l$reactants]] <- paste0(v.expressWithFreeFlux[size.i.k.l$enzyme])
        
        if (paste0(v.expressWithFreeFlux[size.i.k.l$enzyme]) == "NA") stop()
        
        print(paste0("v", size.i.k.l$enzyme, " is expressed as free flux ", S[k, position [size.i.k.l$reactants]]))
      }
      cat("\n")
    }
    
    # replace the double negative -- to a plus sign
    # use c(1, 2) to ensure that the output is a matrix; if S is a vector, it loses its dimension
    S <- apply(S, c(1, 2), str_replace, pattern = "--", replacement = "+")
    
    # print the vector
    tibble(EMU=names(position))
    
    # find the "substrates" of the sub network, i.e., a reactant that is not a product 
    # arrange substrate-containing reactions to the last rows in the EMU vector, 
    # and move their associated fluxes to the most right side in the stoichiometric matrix
    # need to use 'unique()' to remove duplicated 'size.i$reactants'; when different metabolites share the same substrate source (e.g., TAG.W for different NEFA), duplication happens
    # If without 'unique()', then for EMU size 1, Y vector may incorrectly contain non-substrate things; 
    # that is, instead of containing only glycogen, TAG, protein, tracer, i.e., the external input source, it may incorrectly contain intermediates
    substrate.i <- size.i$reactants[!size.i$reactants %in% size.i$products] %>% unique()
    substrate.i
    w <- length(substrate.i)
    
    # indices of columns in S that should be moved to the right, 
    # and rows in EMU vector that should be moved to the bottom
    indices.move <- which(names(position) %in% substrate.i) 
    
    x1 <- !(1:n.col) %in% indices.move # columns/rows not affected
    x2 <-  (1:n.col) %in% indices.move # columns/rows to be rearranged
    
    # !!! If there is only one row in the original full set S, then the cbind can mess up the dimension
    # Need to ensure the cbind result is a single row matrix
    if (nrow(S) == 1) {
      S.rearranged <- c(S[, x1], S[, x2]) %>% matrix(nrow = 1)
    } else{
      S.rearranged <- cbind(S[, x1], S[, x2])
    }
    position.rearranged <- c(position[x1], position[x2])
    
    # divide the S ⋅ EMU = 0 into AX + BY = 0
    A <- S.rearranged       [, 1: (ncol(S.rearranged)-w) ] 
    X <- position.rearranged[  1: (ncol(S.rearranged)-w) ] %>% names()
    
    B <- S.rearranged       [, (ncol(S.rearranged)-w+1):ncol(S.rearranged), drop = FALSE ] # if the matrix output is a single row, keep it as matrix, not a vector
    Y <- position.rearranged[  (ncol(S.rearranged)-w+1):ncol(S.rearranged) ] %>% names()
    
    # flip the signs of all terms in matrix B, such that AX = BY
    apply_minus_signs <- function(x) {
      # check
      # x="2v1+v9-v10"
      if (x != "0") {
        if (str_detect(x, "^\\d{0,2}-")){ # if the string starts with minus sign (w or w/o coefficient), then flip all terms
          x <- x %>% 
            str_replace_all("\\-", "$") %>% # mark the original minus sign
            str_replace_all("\\+", "\\-") %>% # plus turned to minus
            str_replace_all("\\$", "\\+") # original minus turned to plus  
        } else { # if not started with the minus sign, then add the plus sign at the beginning, and then flip all signs
          x <- paste0("+", x) %>% 
            str_replace_all("\\-", "$") %>% # mark the original minus sign
            str_replace_all("\\+", "\\-") %>% # plus turned to minus
            str_replace_all("\\$", "\\+") # original minus turned to plus  
        }
        
      }
      return(x)
    }
    
    # ensure B as a matrix to be applicable for apply function; 
    # sometimes B can be a vector of length 1
    B <- apply(as.matrix(B), MARGIN = c(1, 2), apply_minus_signs)
    
    # Check: Print them out, now with AX = BY
    # A
    # tibble(X)
    # B
    # tibble(Y)
    
    # add a multiplication sign between the coefficient and the v notation
    A <- A %>% as.matrix() %>% apply(MARGIN = 2, str_replace_all, pattern = "(?<=\\d)v", replacement = "*v") 
    B <- B %>% as.matrix() %>% apply(MARGIN = 2, str_replace_all, pattern = "(?<=\\d)v", replacement = "*v") 
    
    # combine like terms with yacas expression. This procedure is rather slow...
    # ensure the A and B are output as matrices. 
    # This is critical: If A and B are single length vector, it loses the dimension property, and would has difficulty working with 'apply()'
    # A <- A %>% as.matrix() %>% apply(MARGIN = c(1,2), function(x) yacas(paste0("Simplify(", x, ")")) %>% as.character() %>% str_remove_all(" ")) %>% as.matrix()
    # B <- B %>% as.matrix() %>% apply(MARGIN = c(1,2), function(x) yacas(paste0("Simplify(", x, ")")) %>% as.character() %>% str_remove_all(" ")) %>% as.matrix()
    
    # in the list AXBY, from the 1st to the last element corresponds to EMU of increasing size
    # the EMU size however does not necessarily correspond to the position index, as the EMU size does not have to be continuous
    l.AXBY[[i]] <- list(A = A, X = X, B = B, Y = Y)
    print(paste("complete iteration", i))
  }
  
  # Check:
  # l.AXBY
  # l.AXBY[[1]]$A
  # l.AXBY[[4]]$X
  # l.AXBY[[1]]$B
  # l.AXBY[[4]]$Y
  
  
  # Calculate derivatives of A and B matrices. 
  # For a given sub-network, the derivatives of A and B are constants. 
  print("Calculating the derivatives of matrices A and B...")
  
  l.deriv_A <- list()
  l.deriv_B <- list()
  
  j <- 1 # each index tracks the derivative of A of a given EMU size with respect to a given flux
  
  for (i in 1:length(l.AXBY)) { # i tracks the index of the EMU of different sizes; note that i is not necessarily the actual EMU size, which may be non-continuous
    
    for (flux.i in paste0("v", index.freeFluxes) ) {
      
      # derivative of matrix A. The yacas here is also a rather slow procedure...
      print(paste("calculating derivatives in matrix A with respect to", flux.i, "| ", "EMU size index", i , "|",   myTracer ))
      
      l.deriv_A[[j]] <- l.AXBY[[i]]$A %>% 
        # !!! as matrix is important; it ensures compatibility with apply() when A is a single length vector
        as.matrix() %>% 
        # !!! note 2: MARGIN = c(1,2) ensures that if the matrix input is a single row, the output maintains the single row matrix dimension 
        apply(MARGIN = c(1,2), 
              # function(x) yacas(paste0("D(", flux.i, ")", x)) %>% as.character() %>% as.numeric ### the slow version
              function(x) {
                if (x == "0") {
                  return(0)
                } else {
                  Ryacas::yac_expr(paste0("D(", flux.i, ")", x)) %>% as.character() %>% as.numeric()
                }
              }
        )
      
      names(l.deriv_A)[j] <- paste0("dA", i, "_", flux.i) # label with name: EMU size, and flux
      
      # derivative of matrix B
      print(paste("calculating derivatives in matrix B with respect to", flux.i, "| ", "EMU size index", i , "|",   myTracer ))
      
      l.deriv_B[[j]] <- l.AXBY[[i]]$B %>% 
        as.matrix() %>% 
        apply(MARGIN = c(1,2), 
              # function(x) yacas(paste0("D(", flux.i, ")", x)) %>% as.character() %>% as.numeric ### the slow version
              function(x) {
                if (x == "0") {
                  return(0)
                } else {
                  Ryacas::yac_expr(paste0("D(", flux.i, ")", x)) %>% as.character() %>% as.numeric()
                }
              }
        ) 
      
      names(l.deriv_B)[j] <- paste0("dB", i, "_", flux.i) # label with name: EMU size, and flux
      
      j <- j + 1
    }
  }
  
  # output the items
  mylist <- list(
    
    # these items should be globally identical across different 13C tracing experiments
    "Stoich" = Stoich, 
    "index.freeFluxes" = index.freeFluxes,
    "m.freeFlux_to_fullSetFlux" = m.freeFlux_to_fullSetFlux,
    "l.deriv_EMU.substrates" = l.deriv_EMU.substrates,
    "m.freeFlux_to_fullSetFlux" = m.freeFlux_to_fullSetFlux,
    "n.EMU_sizes" = n.EMU_sizes,
    "target.EMU" = target.EMU,
    reaction_data = reaction_data,
    
    "l.deriv_A" = l.deriv_A, 
    "l.deriv_B" = l.deriv_B,
    
    # unique to each tracing experiment
    "substrate" = substrate,
    "v.expressWithFreeFlux" = v.expressWithFreeFlux,
    "l.EMU.substrates" = l.EMU.substrates,
    "l.AXBY"    = l.AXBY
  ) 
  mylist %>% return()
}


save.image(file = "2_func.stoicEMU.RData")


# func.stoicEMU(myTracer = "13ChpAAb",       tracerInputRate = 0)
# l.13CGlc       <- func.stoicEMU(myTracer = "13CGlc",       tracerInputRate = 200  * .1)
# l.13CLac       <- func.stoicEMU(myTracer = "13CLac",       tracerInputRate = 380  * .1)
# l.13CGlycerol  <- func.stoicEMU(myTracer = "13CGlycerol",  tracerInputRate = 110  * .1)
# l.13CGln       <- func.stoicEMU(myTracer = "13CGln", tracerInputRate = 85.9 * .1)
# 
# # compile together
# l.stoich.EMU_allTracers <- list("13CGlc"      = l.13CGlc, 
#                                 "13CLac"      = l.13CLac,
#                                 "13CGlycerol" = l.13CGlycerol,
#                                 "13CGln"      = l.13CGln)
# 
# # make objects shared in all tracer experiments available directly from global environment for convenience
# m.freeFlux_to_fullSetFlux <- l.13CGlc$m.freeFlux_to_fullSetFlux
# Stoich                    <- l.13CGlc$Stoich
# index.freeFluxes          <- l.13CGlc$index.freeFluxes
# n.EMU_sizes               <- l.13CGlc$n.EMU_sizes
# 
# save.image(file = "l.stoich.EMU_allTracers.RData")
# 
# 
# #  ----------------------------------- Testing ----------------------------------------
# 
# 
# # the free set of fluxes must be IDENTICAL across all tracing experiments
# all.equal(l.13CGlc$index.freeFluxes, l.13CLac$index.freeFluxes)
# all.equal(l.13CGlc$index.freeFluxes, l.13CGlycerol$index.freeFluxes)
# 
# # the null space matrix mapping free to all set of fluxes must be IDENTICAL across all tracing experiment
# all.equal(l.13CGlc$m.freeFlux_to_fullSetFlux, l.13CLac$m.freeFlux_to_fullSetFlux)
# all.equal(l.13CGlc$m.freeFlux_to_fullSetFlux, l.13CGlycerol$m.freeFlux_to_fullSetFlux)
# 
# # the stoich matrix must be IDENTICAL across all tracing experiments
# all.equal(l.13CGlc$Stoich, l.13CLac$Stoich)
# all.equal(l.13CLac$Stoich, l.13CGlycerol$Stoich)
# 
# 
# # the A and B matrices across 13C experiments must be DIFFERENT
# all.equal(l.13CGlc$l.AXBY, l.13CLac$l.AXBY)
# all.equal(l.13CLac$l.AXBY, l.13CGlycerol$l.AXBY)
# 
# all.equal(l.13CGlc$l.AXBY[[3]]$A, l.13CLac$l.AXBY[[3]]$A)
# all.equal(l.13CGlc$l.AXBY[[3]]$A, l.13CGlycerol$l.AXBY[[3]]$A)
# 
# 
# # The derivatives of A and B over fluxes is the SAME!
# all.equal(l.13CGlc$l.deriv_A, l.13CLac$l.deriv_A)
# all.equal(l.13CGlc$l.deriv_A, l.13CGlycerol$l.deriv_A)
