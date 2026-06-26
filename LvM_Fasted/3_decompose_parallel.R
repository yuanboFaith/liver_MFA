rm(list = ls()) # clearing existing variable  - this function remains active in this script. 

library(tidyverse)
library(parallel)  

# set path to the folder of the current active script
rstudioapi::getActiveDocumentContext()$path %>% dirname() %>% setwd(); getwd()

# load the 'func.stoicEMU' function to perform EMU decomposition
load(file = "2_func.stoicEMU.RData")

# cleaned labeling data (from the 'data' folder)
# load(file = "/Users/boyuan/Desktop/Harvard/Research/Energy Metabolism Atlas/data/cleaned_labeling_data.RData")

time.EMU.Start <- Sys.time(); time.EMU.Start

## Set up the Parallel Environment
# Set up the cluster with four workers (adjust if necessary), and print progress in the console
cl <- makeCluster(5,  outfile = "")

# Make libraries accessible to each working node
clusterEvalQ(cl, {
  library(parallel)  
  library(readxl)
  library(Ryacas)
  library(tidyverse)
})

# Define a wrapper function for calling 'func.stoicEMU' with different parameters
call_func_stoicEMU <- function(myTracer, tracerInputRate) {
  func.stoicEMU(myTracer = myTracer, tracerInputRate = tracerInputRate)
}


# Export all global environment objects to the working nodes
# objects to be used in the parallel computation should be defined above this line
clusterExport(cl, varlist = ls())


# # cleaned labeling - infusion data (from the 'data' folder)
# load(file = "../data/cleaned_labeling_data.RData")
# d.13C.fasted %>% select(Infusate, infuse.nmol.min.g) %>% distinct()
# d.13C.refed  %>% select(Infusate, infuse.nmol.min.g) %>% distinct()


# 'state_fasted.or.refed' is a master control variablel. It may be commented out here and specified in other files
# use this parameter to dictate which state to compute about
# this directly affects 1) the tracer infusion rate and 2) the output .RData file naming; with each state having its own .RData file

state_fasted.or.refed <- "fasted" # "refed" 


# Specify the parameters for each call (each tracing experiment)

if (state_fasted.or.refed == "fasted") {
  # fasted state
  params <- list(
    "13CGlc"      = list(myTracer = "13CGlc",      tracerInputRate = 20),
    "13CLac"      = list(myTracer = "13CLac",      tracerInputRate = 58.8),
    "13CGlycerol" = list(myTracer = "13CGlycerol", tracerInputRate = 10),
    "13CAla"      = list(myTracer = "13CAla",      tracerInputRate = 20),
    "13CGln"      = list(myTracer = "13CGln",      tracerInputRate = 10),
    "13CPalm"     = list(myTracer = "13CPalm",     tracerInputRate = 32),
    "13COle"      = list(myTracer = "13COle",      tracerInputRate = 18),
    "13CLino"     = list(myTracer = "13CLino",     tracerInputRate = 36),
    "13CHB"       = list(myTracer = "13CHB",       tracerInputRate = 4)
  )
}  

# if (state_fasted.or.refed == "refed") {
#   # fed state
#   params <- list(
#     "13CGlc"      = list(myTracer = "13CGlc",      tracerInputRate = 80),
#     "13CLac"      = list(myTracer = "13CLac",      tracerInputRate = 49),
#     "13CGlycerol" = list(myTracer = "13CGlycerol", tracerInputRate = 10),
#     "13CAla"      = list(myTracer = "13CAla",      tracerInputRate = 30),
#     "13CGln"      = list(myTracer = "13CGln",      tracerInputRate = 10),
#     "13CPalm"     = list(myTracer = "13CPalm",     tracerInputRate = 16),
#     "13COle"      = list(myTracer = "13COle",      tracerInputRate = 36),
#     "13CLino"     = list(myTracer = "13CLino",     tracerInputRate = 28.8),
#     "13CHB"       = list(myTracer = "13CHB",       tracerInputRate = 6),
#     
#     # an "artificial tracer" to fit into the computation structure 
#     "13ChpAAb"    = list(myTracer = "13ChpAAb",    tracerInputRate = 0),
#     "13ChpAAc"    = list(myTracer = "13ChpAAc",    tracerInputRate = 0),
#     "13ChpAAd"    = list(myTracer = "13ChpAAd",    tracerInputRate = 0),
#     "13ChpAAe"    = list(myTracer = "13ChpAAe",    tracerInputRate = 0),
#     "13ChpAAf"    = list(myTracer = "13ChpAAf",    tracerInputRate = 0),
#     "13ChpAAj"    = list(myTracer = "13ChpAAj",    tracerInputRate = 0),
#     "13ChpAAk"    = list(myTracer = "13ChpAAk",    tracerInputRate = 0),
#     "13ChpAAl"    = list(myTracer = "13ChpAAl",    tracerInputRate = 0)
#   )
# }



# Perform parallel execution
l.stoich.EMU_allTracers <- parLapply(cl, params, function(p) {
  call_func_stoicEMU(myTracer = p$myTracer, tracerInputRate = p$tracerInputRate)
})

# Stop the cluster after computation
stopCluster(cl)





# make objects shared in all tracer experiments available directly from global environment for convenience
m.freeFlux_to_fullSetFlux <- l.stoich.EMU_allTracers $ '13CGlc' $ m.freeFlux_to_fullSetFlux
Stoich                    <- l.stoich.EMU_allTracers $ '13CGlc' $ Stoich
index.freeFluxes          <- l.stoich.EMU_allTracers $ '13CGlc' $ index.freeFluxes
n.EMU_sizes               <- l.stoich.EMU_allTracers $ '13CGlc' $ n.EMU_sizes
target.EMU                <- l.stoich.EMU_allTracers $ '13CGlc' $ target.EMU 
reaction_data             <- l.stoich.EMU_allTracers $ '13CGlc' $ reaction_data 
v.expressWithFreeFlux     <- l.stoich.EMU_allTracers $ '13CGlc' $ v.expressWithFreeFlux  


# write.csv(Stoich, file = "Stoich.csv")

# Calculate the total net CO2 produced in each organ (expressed as reaction index)
d.CO2 <- reaction_data %>% 
  # extract tissue information, excluding blood compartment
  mutate(reactions = reactions %>% str_remove(".Blood")) %>% 
  mutate(tissue = str_extract(reactions, "\\.[:alpha:]{1,3}") %>% str_remove(".")) %>% 
  filter(! is.na(tissue)) %>% 
  # keep rows containing the CO2 either as reactant or product
  filter(str_detect(reactions, "CO2")) %>% 
  mutate(coef = ifelse(str_detect(reactants, "CO2"), -1, 1)) %>% # 1 for producing CO2, and -1 for consuming CO2 for that reaction
  mutate(index.signed = enzyme * coef) %>% 
  
  # unify organelle-level naming to tissue level
  mutate(tissue = str_remove(tissue, "m$") %>% str_remove("c$")) %>% 
  arrange(tissue) 


# tissue index
i.CO2_Lv <- d.CO2 %>% filter(tissue == "Lv") %>% pull(index.signed)
i.CO2_M  <- d.CO2 %>% filter(tissue == "M")  %>% pull(index.signed)


# index of total CO2 to sink
i.CO2_sink <- reaction_data %>% filter(reactions == "CO2.Blood->sink") %>% pull(enzyme) 



# save the decomposed result separately for each state
if (state_fasted.or.refed == "fasted") {
  save.image(file = "3_list.stoich.EMU_allTracers_fasted.RData")
  
} else if (state_fasted.or.refed == "refed") {
  save.image(file = "3_list.stoich.EMU_allTracers_refed.RData")
}



beepr::beep(2) # mark completion
# system("say Mission Accomplished")
# beepr::beep(11)


time.EMU.End <- Sys.time(); time.EMU.End

sink("log_parallel decompose.txt")
time.EMU.End - time.EMU.Start
sink()

