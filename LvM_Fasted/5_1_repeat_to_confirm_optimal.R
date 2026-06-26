# This script repeats the MFA procedure multiple times
# and select optimal fluxes from the one with the lowest cost

###############################################################################################################################
#  To run this current script, in the file '4_iteration_convergence_parallel':
#  'rm(list = ls())' needs to be COMMENTED OUT 
#  as the loop index 'repeat_index' and 'time.Start.MFA.repeat' needs to be present in the global environment
###############################################################################################################################

library(tidyverse)
rm(list = ls())

time.Start.MFA.repeat <- Sys.time()

# set path to the folder of the current active script
rstudioapi::getActiveDocumentContext()$path %>% dirname() %>% setwd(); getwd()



# # # !!!!!! CAUTION !!!!!!
# 1. Clear the 'repeat_results' folder
#  .Rdata files
list.files(path = "repeat_results", pattern = "\\.RData$", full.names = TRUE) %>% file.remove() 
# .txt log files
list.files(path = "repeat_results", pattern = "\\.txt$", full.names = TRUE) %>% file.remove() 

# 2. Clear all PDF and txt files in the 'plots' folder
pdf_files <- list.files(
  path       =  file.path("plots"),
  # pattern    = "\\.pdf$",
  recursive  = TRUE,   # search subfolders
  full.names = TRUE)

# Remove the PDF files
file.remove(pdf_files)







for (repeat_index in 1:15) {
  # repeat_index=1
  # calculate total CO2 per organ
  message("\nrepeat MFA core procedure to confirm optimal: repeat #", repeat_index )
  
  attempt <- 1
  max_attempts <- 3
  
  repeat {
    tryCatch(
      {
        
        objs_before <- ls() # later will remove all objects generated from '4_MFA_convergence.R'
        
        source("4_MFA_convergence.R")
        
        # rm(l.stoich.EMU_allTracers) # remove the decomposition result to save space
        save(
          cost.iterations, u.initial, d.u.iterations,
          optimal_cost, optimal_free_flux, optimal_allFlux, 
          d.v.final, l.EMU.sim, d.obs.simu, d.flux.CO2, # d.flux.CO2.summary, d.EMU_exchCO2.all,
          file = paste0("./repeat_results/repeat_", repeat_index, ".RData")
        )
        
        # remove the many variables from the global environment, 
        # otherwise they'll be exported to the environment in the many working nodes in the next round
        
        # option 1 - manually specify the objects
        # to_keep <- c("time.Start.MFA.repeat", "repeat_index")
        # rm(list = setdiff(ls(), to_keep))
        
        # option 2 - high efficiency way
        rm(list = setdiff(ls(), objs_before))
        
        # If successful, break out of the repeat loop
        break
      },
      error = function(e) {
        message("\n\nThis repeat was not successful: ", e$message)
        if (attempt >= max_attempts) {
          message(max_attempts, " attempts on MFA repeat have been tried, but all failed!!! OMG!!! What the odds!!!")
          # Return NULL or an appropriate value to indicate failure
          break
        } else {
          message("\n\nRetrying... (Attempt: ", attempt + 1, ")")
        }
      }
    )
    attempt <- attempt + 1
  }
  
  # clear environment, except the starting time, and the repeat index
  # to_keep <- c("time.Start.MFA.repeat", "repeat_index")
  # rm(list = setdiff(ls(), to_keep)) 
}


# record time 
time.End.MFA.repeat <- Sys.time()

# after running the MFA core iteration source script #4, the path is the 'MFA_core' folder
# in the output path, need to change it to the 'repeat_optimal' folder - 'repeat_results' subfolder
sink("./repeat_results/log_repeat_to_confirm_optimal.txt")
time.End.MFA.repeat - time.Start.MFA.repeat
sink()

