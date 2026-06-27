# This script creates functions:
# (1) Pseudo-inverse of a matrix
# (2) Turn matrix rows into elements of a list
# (3) Replace flux notations with actual flux values
# This script is not needed for "Stoic_EMU_Calculator.R", and not needed for '13C_tracing_labeling_data.R'.
# This script is required for the iterative optimization procedure.

rm(list = ls())
library(tidyverse)

# set path to the folder of the current active script
rstudioapi::getActiveDocumentContext()$path %>% dirname() %>% setwd(); getwd()


# (1) Calculate numeric stable inverse
func.pseudo_inverse_svd <- function(A) {
  
  if(nrow(A) == 1 & ncol(A) == 1) {
    return(1/A)
  } else {
    svd_res <- svd(A)
    U <- svd_res$u
    D <- svd_res$d
    V <- svd_res$v
    
    # Tolerance for singular values (machines precision)
    tol <- max(dim(A)) * max(D) * .Machine$double.eps
    
    # Invert only non-zero singular values
    D_inv <- diag(ifelse(D > tol, 1 / D, 0))
    
    # Calculate pseudo-inverse
    A_pseudo_inv <- V %*% D_inv %*% t(U)
    return(A_pseudo_inv)
  }
}
# Check
matrix(c(1, 2, 3, 4, 0, 3, 100, 2330), nrow = 4, ncol = 2) %>% func.pseudo_inverse_svd()




# (2) Turn a matrix array in a list, with each row being an element, row names as the corresponding element names
func.turnMatrixToList <- function(mat){
  l <- list()
  for (i in seq(nrow(mat))) {
    l[i] = mat[i, ] %>% as.data.frame()
  }
  names(l) <- rownames(mat)
  return(l)
}
# Check
mtcars[1:3, 1:3] %>% as.matrix() %>% func.turnMatrixToList()




# (3) Replace a single flux notation with tangible flux value; 
# This function is original defined in `if (optim.step == 1){...` in the iterative optimization script.
# !!! It requires object 'u' in the environment.
# Here the code is able to handle non-unit coefficient
# pattern matching is critical; e.g., 'v1' should not match 'v10'
func.FillwithFluxValue <- function(x, whichFlux = "v10"){
  str_replace_all(x, pattern = paste0("\\b", whichFlux, "\\b"), replacement = u[whichFlux] %>% as.character())
}
# check
# "-(v10+v4+v6)" %>% func.FillwithFluxValue(whichFlux = "v1")





save.image(file = "1_Supplement_functions.RData")