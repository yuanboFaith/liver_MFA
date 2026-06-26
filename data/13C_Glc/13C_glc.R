library(tidyverse)
library(readxl)

rm(list = ls())

# default path to the current folder
rstudioapi::getActiveDocumentContext()$path %>% dirname() %>% setwd()

path_suppl.Data <- "/Users/boyuan/Desktop/Harvard/Research/Energy Metabolism Atlas/data/suppl_data.xlsx"
d.formula      <- read_excel(path_suppl.Data, sheet = "formula")
d.tissueNaming <- read_excel(path_suppl.Data, sheet = "tissueNames")


func.cleanUp <- function(myExcel) {
  #check myExcel="20190222_13C_glucose_refed_M2_M3_serum.xlsm"
  
  # MOUSE ID
  d.ID <- read_excel(path = myExcel, sheet = "ID") %>% 
    select(1:7) %>% rename(mouse = Mouse) %>% relocate(When, Who, .after = last_col()) %>% mutate(When = ymd(When))
  d.ID
  
  # LABELING DATA
  d.area <- read_excel(path = myExcel, sheet = "data")
  d.area
  
  
  # BIG AREA. If a compound's parent peak is too small in a smaple, its labeling would not be reliable; discard the data of all isotopologues
  for ( i in unique(d.area$compound) ){ # loop through compounds
    for (j in 3:ncol(d.area) ) { # loop through samples
      # i="glycine"
      # j="M2_SrTC"
      ij <- filter(d.area, compound == i)[[j]] # vector of areas of the given compound in the given sample
      # if parent peak (first value) is small, e.g., smaller than 1e5, then the labeling is not reliable enough, turn all values to NA
      if ( ij[1] < 1*10^5 ) d.area[d.area$compound == i, ][[j]] <- rep(NA, length(ij))
    }
  }
  
  
  # NATURAL ABUNDANCE CORRECTION
  l.natural_abundance_corrected <- d.area %>% 
    left_join(d.formula, by = "compound") %>% relocate(note, compound, formula) %>% # add formula
    rename(isotopelabel = note) %>%  # do natural abundance correction
    accucor::natural_abundance_correction(resolution = 120000)
  
  
  # TIDY UP
  d.norm <- l.natural_abundance_corrected$Normalized %>% 
    pivot_longer(-c(1:2), names_to = "sample", values_to = "labeling") %>% 
    filter( ! is.na(labeling) ) %>% 
    arrange(Compound, sample, C_Label) %>% 
    separate(sample, into = c("mouse", "tissue"))
  
  # COMBINE WITH ID DATA
  
  # check the mouse ID matches; otherwise report error
  if(setequal( unique(d.norm$mouse), d.ID$mouse) == F) {
    stop(paste("in Excel:", myExcel, ", mouse ID mismatches."))
  }
  
  d.norm.ID <- d.norm %>% left_join(d.ID, by = "mouse") %>% 
    mutate(mouse.when.who = str_c(mouse,"_", When, "_", Who)) # unique mouse ID
  d.norm.ID
  
}

func.cleanUp(myExcel = "20160624_13C_Glc_Refed_M2.xlsm")





# CHECK ALL DATA FILE COMPOUND NAMES, IF PRESENT IN THE FORMULA SHEET
dataFiles <- list.files(pattern = "\\.xlsm$")

compouds.All <- c()

for (a in dataFiles) {
  # check a=dataFiles[1]
  paste("working on file:", a) %>% print()
  compouds.All <- c(compouds.All, read_excel(a)$compound %>% unique()) %>% unique()
}

compouds.All[!compouds.All %in% d.formula$compound] # print which compounds are missing from the formula sheet





# COMPILE ALL DATA FILES
d.13C_glc <- tibble(.rows = 0)

for (a in dataFiles) {
  # check a=dataFiles[1]
  paste("working on file:", a) %>% print()
  d.13C_glc <- bind_rows( d.13C_glc, func.cleanUp(myExcel = a) )
}

head(d.13C_glc)
tail(d.13C_glc)
dim(d.13C_glc)



# CHECK DATA CONTENT
# 1) remove unneeded samples
# In many excels, the tissue is measured twice, and here we keep only one replicate. 
# As the number is many, here we process the duplication in code. 

d.13C_glc$tissue %>% unique() 

d.13C_glc2 <- d.13C_glc %>% 
  # remove 'tissue2', and remove the digit in 'tissue1'
  filter(!str_detect(tissue, "2$")) %>% 
  mutate(tissue = str_remove(tissue, "1$"))
  
d.13C_glc2$tissue %>% unique() 



# Check with Tony's summary data
d.13C_glc2 %>% 
  select(State, Infusate, mouse, `Concentration (mM)`, `Rate (µl/min/g)`, When, Who) %>% 
  distinct() %>% 
  arrange(State, When)


# check duplication
# for each state, infusate, mouse, tissue, compound, C_Label, there should be only one measurement
x <- d.13C_glc2 %>% 
  select(State, mouse.when.who, Compound, tissue, C_Label) %>% 
  duplicated()

sum(x)

save(func.cleanUp, d.13C_glc2, file = "13C_Glc.RData")


d.13C_glc2 %>% filter(mouse.when.who == "M1_2016-06-30_SH") %>% 
  select(State, Infusate, mouse.when.who, tissue, Compound, C_Label) %>% duplicated() %>% sum()

