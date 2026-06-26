library(tidyverse)
library(readxl)

rm(list = ls())

# default path to the current folder
rstudioapi::getActiveDocumentContext()$path %>% dirname() %>% setwd()
getwd()

path_suppl.Data <- "/Users/boyuan/Desktop/Harvard/Research/Energy Metabolism Atlas/data/suppl_data.xlsx"
d.formula      <- read_excel(path_suppl.Data, sheet = "formula")
d.tissueNaming <- read_excel(path_suppl.Data, sheet = "tissueNames")


func.cleanUp <- function(myExcel) {
  #check myExcel="20190427_13C_glc_refed_glc1_glc2_tissue.xlsm"
  
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
      if ( ij[1] < 1e5 ) d.area[d.area$compound == i, ][[j]] <- rep(NA, length(ij))
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



# CHECK ALL DATA FILE COMPOUND NAMES, IF PRESENT IN THE FORMULA SHEET
dataFiles <- list.files(pattern = "\\.xlsm$")

compouds.All <- c()

for (a in dataFiles) {
  # check a=dataFiles[1]
  paste("working on file:", a) %>% print()
  compouds.All <- c(compouds.All, read_excel(a)$compound %>% unique()) %>% unique()
}

# print which compounds are missing from the formula sheet
compouds.All[!compouds.All %in% d.formula$compound] 




# COMPILE ALL DATA FILES
d.13C_Palm <- tibble(.rows = 0)

for (a in dataFiles) {
  # check a=dataFiles[1]
  paste("working on file:", a) %>% print()
  d.13C_Palm <- bind_rows( d.13C_Palm, func.cleanUp(myExcel = a) )
}

head(d.13C_Palm)
tail(d.13C_Palm)
dim(d.13C_Palm)



# Check with Tony's summary data
d.13C_Palm %>% 
  select(State, Infusate, mouse, `Concentration (mM)`, `Rate (µl/min/g)`, When, Who) %>% 
  distinct() %>% 
  arrange(State, When)



# CHECK DATA CONTENT
d.13C_Palm$tissue %>% unique()
d.13C_Palm$Compound %>% unique()

# keep data of serum (steady state) collected from the last time points (more approaching steady state)
d.13C_Palm %>% filter(tissue  %>% str_detect("SrTC")) %>% 
  filter(C_Label != 0) %>% 
  filter(Compound  %in% c("glutamine", "glucose", "lactate", "alanine",  "palmitic acid")) %>%
  ggplot(aes(x = tissue, y = labeling, color = factor(C_Label))) +
  geom_point() +
  geom_line(aes(group = C_Label)) +
  facet_wrap(Compound~mouse.when.who, scales = "free_y")

# glucose labeling differ by two times between SrTC1 and SrTC2. Lac, Gln, and Palm labeling is consistent. 
# As typically done, choose the second time point
d.13C_Palm <- d.13C_Palm %>% 
  filter(! (mouse.when.who %in% c("M3_2019-01-17_SH", "M4_2019-01-17_SH") & tissue == "SrTC1")) %>% 
  mutate(tissue = str_remove(tissue, "\\d$"))

d.13C_Palm



# check duplication
# for each state, infusate, mouse, tissue, compound, C_Label, there should be only one measurement
x <- d.13C_Palm %>% 
  select(State, mouse.when.who, Compound, tissue, C_Label) %>% 
  duplicated()

sum(x)


save(d.13C_Palm, file = "13C_Palm.RData")
