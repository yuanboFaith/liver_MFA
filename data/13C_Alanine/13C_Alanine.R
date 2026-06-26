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
  #check myExcel="20181125_13C_alanine_fasted_M1_M2_serum.xlsm"
  
  # MOUSE ID
  d.ID <- read_excel(path = myExcel, sheet = "ID") %>% 
    select(1:7) %>% rename(mouse = Mouse) %>% relocate(When, Who, .after = last_col()) %>% mutate(When = ymd(When))
  d.ID
  
  # LABELING DATA
  d.area <- read_excel(path = myExcel, sheet = "data")
  d.area
  
  
  # BIG AREA. If a compound's parent peak is too small in a smaple, its labeling would not be reliable; discard the data of all isotopologues
  for ( i in unique(d.area$compound) ){ # loop through compounds
    # print(paste("check compound: ", i))
    
    for (j in 3:ncol(d.area) ) { # loop through samples
      # print(paste("check column (sample): ", j))
      
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
d.13C_Alanine <- tibble(.rows = 0)

for (a in dataFiles) {
  # check a=dataFiles[1]
  paste("working on file:", a) %>% print()
  d.13C_Alanine <- bind_rows( d.13C_Alanine, func.cleanUp(myExcel = a) )
}

head(d.13C_Alanine)
tail(d.13C_Alanine)
dim(d.13C_Alanine)



# Check with Tony's summary data
d.13C_Alanine %>% 
  select(State, Infusate, mouse, `Concentration (mM)`, `Rate (µl/min/g)`, When, Who) %>% 
  distinct() %>% 
  arrange(State, When)



# CHECK DATA CONTENT
d.13C_Alanine$tissue %>% unique()
d.13C_Alanine$Compound %>% unique()

# Check serum data Serum0 to Serum4
d.13C_Alanine %>% filter(tissue %in% paste0("Serum", 0:4)) %>% 
  filter(C_Label != 0) %>% 
  filter(Compound  %in% c("glutamine", "glucose", "lactate", "alanine")) %>%
  ggplot(aes(x = tissue, y = labeling, color = factor(C_Label))) +
  geom_point() +
  geom_line(aes(group = C_Label)) +
  facet_grid(Compound~mouse.when.who, scales = "free_y") +
  theme(axis.text.x = element_text(angle = 50, hjust = 1))

# for ala1_2018-11-25_TT, ala2_2018-11-25_TT, keep the last time point t4
# for ala1_2018_12-03_tt, ala3_2018-12-30_TT, keep the last time point t2
d.13C_Alanine <- d.13C_Alanine %>% 
  filter(! (mouse.when.who %in% c("ala1_2018-12-03_TT", "ala3_2018-12-30_TT") & tissue %in% paste0("Serum", 0:1))) %>% 
  filter(! (mouse.when.who %in% c("ala1_2018-11-25_TT", "ala2_2018-11-25_TT") & tissue %in% paste0("Serum", 0:3)))
# now we get Serum2 and Serum4 left in the dataset, we'll remove the suffix later


# Check Serumt0 and Serumt1
d.13C_Alanine %>% filter(tissue %in% c("Serumt0", "Serumt1")) %>% 
  filter(C_Label != 0) %>% 
  filter(Compound  %in% c("glutamine", "glucose", "lactate", "alanine")) %>%
  ggplot(aes(x = tissue, y = labeling, color = factor(C_Label))) +
  geom_point() +
  geom_line(aes(group = C_Label)) +
  facet_grid(Compound~mouse.when.who, scales = "free_y")

# keep Serumt1; t0 is baseline
d.13C_Alanine <- d.13C_Alanine %>% 
  filter(! (mouse.when.who %in% c("ala2_2018-12-03_TT") & tissue == "Serumt0" )) %>% 
  mutate(tissue = ifelse(tissue == "Serumt1", "Serum", tissue))


# remove digit suffix in tissue names
d.13C_Alanine <- d.13C_Alanine %>% mutate(tissue = str_remove(tissue, "\\d$"))

d.13C_Alanine$tissue %>% unique()



# check duplication
# for each state, infusate, mouse, tissue, compound, C_Label, there should be only one measurement
x <- d.13C_Alanine %>% 
  select(State, mouse.when.who, Compound, tissue, C_Label) %>% 
  duplicated()

sum(x)


save(d.13C_Alanine, file = "13C_Alanine.RData")




