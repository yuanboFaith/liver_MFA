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
      if ( ij[1] < 1e6 ) d.area[d.area$compound == i, ][[j]] <- rep(NA, length(ij))
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

func.cleanUp(myExcel = "20160513_13C_Lac_Fasted_M3_serum.xlsm")




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
d.13C_lac <- tibble(.rows = 0)

for (a in dataFiles) {
  # check a=dataFiles[1]
  paste("working on file:", a) %>% print()
  d.13C_lac <- bind_rows( d.13C_lac, func.cleanUp(myExcel = a) )
}

head(d.13C_lac)
tail(d.13C_lac)
dim(d.13C_lac)



# CHECK DATA CONTENT
# For duplicated measurement (not too many), they have been manually processed in each individual excel. 
d.13C_lac$tissue %>% unique()



func.plt <- function(myData){
  myData %>% 
    filter(Compound %in% c("glucose", "malate",  "succinate", "alanine")) %>% 
    filter(C_Label != 0) %>% 
    ggplot(aes(x = Compound, y = labeling, color = factor(C_Label))) +
    geom_point() +
    facet_wrap(~tissue) +
    theme(axis.text.x = element_text(angle = 60, hjust = 1))
}

d.13C_lac$tissue         %>% unique()
d.13C_lac$mouse.when.who %>% unique()

# all mice ID
x <- (d.13C_lac %>% filter(tissue %>% str_detect("[12]$")))$mouse.when.who  %>% unique()

d.13C_lac %>% filter(mouse.when.who == x[1]) %>% func.plt()
d.13C_lac <- d.13C_lac %>% mutate(tissue = ifelse(mouse.when.who == x[1], str_remove(tissue, "[12]$"), tissue))


d.13C_lac %>% filter(mouse.when.who == x[2]) %>% func.plt()
d.13C_lac <- d.13C_lac %>% mutate(tissue = ifelse(mouse.when.who == x[2], str_remove(tissue, "[12]$"), tissue))


d.13C_lac %>% filter(mouse.when.who == x[3]) %>% func.plt()
d.13C_lac <- d.13C_lac %>% mutate(tissue = ifelse(mouse.when.who == x[3], str_remove(tissue, "[12]$"), tissue))


d.13C_lac %>% filter(mouse.when.who == x[4]) %>% func.plt()
# most tissues are duplicated, with very similar result. Br2 has gluc data, but not in Br1. SrTC1 and SrTC2 very similar result. 
# keep tissue-2 for all tissues; keep TA1 and BAT1 - they don't have replicates
d.13C_lac <- d.13C_lac %>% 
  filter(! (mouse.when.who == x[4] & tissue %in% c("Br1", "Ht1", "Kd1", "Lv1", "Pc1", "Q1", "SI1", "Sp1", "SrTC1") )) %>% 
  mutate(tissue = ifelse(mouse.when.who == x[4], str_remove(tissue, "[12]$"), tissue))


d.13C_lac %>% filter(mouse.when.who == x[5]) %>% func.plt()
d.13C_lac <- d.13C_lac %>% 
  mutate(tissue = ifelse(mouse.when.who == x[5], str_remove(tissue, "[2]$"), tissue)) # only serum has digit suffix


d.13C_lac %>% filter(mouse.when.who == x[6]) %>% func.plt()
# for duplicates: Ht1 and Ht2, Q1 and Q2 - the result is similar; keep one by removing -2
d.13C_lac <- d.13C_lac %>% 
  filter(!(mouse.when.who == x[6] & tissue %in% c("Ht2", "Q2"))) %>% 
  mutate(tissue = ifelse(mouse.when.who == x[6], str_remove(tissue, "[12]$"), tissue))


d.13C_lac %>% filter(mouse.when.who == x[7]) %>% func.plt()
# for duplicates: Ht1 and Ht2, Q1 and Q2 - the result is similar; keep one by removing -2
d.13C_lac <- d.13C_lac %>% 
  filter(!(mouse.when.who == x[7] & tissue %in% c("Ht2", "Q2"))) %>% 
  mutate(tissue = ifelse(mouse.when.who == x[7], str_remove(tissue, "[12]$"), tissue))


d.13C_lac %>% filter(mouse.when.who == x[8]) %>% func.plt()
# for duplicates: Ht1 and Ht2, Q1 and Q2 - the result is similar; keep one by removing -2
d.13C_lac <- d.13C_lac %>% 
  filter(!(mouse.when.who == x[8] & tissue %in% c("Ht2", "Q2"))) %>% 
  mutate(tissue = ifelse(mouse.when.who == x[8], str_remove(tissue, "[12]$"), tissue))


d.13C_lac %>% filter(mouse.when.who == x[9]) %>% func.plt()
# for duplicates: Ht1 and Ht2, Q1 and Q2 - the result is similar; keep one by removing -2
d.13C_lac <- d.13C_lac %>% 
  filter(!(mouse.when.who == x[9] & tissue %in% c("Ht2", "Q2"))) %>% 
  mutate(tissue = ifelse(mouse.when.who == x[9], str_remove(tissue, "[12]$"), tissue))

d.13C_lac$tissue %>% unique()


# check duplication
x <- d.13C_lac %>% 
  select(State, mouse.when.who, Compound, tissue, C_Label) %>% 
  duplicated()

sum(x)
d.13C_lac[x, ]$mouse.when.who %>% unique()


save(d.13C_lac, file = "13C_lac.RData")



# Check with Tony's summary data
d.13C_lac %>% 
  select(State, Infusate, mouse, `Concentration (mM)`, `Rate (µl/min/g)`, When, Who) %>% 
  distinct() %>% 
  arrange(State, When)


