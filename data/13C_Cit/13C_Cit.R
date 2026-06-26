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
      # Note: citrate in blood signal can be quite low, so here we have to use lower cutoff to keep citrate labeling in blood from being out
      if ( ij[1] < 1*(10^5) ) d.area[d.area$compound == i, ][[j]] <- rep(NA, length(ij))
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
d.13C_Cit <- tibble(.rows = 0)

for (a in dataFiles) {
  # check a=dataFiles[1]
  paste("working on file:", a) %>% print()
  d.13C_Cit <- bind_rows( d.13C_Cit, func.cleanUp(myExcel = a) )
}

head(d.13C_Cit)
tail(d.13C_Cit)
dim(d.13C_Cit)

d.13C_Cit$tissue %>% unique()


# Check with Tony's summary data
d.13C_Cit %>% 
  select(State, Infusate, mouse, `Concentration (mM)`, `Rate (µl/min/g)`, When, Who) %>% 
  distinct() %>% 
  arrange(State, When)
# GOOD!


# CHECK DATA CONTENT
d.13C_Cit$tissue %>% unique()
d.13C_Cit$Compound %>% unique()

# check data of serum: keep the samples from the last time point 
d.13C_Cit %>% 
  filter(str_detect(tissue, "erum") | str_detect(tissue, "srm")) %>% 
  filter(C_Label != 0) %>% 
  filter(Compound  %in% c("glutamine", "glucose", "lactate", "citrate")) %>%
  ggplot(aes(x = tissue, y = labeling, color = factor(C_Label))) +
  geom_point() +
  geom_line(aes(group = C_Label)) +
  facet_grid(Compound~mouse.when.who, scales = "free") +
  theme(axis.text.x = element_text(angle = 50, hjust = 1)) 

d.13C_Cit <- d.13C_Cit %>% 
  filter(! (mouse.when.who == "M1_2018-01-09_CB" & tissue %in% c("Serum1"))) %>% 
  filter(! (mouse.when.who == "M1_2018-12-21_CB" & tissue %in% c("srm120min"))) %>% 
  filter(! (mouse.when.who == "M1_2019-01-17_CB" & tissue %in% c("serum140min"))) %>% 
  filter(! (mouse.when.who == "M2_2018-01-09_CB" & tissue %in% c("Serum1"))) %>% 
  filter(! (mouse.when.who == "M2_2018-12-21_CB" & tissue %in% c("srm115min"))) %>% 
  mutate(tissue = str_remove(tissue, "min$") %>% str_remove("\\d{1,3}$")) # remove ending "min", and the ending 1 to 3 digits



# check duplication : found two mice containing duplicates / replicates in tissue "srm"
# for each state, infusate, mouse, tissue, compound, C_Label, there should be only one measurement
x <- d.13C_Cit %>% 
  select(State, mouse.when.who, Compound, tissue, C_Label) %>% 
  duplicated()

sum(x)



# Check kidney data number of replicates
filter(d.13C_Cit, State == "fasted")$tissue %>% unique()
d.13C_Cit %>% filter(State == "fasted") %>%  filter(tissue %in% c("Kd") ) %>% pull(mouse.when.who) %>% unique()
d.13C_Cit %>% filter(State == "fasted") %>%  filter(tissue %in% c("kd") ) %>% pull(mouse.when.who) %>% unique()



# Check kidney's M+6 citrate value, without one mouse' data clearly wrong
d.13C_Cit %>% 
  # filter(State == "fasted") %>% 
  filter(tissue %in% c("Kd", "kd", "kid")) %>% 
  filter(C_Label != 0) %>% 
  filter(Compound == "citrate") %>% 
  # filter(Compound == "succinate") %>% 
  ggplot(aes(x = C_Label, y = labeling, color = mouse.when.who)) +
  geom_point() +
  geom_line(aes(group = mouse.when.who)) +
  facet_wrap(~State) +
  scale_y_continuous(breaks = seq(0, 1, .1))
  
# remove one mouse's kidney citrate data, as the M+6 labeling is clearly wrong (peak value being zero)
d.13C_Cit <- d.13C_Cit %>% filter(! (mouse.when.who == 'M1_2018-12-21_CB' & tissue == "kd" & Compound == "citrate"))



save(d.13C_Cit, file = "13C_Cit.RData")
