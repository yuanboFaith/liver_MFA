rm(list=ls())

library(rebus)
library(tidyverse)
library(readxl)
library(accucor)

theme_set(
  theme_classic(base_size = 16) +
    theme(legend.position = "none", 
          # panel.border = element_rect(colour = "black", fill = NA, linewidth = .3),
          strip.background = element_blank(),
          axis.title.y = element_text(margin = margin(r = 10, unit = "pt")),
          axis.title = element_text(face = "bold", colour = "black"),
          axis.text = element_text(colour = "black"),
          strip.text = element_text(size = 12, face = "bold")) 
)


rstudioapi::getActiveDocumentContext()$path %>% dirname() %>% setwd(); getwd()






# $$-$$-$$-$$-$$-$$-$$-$$-$$-$$-$$-$$-$$-$$-$$-$$-$$-$$-$$-$$-$$-$$-$$-$$-$$-$$-$$-$$-$$-$$-$$-$$-$$-$$-$$-$$-$$-$$-$$-$$-$$-$$-$$-$$-$$-$$-$$-$$-

# FASTED PULSE

myPath <- "../Pulse Fasted/pulse data.xlsx"

# serum 
func.process <- function(whichSheet = "serum"){
  
  x <- (read_excel(myPath, sheet = whichSheet) %>% 
          natural_abundance_correction(resolution = 120000)) $ Normalized
  
  x %>% pivot_longer(-c(1, 2), names_to = "mouse", values_to = "labeling") %>% 
    separate(mouse, into = c("tissue", "mouse")) %>% 
    return()
}

x1 <- func.process(whichSheet = "serum");   x1 # batch 1 and 2
x2 <- func.process(whichSheet = "serum_2"); x2 # batch 3 (serum collected during the pulse stage of a targeted 'chase' experiment)

d.serum <- bind_rows(x1, x2) 


# liver TCA intermediates
d.liver.TCA <-  func.process(whichSheet = "liver TCA")

# liver TCA TAG C18：2
d.liver.TAG <- func.process(whichSheet = "liver TAG saponified")


# combine serum and liver data
d.labeling <- bind_rows(d.serum, d.liver.TCA, d.liver.TAG)


# add sample ID
d.id <- read_excel(myPath, sheet = "sampleID") %>% select(batch, mouse, `time (h)`)
d.labeling <- d.labeling %>% left_join(d.id)
d.labeling

# remove blood collected from the second batch (cardiac puncture, time delay vs tail blood during infusion)
d.labeling.Fast.Pulse <- d.labeling %>% filter(! (tissue == "serum" & batch == 2))

d.labeling.Fast.Pulse$`time (h)` %>% unique() %>% round(2) %>%  sort()
d.labeling.Fast.Pulse$tissue %>% unique()


d.labeling.Fast.Pulse.clean <- d.labeling.Fast.Pulse %>% 
  filter(`time (h)` != 0) %>%  # remove the t=0 mouse for plotting purpose
  filter(`time (h)` >= 2.5) %>%    # select later time points 2-12 hours 
  
  # keep only mice that has both serum and liver data
  group_by(mouse) %>%
  filter(any(tissue == "serum") & any(tissue == "LIVER")) %>%
  ungroup() %>% 
  
  # remove failed animal
  filter(mouse != "G") # serious leaking from the red pinport found at the end of infusion


d.labeling.Fast.Pulse.clean %>% select(Compound, tissue) %>% table()


# quick check via average labeling
d.labeling.Fast.Pulse.clean %>%
  group_by(mouse, Compound, tissue, `time (h)`) %>% 
  summarise(lab = sum(labeling * C_Label / max(C_Label))) %>% 
  ggplot(aes(x = `time (h)`, y = lab, color = mouse))+
  geom_point() + facet_grid(tissue~Compound) +
  # scale_y_log10() + annotation_logticks(sides = "l") +
  theme_bw()




d.labeling.Fast.Pulse.clean$mouse %>% unique()




# $$-$$-$$-$$-$$-$$-$$-$$-$$-$$-$$-$$-$$-$$-$$-$$-$$-$$-$$-$$-$$-$$-$$-$$-$$-$$-$$-$$-$$-$$-$$-$$-$$-$$-$$-$$-$$-$$-$$-$$-$$-$$-$$-$$-$$-$$-$$-$$-

# Fasted CHASE

myPath <- "chase data.xlsx"

# serum 
func.process <- function(whichSheet = "serum"){
  
  x <- (read_excel(myPath, sheet = whichSheet) %>% 
          natural_abundance_correction(resolution = 120000)) $ Normalized
  
  x %>% pivot_longer(-c(1, 2), names_to = "mouse", values_to = "labeling") %>% 
    return()
}

x1 <- func.process(whichSheet = "chase_serum_middle_point");   x1 
x2 <- func.process(whichSheet = "chase_serum_middle_point_2"); x2

x3 <- func.process(whichSheet = "chase_serum_SAC");            x3 
x4 <- func.process(whichSheet = "liver_Free_Metabolites_SAC"); x4
x5 <- func.process(whichSheet = "liver_TAG_SAC");              x5 


d.labeling.0 <- bind_rows(x1, x2, x3, x4, x5)

d.labeling.1 <- d.labeling.0 %>% 
  mutate(sample = mouse, 
         tissue = str_extract(sample, "[a-zA-Z]{1,10}"),
         timepoint = str_extract(sample, "(?<=_)[a-zA-Z]{1,10}"),
         mouse = str_extract(sample, "\\d$"))


# combine with sample ID
d.id <- read_excel(myPath, sheet = "sample.id") %>% select(sample, `time (h)`)
d.labeling.Fast.Chase <- d.labeling.1 %>% left_join(d.id)

# clean up
d.labeling.Fast.Chase$timepoint %>% unique()

d.labeling.Fast.Chase$tissue %>% unique()


# Clean up
d.labeling.Fast.Chase.clean <- d.labeling.Fast.Chase %>% 
  # keep only mice that has both serum and liver data
  group_by(mouse) %>%
  filter(any(tissue == "serum") & any(tissue == "liver")) %>%
  ungroup() 

d.labeling.Fast.Chase.clean %>% select(Compound, tissue) %>% table()
d.labeling.Fast.Chase.clean


# quick check via average labeling
d.labeling.Fast.Chase.clean %>%
  group_by(mouse, Compound, tissue, `time (h)`) %>% 
  summarise(lab = sum(labeling * C_Label / max(C_Label))) %>% 
  ggplot(aes(x = `time (h)`, y = lab, color = mouse))+
  geom_point() + facet_grid(tissue~Compound) +
  scale_y_log10() + annotation_logticks(sides = "l") +
  theme_bw()








# $$-$$-$$-$$-$$-$$-$$-$$-$$-$$-$$-$$-$$-$$-$$-$$-$$-$$-$$-$$-$$-$$-$$-$$-$$-$$-$$-$$-$$-$$-$$-$$-$$-$$-$$-$$-$$-$$-$$-$$-$$-$$-$$-$$-$$-$$-$$-$$-

# COMBINE FASTED AND REFED

x.Fasted.pulse <- d.labeling.Fast.Pulse.clean %>% select(-c(`time (h)`, batch)) # mouse id as big letters

x.Fasted.chase <- d.labeling.Fast.Chase.clean %>% 
  select(-c(timepoint, `time (h)`, sample)) %>% 
  relocate(Compound, C_Label, tissue, mouse, labeling) %>% 
  # convert mouse id to small letters for chase
  mutate(mouse = letters[as.integer(mouse)]) 

# ensure the mouse id of the two datasets are unique ()
x.Fasted.pulse$mouse %>% unique() %>% sort()
x.Fasted.chase$mouse %>% unique() %>% sort()

x.Fasted.chase <- x.Fasted.chase %>% 
  
  # there is a slight M+1 labeling not removed completely from natural abundance correction
  # due to the large n available, here this mouse is simply removed to make code faster
  filter(mouse != "e") %>% 
  
  # these two mice present data structure issue (error in MFA: one node produced an error: 'names' attribute [68] must be the same length as the vector [39])
  # due to the large n we have, here these two mice data are thrown away
  filter(! mouse %in% c("i", "h")) %>% 
  mutate(mouse = str_replace(mouse, "a", "O")) %>% 
  mutate(mouse = str_replace(mouse, "b", "P")) %>% 
  mutate(mouse = str_replace(mouse, "c", "Q")) %>% 
  mutate(mouse = str_replace(mouse, "d", "R")) %>% 
  mutate(mouse = str_replace(mouse, "f", "S")) %>% 
  mutate(mouse = str_replace(mouse, "g", "T")) 


d.Fast.kinetics <- bind_rows(x.Fasted.pulse, x.Fasted.chase)
# d.Fast.kinetics <- bind_rows(x.Fasted.pulse)



# # keep the most relevant compounds
# d.Fast.kinetics <- d.Fast.kinetics %>% 
#   filter(Compound %in% c("Malate", "Succinic acid", "TAG C18:2", "C18:2"))






# $$-$$-$$-$$-$$-$$-$$-$$-$$-$$-$$-$$-$$-$$-$$-$$-$$-$$-$$-$$-$$-$$-$$-$$-$$-$$-$$-$$-$$-$$-$$-$$-$$-$$-$$-$$-$$-$$-$$-$$-$$-$$-$$-$$-$$-$$-$$-$$-

# Organize data to fit for MFA
d.Fast.kinetics2 <- d.Fast.kinetics %>% 
  mutate(tissue = str_replace(tissue, "serum", "Blood")) %>% 
  mutate(tissue = str_replace(tissue, "liver", "Lv")) %>% 
  mutate(tissue = str_replace(tissue, "LIVER", "Lv"))


# simplify the names of compounds
d.metab_names <- tribble(
  ~Compound,               ~abbr,   
  "Glucose",               "Glc",   
  "3-Hydroxybutyric acid", "HB",    
  "Lactate",               "Lac",   
  "Malate",                "Mal",   
  "Succinic acid",         "Suc",   
  "Alanine",               "Ala",   
  "Glutamine",             "Gln",
  "C18:2",                 "Lino",
  "TAG C18:2",             "TAGLino")

CC <- d.metab_names$Compound
aa <- d.metab_names$abbr
names(aa) <- CC; aa


sources <- c("Glucose", "Lactate", "Alanine", "Glutamine", "C18:2", "3-Hydroxybutyric acid")

d.dataForMFA.TAGkinetics.Fasted <- d.Fast.kinetics2 %>% 
  
  # select relevant compounds; e.g., lactate appear in both liver and serum, but use the serum data
  filter( (tissue == "Lv" & Compound %in% c("Malate", "Succinic acid", "TAG C18:2")) | 
            tissue == "Blood" & Compound %in% sources) %>% 
  
  # rename compounds using the unified abbreviations
  mutate(Compound = aa[Compound]) %>% 
  
  # rename mouse id
  rename(m.id = mouse) %>% 
  
  # standard deviation
  group_by(tissue, Compound, C_Label) %>%  
  mutate(labeling.sd = sd(labeling)) %>% ungroup() %>% 
  mutate(labeling.sd =  replace_na(labeling.sd, 0)) %>% 
  
  mutate(State = "fasted",
         Infusate = paste0("13CLinoKin", m.id)) %>% 
  mutate(mouse.when.who = paste0(m.id, "_2025Sep_BY")) %>% 
  mutate(infuse.nmol.min.g = NA) %>% 
  mutate(Compound.tissue          = paste0(Compound, ".", tissue)) %>% 
  
  

    # Simplify TAG-linoleate and blood linoleate's M+18 as M+2 (as assembled repeating acetyl-CoA units)
  filter(! (Compound %in% c("Lino", "TAGLino" ) & C_Label %in% 2:17 )) %>%  # remove M+2 - M+17
  mutate(C_Label = ifelse(Compound %in% c("Lino", "TAGLino" ) & C_Label == 18, 2, C_Label)) %>% 
  
  # calculate C -max; manually adjust SCFA C max (not considering derivatization)
  group_by(Compound) %>% mutate(C.max = max(C_Label)) %>%  # C max
  
  # rowwise important to add the correct C-max
  rowwise() %>% 
  mutate(Compound.tissue.seq      = str_c(Compound.tissue, "_", str_c(1:C.max, collapse = ""))) %>% ungroup() %>% 
  mutate(Compound.tissue.seq_m.id = str_c(Compound.tissue.seq, "|", m.id)) %>% 
  
  relocate(Compound, C_Label, tissue, labeling, labeling.sd, State, Infusate, mouse.when.who, 
           infuse.nmol.min.g, Compound.tissue, m.id, C.max, Compound.tissue.seq, Compound.tissue.seq_m.id) %>% 
  
  
  # each mouse need to have the Blood Lino and Liver TAGLino
  group_by(m.id) %>% 
  filter(any(Compound == "Lino") & any(Compound == "TAGLino")) %>% 
  ungroup()


d.dataForMFA.TAGkinetics.Fasted$m.id %>% unique()
d.dataForMFA.TAGkinetics.Fasted$Compound %>% unique()


  



# a quick visual check
d.dataForMFA.TAGkinetics.Fasted %>% 
  filter(Compound %in% c("Lino", "TAGLino" )) %>% 
  filter(C_Label != 0) %>% 
  ggplot(aes(x = Compound, y = labeling, fill = as.character(C_Label))) +
  geom_col() + facet_wrap(~m.id) + theme(legend.position = "right")




save(d.dataForMFA.TAGkinetics.Fasted, 
     file = "../../data/data_section_Fasted_TAGkinetics.RData") 


