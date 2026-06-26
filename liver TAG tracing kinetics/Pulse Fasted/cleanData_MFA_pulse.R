rm(list=ls())

library(tidyverse)

rstudioapi::getActiveDocumentContext()$path %>% dirname() %>% setwd(); getwd()
load("liver_kinetics_pulse.RData")



#  OPTIONAL: quickly check Tony's malate labeling with Lino infusion
load("../../data/cleaned_labeling_data.RData")
d.13C.fasted %>% filter(Infusate == "13CLino" & tissue == "Lv" & Compound == "Mal") %>% 
  group_by(mouse.when.who) %>% 
  summarise(L.mal = sum(C_Label / 4 * labeling))
# so Tony's and my infusion rate is the same, but Tony's serum lino labeling is almost twice lower...lipolysis flux may not be determined consistently. 
# What is great though is that Tony's liver malate labeling is also two times lower than mine. Good! For liver this is very consistent regarding hepatic use of linoleate




C.old <- d.labeling$Compound %>% unique() %>% sort(); C.old
C.new <- c("HB", "Ala", "aKG", "Palm", "Ole", "Lino", "Glc", "Gln", "Lac", "Mal", "Suc", "TAGlino")
names(C.new) <- C.old

T.old <- d.labeling$tissue %>% unique() %>% sort(); T.old
T.new <- c("Lv", "Blood")
names(T.new) <- T.old

# these names should be used consistently across MFA modeling
# "Cit"      "Mal"      "Suc"      "Ala"      "CO2"      "Glc"      "Gln"      "HB"       "Lac"      "aKG"      "Lino"     "Ole"      "Palm"     "Glycerol"
# "Bat"   "Blood" "H"     "I"     "K"     "Ln"    "Lv"    "M"     "P"     "S"     "gW"    "iW"    "Br"   

d.labeling2 <- d.labeling %>% 
  filter(mouse %in% c("A", "B", "C", "D", "E", "F", "H", "I")) %>% 
  mutate(Compound = C.new[ Compound ]) %>% 
  mutate(tissue   = T.new[ tissue   ]) %>% 
  mutate(State = "fasted", 
         Infusate = "13CLino", 
         mouse.when.who = str_c(mouse, "_", "2025-8-X_BY"),
         infuse.nmol.min.g = 36 * 3 / 27, # 36 mM at 3 uL/min/animal; in equivalent C2 units; 27 g BW
         Compound.tissue = paste0(Compound, ".", tissue)) %>% 
  rename(m.id = mouse) 




# Fatty acids are barely labeled unless infused; only labeling data of FAs in blood is used
# turn fatty acids into C2 equivalents: consider M+16 / m+18 as M+2, M+1 as 0, and M+0 is unchanged
# plot fatty acids
f.plt_FAs <- function(myData, myState){
  myData %>% 
    filter(State == myState) %>% 
    filter(Infusate %in% c("13CLino")) %>% 
    filter(Compound %in% c("Lino", "Palm", "Ole")) %>% 
    filter(tissue == "Blood") %>% 
    ggplot(aes(x = mouse.when.who , y = labeling, fill = factor(C_Label))) +
    geom_col(position = "stack", color = "black") +
    geom_text(aes(label = C_Label), position = position_stack(vjust = .5), size = 3, color = "black") +
    scale_fill_brewer(palette = "Set2") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_y_continuous(breaks = seq(0, 1, .05), expand = expansion(mult = c(0, 0.02))) +
    facet_wrap(~Compound)
}
f.plt_FAs(myData = d.labeling2, myState = "fasted")



f.FAs_2C <- function(whichFA, C){
  d.labeling2 %>% 
    filter(Compound == c(whichFA)) %>% 
    filter(C_Label %in% c(0, 1, C)) %>% 
    mutate(C_Label = ifelse(C_Label == C, 2, C_Label)) # turn (M+full label) to M+2
}

x.pal <-  f.FAs_2C(whichFA = "Palm", C = 16) 
x.ole <-  f.FAs_2C(whichFA = "Ole",  C = 18) 
x.lino <- f.FAs_2C(whichFA = "Lino", C = 18) 

# x.TAGpalm <- f.FAs_2C(whichFA = "TAGpalm", C = 16) 
# x.TAGole  <- f.FAs_2C(whichFA = "TAGole",  C = 18) 
x.TAGlino <- f.FAs_2C(whichFA = "TAGlino", C = 18) 


d.labeling3 <- d.labeling2 %>% 
  filter(! Compound %in% c("Palm", "Ole", "Lino", "TAGlino")) %>% 
  bind_rows(x.pal, x.ole, x.lino, x.TAGlino) %>% 
  
  # adjust the FAs infusion rate to C2 equivalent
  mutate(infuse.nmol.min.g = ifelse(Infusate == "13CLino",  infuse.nmol.min.g   * 18/2, infuse.nmol.min.g)) %>% 
  # arrange in order
  arrange(State, Infusate, mouse.when.who, tissue, Compound) %>% 
  
  # add C max
  group_by(Compound) %>% 
  mutate(C.max = max(C_Label)) %>% 
  ungroup()

d.labeling3

f.plt_FAs(myData = d.labeling3, myState = "fasted")

d.labeling3$Compound %>% unique()
d.labeling3$C.max %>% unique()




# continue to make the data format consistent for MFA
d.labeling4 <- d.labeling3 %>% 
  rowwise() %>% mutate(Compound.tissue.seq = paste0(1:C.max, collapse = "")) %>% ungroup() %>% 
  mutate(Compound.tissue.seq = paste0(Compound.tissue, "_", Compound.tissue.seq)) %>% 
  mutate(Compound.tissue.seq_m.id = paste0(Compound.tissue.seq, "|", m.id)) %>% 
  
  # calculate labeling sd
  group_by(State, Infusate, tissue, Compound, C_Label) %>% 
  mutate(labeling.sd = sd(labeling) %>% round(7), .after = labeling) %>% 
  ungroup() %>% 
  
  # reorder
  relocate(
    Compound, C_Label, tissue, labeling, labeling.sd, State, Infusate, mouse.when.who, 
    infuse.nmol.min.g, Compound.tissue, m.id, C.max, Compound.tissue.seq, Compound.tissue.seq_m.id)



# the 13C linoleate benchmark infusion rate is 4 nmol/min/g BW used in MFA
# the infusion rate happens to be identical in the kinetics assay; no normalization is needed
# just need to update the infusion rate to AcCoA C2 equivalent flux
d.labeling4 %>% filter(Infusate == "13CLino" & Compound == "Lino")


# Adjust the FAs infusion rate to C2 equivalent
d.labeling5 <- d.labeling4 %>% 
  mutate(infuse.nmol.min.g = ifelse(Infusate == "13CLino", infuse.nmol.min.g * 18/2, infuse.nmol.min.g))


# combine the pulse kinetics labeling with the original MFA labeling dataset
d.13C.fasted.withPulse <- d.labeling5 %>% bind_rows(d.13C.fasted)

save(d.13C.fasted.withPulse, file = "allData_MFA-pulse.RData")



d.13C.fasted.withPulse$Compound %>% unique()
d.13C.fasted.withPulse$Compound.tissue.seq %>% unique()


