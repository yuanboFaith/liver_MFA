rm(list = ls())
library(tidyverse)
library(readxl)
(rstudioapi::getActiveDocumentContext())$path %>% dirname() %>% setwd(); getwd()

theme_set(theme_bw(base_size = 16) + 
            theme(strip.background = element_blank(),
                  panel.spacing.y = unit(20, "pt"),
                  strip.text = element_text(size = 16, face = "bold"),
                  axis.title = element_text(face = "bold")))

# Liver + serum data
d.liver.serum <- read.csv("inulin tracing liver and serum labeling.csv") %>% as_tibble()

d.liver.serum.backgroundSubtracted <- d.liver.serum %>% 
  # select useful columns
  select(isotopeLabel, compound, formula, Blank_1:last_col()) %>% 
  # remove background
  rowwise() %>% 
  mutate(across(Blank_1:last_col(), ~. - max(c(Blank_1, Blank_2)))) %>% 
  mutate(across(Blank_1:last_col(), ~ ifelse(. >0, ., 0)  )) %>% 
  # remove blank
  select(-contains("Blank")) 

l <- d.liver.serum.backgroundSubtracted %>% 
  accucor::natural_abundance_correction(resolution = 120000)

d.normalized <- l$Normalized  

# tidy up
d.labeling.tidy.liver.serum. <- d.normalized %>% 
  pivot_longer(-c(Compound, C_Label), names_to = "tissue", values_to = "labeling") %>% 
  separate(tissue, into = c("tissue", "m.id")) 






# acetate labeling
d.acetate <- read.csv("SCFA cecum feces.csv") %>% as_tibble()

d.acetate.backgroundSubtracted <- d.acetate %>% 
  mutate(compound = str_remove(compound, pattern = "derivatized-")) %>% 
  select(isotopeLabel, compound, formula, blank_1:last_col()) %>% 
  # remove background
  rowwise() %>% 
  mutate(across(blank_1:last_col(), ~. - max(c(blank_1, blank_2, blank_3, blank_1_r, blank_2_r, blank_3_r)))) %>% 
  mutate(across(blank_1:last_col(), ~ ifelse(. >0, ., 0)  )) %>% 
  # remove blank
  select(-contains("blank")) 

# remove higher than C.max labeling (artifact or background)
d.acetate.beforeDeriv <- d.acetate.backgroundSubtracted %>% 
  mutate(No. = str_extract(isotopeLabel, pattern = "\\d$"), .before = 2) %>%
  mutate(No. = ifelse(is.na(No.), "0", No.) %>% as.integer()) %>%
  filter( (compound == "Acetate" & No. <=2 ) | (compound == "Propionate" & No. <=3 ) | (compound == "Butyrate" & No. <=4 )) %>%
  
  select(-c(No.))


# natural abundance correction

l2 <- d.acetate.beforeDeriv %>% 
  accucor::natural_abundance_correction(resolution = 120000)

d.normalized2 <- l2$Normalized  %>% 
  # this step is critical! Remove the C number introduced from derivatization
  filter( (Compound == "Acetate" & C_Label <=2 ) | (Compound == "Propionate" & C_Label <=3 ) | (Compound == "Butyrate" & C_Label <=4 ))


# tidy up
d.labeling.tidy.acetate <- d.normalized2 %>% 
  pivot_longer(-c(Compound, C_Label), names_to = "tissue", values_to = "labeling") %>% 
  separate(tissue, into = c("tissue", "m.id"))

# check labeling pattern to confirm higher than C max labeling is negligible
d.labeling.tidy.acetate %>% ggplot(aes(x = C_Label, y = labeling)) + geom_point() + facet_wrap(~Compound)


# combine liver, serum and acetate data
d.tidy.combined <- rbind(d.labeling.tidy.liver.serum., d.labeling.tidy.acetate) %>% 
  # convert to labeling %
  mutate(label.pct = labeling * 100)

# update mouse id to alphabetical letters
named.m.id <- letters[1:7]
names(named.m.id) <- 1:7
d.tidy.renameMouseID <- d.tidy.combined %>% mutate(m.id = named.m.id[m.id])



# combine with mouse id
d.m.id <- read_excel("mouse id.xlsx") %>%
  mutate(m.id = as.character(m.id)) %>% 
  select(-`time (h) exact`)

d.tidy <- d.tidy.renameMouseID %>% left_join(d.m.id, by = "m.id")





# remove cecum 4 and feces 5 & 6 as their parent SCFAx3 are 100 times lower than the others
# (though the labeling pattern is still rather similar)
# probably the derivatization step is not successful
d.tidy <- d.tidy %>% 
  # filter(! (tissue == "caecum" & m.id %in%  4)) %>% 
  filter(! (tissue == "feces" & m.id %in% c(5)))



# calcualte average carbon labeling
d.average.C.label <- d.tidy %>% 
  group_by(Compound) %>% 
  mutate(C_max = max(C_Label)) %>% 
  group_by(tissue, m.id, Compound, `time (h)`) %>% 
  mutate(enrich = C_Label / C_max * label.pct) %>% 
  summarise(label.pct = sum(enrich)) %>% 
  mutate(tissue_compound = str_c(tissue, "-", Compound))

d.average.C.label







# # plot
# d.average.C.label %>% 
#   ggplot(aes(x = Compound, y =  average.C.labeling)) +
#   geom_point() +
#   facet_wrap(~tissue) +
#   theme(axis.text.x = element_text(angle = 60, hjust = 1))
# 
# 
# 
# d.average.C.label %>% 
#   filter(tissue != "serum") %>% 
#   filter(Compound %in% c("Malate", "Succinic acid", "Acetate", "Propionate", "Butyrate")) %>% 
#   ggplot(aes(x = tissue_compound, y = average.C.labeling, color = m.id)) +
#   geom_point() +
#   geom_line(aes(group = m.id)) +
#   theme_classic(base_size = 14) +
#   theme(axis.text.x = element_text(angle = 40, hjust = 1)) +
#   labs(x = NULL)


SCFA <- c("Acetate", "Propionate", "Butyrate")



# caecum and feces have similar labeling
d.average.C.label %>% 
  filter(tissue %in% c("caecum", "feces")) %>% 
  filter(Compound %in% SCFA ) %>% 
  ggplot(aes(x = tissue, y = label.pct, color = Compound)) +
  facet_grid(`time (h)` ~ Compound, scales = "free") +
  geom_point() +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, .1))) +
  labs (y = "labeling %", x = NULL) +
  theme(legend.position = "none")




# compare 3 h vs 8 h
d.average.C.label %>% 
  filter(tissue %in% c("caecum")) %>% 
  filter(Compound %in% SCFA ) %>% 
  ggplot(aes(x = `time (h)`, y = label.pct, color = Compound)) +
  facet_grid(. ~ Compound, scales = "free") +
  geom_point() +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, .1))) +
  labs (y = "labeling %", x = NULL) +
  theme(legend.position = "none")




# compare 3 h vs 8 h
# SCFA in caecum
d.average.C.label %>% 
  filter(tissue %in% c("caecum")) %>% 
  mutate(Compound = factor(Compound, levels = c(SCFA, "Malate", "Succinic acid"))) %>% 
  ggplot(aes(x = `time (h)`, y = label.pct, color = Compound)) +
  facet_wrap(. ~ Compound, nrow = 1) +
  geom_point() +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, .1))) +
  labs (y = "labeling %", x = NULL) +
  theme(legend.position = "none")

# TCA in liver
d.average.C.label %>% 
  filter(tissue %in% c("liver")) %>% 
  filter(Compound %in% c("Malate", "Succinic acid")) %>% 
  ggplot(aes(x = `time (h)`, y = label.pct, color = Compound)) +
  facet_wrap(. ~ Compound, nrow = 1) +
  # geom_point() +
  geom_text(aes(label = m.id)) +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, .1))) +
  labs (y = "labeling %", x = NULL) +
  theme(legend.position = "none")


# normalized liver labeling
d.SCFA.label <- d.average.C.label %>% 
  filter(tissue == "caecum") %>% 
  group_by(`time (h)`, m.id) %>% 
  summarise(SCFA.label.pct = mean(label.pct))


d.TCA.label <- d.average.C.label %>% 
  filter(tissue == "liver") %>% 
  filter(Compound %in% c("Malate", "Succinic acid")) %>% 
  group_by(`time (h)`, m.id, Compound) %>% 
  summarise(TCA.label.pct = mean(label.pct))

d.SCFA.energy <- left_join(d.SCFA.label, d.TCA.label) %>% 
  mutate(TCA.label.pct.normalized = TCA.label.pct * 100 / SCFA.label.pct)

d.SCFA.energy %>% 
  ggplot(aes(x = Compound, y = TCA.label.pct.normalized, colour = `time (h)`)) +
  # geom_text(aes(label = m.id)) +
  geom_point() + 
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, .1))) 






# organize data into structure suitable for MFA
d.tidy$Compound %>% unique()

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
  "Acetate",               "Acetate",
  "Propionate",            "Propionate",
  "Butyrate",              "Butyrate")   

CC <- d.metab_names$Compound
aa <- d.metab_names$abbr
names(aa) <- CC; aa


d.dataForMFA.SCFA <- d.tidy %>% 
  
  # rename tissue
  filter(tissue %in% c("liver", "caecum" , "serum")) %>% 
  mutate(tissue = str_replace(tissue, "liver", replacement = "Lv"),
         tissue = str_replace(tissue, "serum", replacement = "Blood"),
         tissue = str_replace(tissue, "caecum", replacement = "hp")) %>% 
  
  filter(! Compound %in% c("C16:0", "C18:1", "C18:2")) %>% 
  
  # simplify compound names
  mutate(Compound = str_remove(Compound, "derivatized-")) %>% 
  mutate(Compound = aa[Compound]) %>% 
  # not all compounds have their derivatization listed in the named vector, e.g., valine, and their names will be turned to NA
  # they are not useful to provide more information for MFA, so here for simplicity, they are removed
  filter(! is.na(Compound)) %>% 
  
  # standard deviation
  group_by(tissue, Compound, C_Label) %>%  
  mutate(labeling.sd = sd(labeling)) %>% ungroup() %>% 
  
  # calculate C -max; manually adjust SCFA C max (not considering derivatization)
  group_by(Compound) %>% mutate(C.max = max(C_Label)) %>%  # C max
  mutate(C.max = ifelse(Compound == "Acetate",    2, C.max)) %>% 
  mutate(C.max = ifelse(Compound == "Propionate", 3, C.max)) %>% 
  mutate(C.max = ifelse(Compound == "Butyrate",   4, C.max)) %>% 
  
  
  mutate(State = "refed",
         Infusate = paste0("13ChpSCFA", m.id)) %>% 
  mutate(mouse.when.who = paste0(m.id, "_2025Dec_BY")) %>% 
  mutate(infuse.nmol.min.g = NA) %>% 
  mutate(Compound.tissue          = paste0(Compound, ".", tissue)) %>% 
  
  # rowwise important to add the correct C-max
  rowwise() %>% 
  mutate(Compound.tissue.seq      = str_c(Compound.tissue, "_", str_c(1:C.max, collapse = ""))) %>% ungroup() %>% 
  mutate(Compound.tissue.seq_m.id = str_c(Compound.tissue.seq, "|", m.id)) %>% 
  
  relocate(Compound, C_Label, tissue, labeling, labeling.sd, State, Infusate, mouse.when.who, 
           infuse.nmol.min.g, Compound.tissue, m.id, C.max, Compound.tissue.seq, Compound.tissue.seq_m.id) %>% 
  
  select(-c(`Start feeding`, "SAC", "time (h)", label.pct)) %>% 
  
  # use only time point of 3h
  filter(m.id %in% c("a", "b", "c", "d")) # a-d for 3 h, and e, f, and g for 8 h


d.dataForMFA.SCFA



# plot isotopologue distribution (before simplified into averaged carbon labeling)
p.SCFA.isotopologues <- d.dataForMFA.SCFA %>% 
  filter(Compound %in% c("Acetate", "Butyrate", "Propionate")) %>% 
  filter(tissue == "hp") %>% 
  filter(C_Label <=4) %>% 
  ggplot(aes(x = m.id, y = labeling, fill = as.character(C_Label))) +
  geom_col(color = "black", position = position_stack(reverse = T)) +
  geom_text(aes(label = paste0("M+", C_Label)), position = position_stack(reverse = T, vjust = .63), size = 3.5) +
  geom_text(aes(label = labeling %>% round(2)), position = position_stack(reverse = T, vjust = .37), size = 3.5) +
  facet_wrap(~Compound) +
  scale_y_continuous(expand = expansion(mult = c(0, 0))) +
  theme(panel.border = element_blank()) +
  coord_cartesian(ylim = c(.8, 1)) 
p.SCFA.isotopologues






# Clean up the average carbon labeling to plug into the `d.dataForMFA.SCFA`
d.plugIn <- d.average.C.label %>% 
  ungroup() %>% 
  filter(Compound %in% SCFA) %>% 
  filter(tissue == "caecum") %>% 
  filter(`time (h)` == "3 h") %>% 
  select(m.id, Compound, label.pct) %>% 
  
  # convert to fraction and rename in accordance to 'd.dataForMFA.SCFA'
  mutate(label.pct = label.pct / 100) %>% 
  rename(labeling = label.pct) %>% 
  
  # calculate labeling sd
  group_by(Compound) %>% 
  mutate(labeling.sd = sd(labeling)) %>% 
  arrange(Compound) %>% 
  
  # add C. max for which we have the average labeling
  mutate(C.max = ifelse(Compound == "Acetate", 2, ifelse (Compound == "Butyrate", 4, 3)))

d.plugIn



# Replace the isotopologue c max with average labeling  
x <- d.dataForMFA.SCFA %>% filter( Compound %in% SCFA) # SCFA
y <- d.dataForMFA.SCFA %>% filter(!Compound %in% SCFA) # others

# SCFA c max using average labeling
x.SCFA.Cmax <- x %>% filter(C_Label == C.max) %>% select(-c(labeling, labeling.sd)) %>% left_join(d.plugIn) %>% 
  relocate(Compound, C_Label, tissue, labeling, labeling.sd, State, Infusate, mouse.when.who, 
           infuse.nmol.min.g, Compound.tissue, m.id, C.max, Compound.tissue.seq, Compound.tissue.seq_m.id)
  
# SCFA rest being M+0 
x.SCFA.C0 <- x.SCFA.Cmax %>% mutate(labeling = 1 - labeling, C_Label = 0)
  
x <- bind_rows(x.SCFA.Cmax, x.SCFA.C0)  

# reconstruct back
d.dataForMFA.SCFA <- bind_rows(x, y)



# plot isotopologue distribution
p.SCFA.averaged <- d.dataForMFA.SCFA %>% 
  filter(Compound %in% c("Acetate", "Butyrate", "Propionate")) %>% 
  filter(tissue == "hp") %>% 
  filter(C_Label <=4) %>% 
  ggplot(aes(x = m.id, y = labeling, fill = as.character(C_Label))) +
  geom_col(color = "black", position = position_stack(reverse = T)) +
  geom_text(aes(label = paste0("M+", C_Label)), position = position_stack(reverse = T, vjust = .63), size = 3.5) +
  geom_text(aes(label = labeling %>% round(2)), position = position_stack(reverse = T, vjust = .37), size = 3.5) +
  facet_wrap(~Compound) +
  scale_y_continuous(expand = expansion(mult = c(0, 0))) +
  theme(panel.border = element_blank()) +
  coord_cartesian(ylim = c(.8, 1)) 
p.SCFA.averaged




# Simplify Butyrate's M+4 as M+2
# this is in particularly convenient to calculate acetate and butyrate summation contribution confidence interval
d.dataForMFA.SCFA <- d.dataForMFA.SCFA %>% 
  filter(! (Compound == "Butyrate" & C_Label %in% c(2, 3))) %>%  # remove M+2 and M+3
  mutate(C_Label = ifelse(Compound == "Butyrate" & C_Label == 4, 2, C_Label))


p.SCFA.averaged.simplified <- d.dataForMFA.SCFA %>% 
  filter(Compound %in% c("Acetate", "Butyrate", "Propionate")) %>% 
  filter(tissue == "hp") %>% 
  ggplot(aes(x = m.id, y = labeling, fill = as.character(C_Label))) +
  geom_col(color = "black", position = position_stack(reverse = T)) +
  geom_text(aes(label = paste0("M+", C_Label)), position = position_stack(reverse = T, vjust = .63), size = 3.5) +
  geom_text(aes(label = labeling %>% round(2)), position = position_stack(reverse = T, vjust = .37), size = 3.5) +
  facet_wrap(~Compound) +
  scale_y_continuous(expand = expansion(mult = c(0, 0))) +
  theme(panel.border = element_blank()) +
  coord_cartesian(ylim = c(.8, 1)) 
p.SCFA.averaged.simplified


cowplot::plot_grid(p.SCFA.isotopologues, p.SCFA.averaged, p.SCFA.averaged.simplified, nrow = 3)


save(d.dataForMFA.SCFA, file = "../data/data_section_SCFA.RData") 


