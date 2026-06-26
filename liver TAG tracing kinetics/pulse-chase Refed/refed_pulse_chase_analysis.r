rm(list=ls())

library(rebus)
library(tidyverse)
library(readxl)
library(accucor)

theme_set(
  theme_classic(base_size = 16) +
    theme(
      # panel.border = element_rect(colour = "black", fill = NA, linewidth = .3),
      strip.background = element_blank(),
      axis.title.y = element_text(margin = margin(r = 10, unit = "pt")),
      axis.title = element_text(face = "bold", colour = "black"),
      axis.text = element_text(colour = "black"),
      strip.text = element_text(size = 12, face = "bold")) 
)


rstudioapi::getActiveDocumentContext()$path %>% dirname() %>% setwd(); getwd()


myPath <- "refed labeling data.xlsx"

# serum 
func.process <- function(whichSheet = "serum"){
  
  x <- (read_excel(myPath, sheet = whichSheet) %>% 
          natural_abundance_correction(resolution = 120000)) $ Normalized
  
  x %>% pivot_longer(-c(1, 2), names_to = "mouse", values_to = "labeling") %>% 
    arrange(Compound, mouse, C_Label) %>% 
    return()
}

x1 <- func.process(whichSheet = "serum");   x1 
x2 <- func.process(whichSheet = "liver");   x2
x3 <- func.process(whichSheet = "saponification"); x3 



# check saponification samples 5 uL vs 1 uL injection volume
# ensure that the detector is not saturated, and that the M+18 labeling is the same regardless of the injection volume
x3 %>% filter(C_Label == 18) %>% 
  filter(!str_detect(mouse, "K")) %>%  # remove the artificial t=0 marker
  mutate(inj = rep(c(1, 5), each = 10),
         mouse = str_extract(mouse, "[A-Z]$")) %>% 
  ggplot(aes(x = mouse, y = labeling, color = factor(inj))) + 
  geom_point(size = 4, shape = 21, stroke = 1) +
  theme(legend.position = c(.9, .2))

# 1 uL & 5 uL inj. volume result is very consistent. Select only the 5 uL injection result
x3.5uL <- x3 %>% filter(!str_detect(mouse, "1uL"))



# clean up the data structure
d.labeling.0 <- bind_rows(x1, x2, x3.5uL)
d.labeling.0$mouse %>% unique()

d.labeling.1 <- d.labeling.0 %>% 
  mutate(sample = mouse) %>%  
  mutate(mouse = str_extract(sample, "(?<=_)[A-Z]$")) %>% 
  mutate(tissue = str_extract(sample, "[a-zA-Z]{1,10}(?=_[A-Z])")) %>% 
  mutate(treatment = ifelse(str_detect(sample, "sap"), "sap", "metab")) %>% 
  # correct the tissue for the saponified group
  mutate(tissue = ifelse(treatment == "sap", "Liver", tissue))

d.labeling.1




# calculate average labeling
d.labeling.2 <- d.labeling.1 %>% 
  group_by(mouse, tissue, Compound, sample) %>% 
  mutate(C.max = max(C_Label),
         labeling.weighted = C_Label / C.max * labeling) %>% 
  summarise(labeling = sum(labeling.weighted))




# clean up sample ID and combine with labeling data
# PULSE time counted from START of infusion
d.id.pulse <- read_excel(myPath, sheet = "sampleID_pulse") %>% select(mouse, `time (h)`)
# CHASE time counted from END of infusion
d.id.chase <- read_excel(myPath, sheet = "sampleID_chase") %>% select(mouse, `time (h)`) 
d.id <- bind_rows(d.id.pulse, d.id.chase)

# visualize: change time counting from the start of the infusion
d.id.chase2 <- d.id.chase %>% mutate(`time (h)` = `time (h)` + 8) # time counted from the start of the infusion
d.id.all <- bind_rows(
  d.id.pulse  %>% mutate(method = "pulse"), 
  d.id.chase2 %>% mutate(method = "chase")) %>% 
  arrange(`time (h)`)


d.labeling.3 <- d.labeling.2 %>% left_join(d.id.all, by = "mouse")
d.labeling.3$mouse %>% unique()

# visualize

##  'sources', 'targets', 'myColors' are borrowed from 'liver_kinetics_pulse.RData'
load("../Pulse Fasted/liver_kinetics_pulse.RData")

func.plt <- function(md = d.labeling.3){
  md %>%
    mutate(`time (h)` = floor(`time (h)`)) %>% 
    ggplot(aes(x = `time (h)`, y = labeling, color = Compound)) +
    annotate(geom = "rect", xmin = 0, xmax = 8, ymin = -Inf, ymax = Inf, fill = "purple", alpha = .04) +
    annotate(geom = "rect", xmin = 8, xmax = 9.5, ymin = -Inf, ymax = Inf, fill = "green3", alpha = .05) +
    geom_vline(xintercept = 8, linewidth = .5, linetype = "dashed", alpha = .4) +
    # ggbeeswarm::geom_quasirandom(size = 2, shape = 19, varwidth = .1, alpha = .7) +
    stat_summary(fun = mean, aes(group = Compound), geom = "line") +
    geom_text(aes(label = mouse), size = 5) +
    # facet_wrap(~Compound, scales = "free") +
    scale_x_continuous(breaks = seq(0, 9, 1), expand = expansion(mult = c(0, .03)))+
    scale_y_continuous(expand = expansion(mult = c(0, .03)), n.breaks = 7,
                       labels = ~.x * 100,
                       name = "labeling (%)")  +
    coord_cartesian(xlim = c(0, NA), ylim = c(0, NA)) +
    theme(legend.position = "right") 
  # stat_smooth(
  #   method = "nls",
  #   formula = y ~ a * exp(-b * x) + c,
  #   method.args = list(start = list(a = .05, b = .05, c = 30)),
  #   se = FALSE,
  #   fullrange = T,
  #   # color = "blue"
  # )
}


# source labeling
p1 <- d.labeling.3 %>% 
  filter( (tissue == "Serum" & Compound %in% sources) | 
            (tissue == "Liver" & Compound == "TAG C18:2") ) %>%
  func.plt() +
  scale_color_manual(values = myColors) 
p1    


# TCA labeling
p2 <- d.labeling.3 %>% 
  filter(tissue == "Liver" & Compound %in% targets) %>% func.plt() 
p2

# combine
cowplot::plot_grid(p1 + theme(legend.position = "none"), 
                   ggplot() + theme_void(),
                   p2 + theme(legend.position = "none"), 
                   nrow = 1, 
                   rel_widths = c(1, .1, 1))

ggsave(filename = "refed pulse - chase.pdf", height = 3.4, width = 7.3)







# Organize data to fit for MFA
d.labeling.1$Compound %>% unique()


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


d.labeling.1 %>% select(Compound, tissue, treatment) %>% table()
d.labeling.1$tissue %>% unique()

d.dataForMFA.TAGkinetics <- d.labeling.1 %>% 
  
  # select relevant compounds; e.g., lactate appear in both liver and serum, but use the serum data
  filter( (tissue == "Liver" & Compound %in% c("Malate", "Succinic acid", "TAG C18:2")) | 
           tissue == "Serum" & Compound %in% c(sources) ) %>% 
  # remove redundant columns
  select(-c(sample, treatment)) %>% 
  
  # rename compounds using the unified abbreviations
  mutate(Compound = aa[Compound]) %>% 
  
  # rename mouse id
  rename(m.id = mouse) %>% 
  
  # rename tissue
  mutate(tissue = str_replace(tissue, "Liver", replacement = "Lv"),
         tissue = str_replace(tissue, "Serum", replacement = "Blood")) %>% 
  
  # standard deviation
  group_by(tissue, Compound, C_Label) %>%  
  mutate(labeling.sd = sd(labeling)) %>% ungroup() %>% 
  mutate(labeling.sd =  replace_na(labeling.sd, 0)) %>% 
  
  
  mutate(State = "refed",
         Infusate = paste0("13CLinoKin", m.id)) %>% 
  mutate(mouse.when.who = paste0(m.id, "_2025Dec_BY")) %>% 
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
  
  filter(m.id != "K") # K is the time = 0 origin marker for plotting purpose
  






# a quick visual check
d.dataForMFA.SCFA %>% 
  filter(Compound %in% c("Lino", "TAGLino" )) %>% 
  ggplot(aes(x = Compound, y = labeling, fill = as.character(C_Label))) +
  geom_col() + facet_wrap(~m.id)


d.dataForMFA.TAGkinetics$m.id %>% unique()



save(d.dataForMFA.TAGkinetics, file = "../../data/data_section_TAGkinetics.RData") 



d.dataForMFA.TAGkinetics %>% filter(Compound == "Lino") %>% tail()
d.dataForMFA.TAGkinetics$Compound %>% unique()





# # # calculate circulating C18:2 and liver TAG-C18:2 contribution to TCA cycle # --------------
# d.labeling.4 <- d.labeling.3 %>% 
#   filter((tissue == "Serum" & Compound %in% sources) |
#            (tissue == "Liver" & Compound %in% c("TAG C18:2", "Malate"))) %>% 
#   filter(mouse != "K") %>%  # remove t=0 time marker
#   ungroup() %>% 
#   select(-c(sample, tissue)) %>% 
#   pivot_wider(names_from = Compound, values_from = labeling)
# 
# 
# # Load the carbon-wise EMU labeling (find pyruvate, glutamine, and CO2 carbon - specific labeling)
# env <- new.env()
# load("../../LvM_Refed/5_2_optimal_solution.RData", envir = env)
# ls(env)  # List all objects in the file
# l.EMU.sim <- env$l.EMU.sim  # Extract just the list of simulated carbon-wise labeling
# 
# 
# 
# # liver pyruvate labeling (equal to blood lactate labeling)
# Pyr.Lv_1 <- l.EMU.sim $ `13CLino` $ Pyr.Lvc_1[2] # PHD, pyruvate 2-3 enters acetyl-CoA
# Pyr.Lv_2 <- l.EMU.sim $ `13CLino` $Pyr.Lvc_2[2] 
# Pyr.Lv_3 <- l.EMU.sim $ `13CLino` $Pyr.Lvc_3[2] 
# Pyr.Lv.mean <- mean(Pyr.Lv_1, Pyr.Lv_2, Pyr.Lv_3)
# Pyr.Lv_23   <- mean(Pyr.Lv_2, Pyr.Lv_3)
# 
# # blood glutamine labeling
# Gln.Blood_1 <- l.EMU.sim $ `13CLino` $Gln.Blood_1[2] # Gln -> aKG -> Suc: 1st carbon lost to CO2
# Gln.Blood_2 <- l.EMU.sim $ `13CLino` $Gln.Blood_2[2]
# Gln.Blood_3 <- l.EMU.sim $ `13CLino` $Gln.Blood_3[2]
# Gln.Blood_4 <- l.EMU.sim $ `13CLino` $Gln.Blood_4[2]
# Gln.Blood_5 <- l.EMU.sim $ `13CLino` $Gln.Blood_5[2]
# Gln.Blood_mean <- mean(Gln.Blood_1, Gln.Blood_2, Gln.Blood_3, Gln.Blood_4, Gln.Blood_5)
# Gln.Blood_2345 <- mean(Gln.Blood_2, Gln.Blood_3, Gln.Blood_4, Gln.Blood_5)
# 
# # CO2 labeling (assume proportional to blood linoleate labeling)
# CO2.Blood <- l.EMU.sim $ `13CLino` $CO2.Blood_1[2]
# 
# # blood linoleate labeling
# Lino.Blood_1 <- l.EMU.sim $ `13CLino` $Lino.Blood_1[2]
# Lino.Blood_2 <- l.EMU.sim $ `13CLino` $Lino.Blood_2[2]
# Lino.Blood_mean <- mean(Lino.Blood_1, Lino.Blood_2)
# 
# 
# 
# 
# # annotate fluxes with reaction marker (enzymes) names 
# # fluxes are in the unit of molecules of products produced, nmol /min/ g body weight
# 
# # load enzyme index from the fasted folder (shared by both fasting and refed reduced modeling)
# d.enz <- read_excel(path = "../../LvM_Fasted/LvM.xlsx") %>%  
#   filter(!is.na(Marker)) %>% 
#   select(enzyme, reactions, Marker)
# 
# optimal_keyFlux <- tibble(flux = c(env$optimal_allFlux)) %>% 
#   mutate(enzyme = 1:nrow(.)) %>% # flux created 
#   right_join(d.enz, by = "enzyme")   # augment with name of enzymes
# 
# PEPCK <- optimal_keyFlux[optimal_keyFlux$Marker == "PEPCK", ]$flux ; PEPCK
# CS    <- optimal_keyFlux[optimal_keyFlux$Marker == "CS", ]   $flux ; CS
# ME    <- optimal_keyFlux[optimal_keyFlux$Marker == "ME", ]   $flux ; ME
# PEPCK.CS.ME <- PEPCK + CS + ME
# 
# PC    <- optimal_keyFlux[optimal_keyFlux$Marker == "PC", ]   $flux ; PC
# PDH   <- optimal_keyFlux[optimal_keyFlux$Marker == "PDH", ]  $flux ; PDH
# GLS   <- optimal_keyFlux[optimal_keyFlux$Marker == "GLS", ]  $flux ; GLS
# 
# 
# 
# # calculate the time point-specific partial-molecular average labeling
# d.pulse.chase.atomFlux <- d.labeling.4 %>% 
#   mutate(Glutamine_2345 = Gln.Blood_2345 / Gln.Blood_mean * Glutamine,
#          CO2 = CO2.Blood / Lino.Blood_mean * `C18:2`,
#          Pyr = Lactate,
#          Pyr_23 = Pyr.Lv_23 / Pyr.Lv.mean * Pyr) %>% 
#   mutate(b = 2 * Malate * (PEPCK.CS.ME) - PDH * Pyr_23 - 1.5 * PC * Pyr - .5 * PC * CO2 - 2 * GLS * Glutamine_2345 )
# 
# 
# d.pulse.chase.atomFlux <- d.pulse.chase.atomFlux %>% filter(mouse != "A")
# 
# mdl <- lm(data = d.pulse.chase.atomFlux, formula = b ~ `TAG C18:2` + `C18:2` + 0 )
# mdl
# 
# summary(mdl) # print the standard error of the coefficient
# 
# # Check fitting
# d.pulse.chase.atomFlux.fitted <- d.pulse.chase.atomFlux %>% 
#   mutate(fitted = mdl$fitted.values)
# 
# d.pulse.chase.atomFlux.fitted %>% 
#   ggplot(aes(x = b, y = fitted, color = `time (h)`)) + 
#   #  geom_point(size = 3) +
#   geom_text(aes(label = mouse)) +
#   geom_abline(slope = 1, intercept = 0) +
#   coord_cartesian(xlim = c(0, NA), ylim = c(0, NA)) +
#   scale_color_distiller(palette = "Spectral")
# 
# 
# 
# # extract NEFA oxidation fluxes from reduced model, without TAG in the model yet
# # this data is based on Tony's fatty acid flux data
# # use the ratio of Palm : Ole: Lino from Tony's data to extrapolate 
# # Palm and Ole NEFA and TAG 's contribution, based on Bo's kinetics tracing from 36 mM linoleate
# f.Palm.SH   <- optimal_keyFlux[optimal_keyFlux$Marker == "NEFA-Palm", ] $ flux ; f.Palm.SH
# f.Ole.SH    <- optimal_keyFlux[optimal_keyFlux$Marker == "NEFA-Ole",  ] $ flux ; f.Ole.SH
# f.Lino.SH   <- optimal_keyFlux[optimal_keyFlux$Marker == "NEFA-Lino", ] $ flux ; f.Lino.SH
# f.others.SH <- optimal_keyFlux[optimal_keyFlux$Marker == "others",    ] $ flux ; f.others.SH
# 
# 
# # use Bo's data to calculate the actual fatty acid & TAG beta oxidation flux
# f.TAG_Lino  <- mdl$coefficients[1]; f.TAG_Lino  
# f.NEFA_Lino <- mdl$coefficients[2]; f.NEFA_Lino 
# 
# f.TAG_Ole  <- f.Ole.SH / f.Lino.SH  * f.TAG_Lino ; f.TAG_Ole
# f.TAG_Palm <- f.Palm.SH / f.Lino.SH * f.TAG_Lino ; f.TAG_Palm
# 
# f.NEFA_Ole  <- f.Ole.SH / f.Lino.SH  * f.NEFA_Lino ; f.NEFA_Ole
# f.NEFA_Palm <- f.Palm.SH / f.Lino.SH * f.NEFA_Lino ; f.NEFA_Palm
# 
# f.PDH <- PDH # calculated above
# 
# f.others <- CS - (f.TAG_Lino + f.TAG_Ole + f.TAG_Palm) - (f.NEFA_Lino + f.NEFA_Ole + f.NEFA_Palm) - f.PDH
# 
# # error
# e.TAG_Lino  <- (summary(mdl)$coefficients)[1, 2] ; e.TAG_Lino
# e.NEFA_Lino <- (summary(mdl)$coefficients)[2, 2] ; e.NEFA_Lino
# 
# e.TAG_Ole   <- f.Ole.SH  / f.Lino.SH  * e.TAG_Lino ; e.TAG_Ole
# e.TAG_Palm  <- f.Palm.SH / f.Lino.SH  * e.TAG_Lino ; e.TAG_Palm
# 
# e.NEFA_Ole   <- f.Ole.SH  / f.Lino.SH  * e.NEFA_Lino ; e.NEFA_Ole
# e.NEFA_Palm  <- f.Palm.SH / f.Lino.SH  * e.NEFA_Lino ; e.NEFA_Palm
# 
# # error of PDH
# load("../../LvM_Refed/6_CI_refed.RData")
# CI.PDH <- d.CI.bounds.reactions.lv[d.CI.bounds.reactions.lv$reactions == "Pyr.Lvm->AcCoA.Lv+CO2.Blood", ] # C.I. of PDH flux in liver
# e.PDH <- (CI.PDH$R - CI.PDH$L) /2 / 2 # (top - bottom) / 2 to be CI of one side, take half as the SEM 
# 
# # error of others
# e.others <- sqrt(e.TAG_Lino^2 + e.TAG_Ole^2 + e.TAG_Palm^2 + e.NEFA_Lino^2 + e.NEFA_Ole^2 + e.NEFA_Palm^2 + e.PDH^2)
# 
# # plot contribution sources to liver acetyl-CoA
# d.liver.energy <- tibble(
#   source = c("NEFA-Palmitate", "NEFA-Oleate", "NEFA-Linoleate", "TAG-Palmitate", "TAG-Oleate", "TAG-Linoleate",     "others",      "Carbs"),
#   flux   = c(f.NEFA_Palm,         f.NEFA_Ole,     f.NEFA_Lino,       f.TAG_Palm,    f.TAG_Ole,      f.TAG_Lino,     f.others,        f.PDH),
#   error  = c(e.NEFA_Palm,         e.NEFA_Ole,     e.NEFA_Lino,       e.TAG_Palm,    e.TAG_Ole,      e.TAG_Lino,     e.others,        e.PDH)
# ) %>% 
#   mutate(source = factor(source, levels = .$source, ordered = T)) %>% 
#   arrange(desc(source)) %>% 
#   mutate(errY = cumsum(flux))
# 
# p <- d.liver.energy %>% 
#   ggplot(aes(x = 1, y = flux, fill = source, alpha = source)) +
#   geom_col(color = "black") +
#   # coord_polar(theta = "y") +
#   # scale_x_continuous(limits = c(0, 1.5)) +
#   # theme_void() +
#   scale_fill_manual(values = c(
#     "Carbs" = "tomato",
#     "NEFA-Palmitate" = "purple", "NEFA-Oleate" = "orange",  "NEFA-Linoleate" = "skyblue2", 
#     "TAG-Palmitate" = "purple",  "TAG-Oleate"  = "orange",  "TAG-Linoleate"  = "skyblue2", 
#     "others" = "grey")
#   ) +
#   geom_errorbar(aes(ymax = errY, ymin = errY - error), width = .2, alpha = 1) +
#   scale_alpha_manual(values = c(1, 1, 1, .2, .2, .2, 1, 1)) +
#   scale_y_continuous(expand = expansion(mult = c(0, .1)), 
#                      breaks = seq(0, 1000, 100)) +
#   scale_x_continuous(expand = expansion(add = .3)) +
#   theme(axis.ticks.x = element_blank(),
#         axis.text.x = element_blank(),
#         axis.text = element_text(colour = "black"),
#         axis.title.x = element_blank()) +
#   theme(legend.position = "none")
# p
# 
# ggsave("liver energy source.pdf", height = 4, width = 2)
# 
# 
# 
# # add label
# p +
#   geom_text(aes(label = paste(source, " - ", round(flux/sum(flux) * 100, 1), "%"), 
#                 y = flux, 
#                 x = .7), 
#             position = position_stack(vjust = .5), 
#             color = "black", alpha = 1, hjust = 0)  

