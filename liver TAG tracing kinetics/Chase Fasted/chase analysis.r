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
d.labeling.2 <- d.labeling.1 %>% left_join(d.id)



# calculate average labeling
d.labeling.3 <- d.labeling.2 %>% 
  group_by(tissue, Compound, sample, `time (h)`) %>% 
  mutate(C.max = max(C_Label),
         labeling.weighted = C_Label / C.max * labeling) %>% 
  summarise(labeling = sum(labeling.weighted))






# plots

##  'sources', 'targets' are borrowed from 'liver_kinetics_pulse.RData'
load("../Pulse Fasted/liver_kinetics_pulse.RData")


func.plt.chase <- function(md = d.average.C.labeling){
  md %>% 
    ggplot(aes(x = `time (h)`, y = labeling, color = Compound)) +
    geom_point(size = 2) +
    # geom_text(aes(label = mouse), size = 5) +
    # facet_wrap(~Compound, scales = "free") +
    scale_x_continuous(breaks = seq(0, 2, .5))  +
    scale_y_continuous(expand = expansion(mult = c(0, .03)), n.breaks = 7,
                       labels = ~.x * 100,
                       name = "labeling (%)")  +
    coord_cartesian(xlim = c(0, NA), ylim = c(0, NA)) +
    theme(legend.position = "right",
          axis.text = element_text(size = 17)) +
    stat_smooth(
      method = "nls",
      formula = y ~ a * exp(-b * x) + c,
      method.args = list(start = list(a = .05, b = .05, c = 30)),
      se = FALSE,
      fullrange = T, 
      # color = "blue"
    )
}


# source labeling
p1 <- d.labeling.3 %>% 
  filter((tissue == "serum" & Compound %in% sources) |
           (tissue == "liver" & Compound == "TAG C18:2")) %>%
  func.plt.chase() +
  
  # add mouse label
  # ggrepel::geom_text_repel(aes(label = sample %>% str_extract("\\d{1,2}")), size = 8, min.segment.length = 0) +
  
  # add regression to TAG C18:2
  stat_smooth(
    data = d.labeling.3 %>% filter(Compound == "TAG C18:2"),
    method = "nls",
    formula = y ~ a * exp(-b * x),
    method.args = list(start = list(a = .1, b = .1)),
    se = FALSE,
    fullrange = T, 
    # color = "blue"
  ) +
  scale_color_manual(values = myColors) 

p1


# TCA labeling
p2 <- d.labeling.3 %>% 
  filter(tissue == "liver" & Compound %in% targets) %>% func.plt.chase() 


# combine
cowplot::plot_grid(p1 + theme(legend.position = "none"), 
                   ggplot() + theme_void(),
                   p2 + theme(legend.position = "none"), 
                   nrow = 1, 
                   rel_widths = c(1, .1, 1))

ggsave(filename = "chase.pdf", height = 4.5, width = 12)



# A <- c(5, 12, 12, 12, 8, 12, 6, 2) %>% matrix(nrow = 4, byrow = T)
# b = c(1.5, 2.5, 2, .6)
# result <- lm.fit(A, b)
# x <- result$coefficients
# x
# 
# A %*% x


# output 
d.chase <- d.labeling.3
d.chase$Compound %>% unique()


# Combine pulse and chase data together
# chase data
d.chase.clean <- d.chase %>% 
  filter(sample %>% str_detect("SAC")) %>% 
  filter(Compound %in% c("C18:2", "TAG C18:2", "Lactate", "Malate", "Glutamine", "3-Hydroxybutyric acid")) %>% ungroup() %>% 
  mutate(mouse = str_extract(sample, "\\d$")) %>% 
  select(-c(sample, tissue)) %>% 
  pivot_wider(names_from = Compound, values_from = labeling) %>% 
  mutate(method = "chase")


# pulse data
d.pulse.clean <- d.pulse %>% 
  filter(Compound %in% c("C18:2", "TAG C18:2", "Lactate", "Malate", "Glutamine", "3-Hydroxybutyric acid")) %>% ungroup() %>% 
  select(-tissue) %>% 
  pivot_wider(names_from = Compound, values_from = labeling) %>% 
  select(-batch) %>% relocate('time (h)') %>% 
  mutate(method = "pulse")
d.pulse.clean <- d.pulse.clean[complete.cases(d.pulse.clean), ]

# stack chase - pulse data
d.pulse.chase <- rbind(d.pulse.clean, d.chase.clean)





# Load the carbon-wise EMU labeling (find pyruvate, glutamine, and CO2 carbon - specific labeling)
env <- new.env()
load("../../LvM_Fasted/5_2_optimal_solution.RData", envir = env)
ls(env)  # List all objects in the file
l.EMU.sim <- env$l.EMU.sim  # Extract just the list of simulated carbon-wise labeling



# liver pyruvate labeling (equal to blood lactate labeling)
Pyr.Lv_1 <- l.EMU.sim $ `13CLino` $ Pyr.Lvc_1[2] # PHD, pyruvate 2-3 enters acetyl-CoA
Pyr.Lv_2 <- l.EMU.sim $ `13CLino` $ Pyr.Lvc_2[2] 
Pyr.Lv_3 <- l.EMU.sim $ `13CLino` $ Pyr.Lvc_3[2] 
Pyr.Lv.mean <- mean(Pyr.Lv_1, Pyr.Lv_2, Pyr.Lv_3)
Pyr.Lv_23   <- mean(Pyr.Lv_2, Pyr.Lv_3)

# liver citrate #6 labeling to lose to CO2
Cit.Lv_6 <- l.EMU.sim $ `13CLino` $ Cit.Lv_6[2] 
# OAA.Lv_1 <- l.EMU.sim $ `13CLino` $ OAA.Lv_1[2] 


# malate labeling
Mal.Lv_1 <- l.EMU.sim $ `13CLino` $ Mal.Lv_1[2]  ; Mal.Lv_1
Mal.Lv_2 <- l.EMU.sim $ `13CLino` $ Mal.Lv_2[2]  ; Mal.Lv_2 
Mal.Lv_3 <- l.EMU.sim $ `13CLino` $ Mal.Lv_3[2]  ; Mal.Lv_3
Mal.Lv_4 <- l.EMU.sim $ `13CLino` $ Mal.Lv_4[2]  ; Mal.Lv_4
Mal.Lv.mean <- mean(Mal.Lv_1, Mal.Lv_2, Mal.Lv_3, Mal.Lv_4)

# aKG labeling
aKG.Lv_1 <- l.EMU.sim $ `13CLino` $ aKG.Lv_1[2]  ; aKG.Lv_1
aKG.Lv_2 <- l.EMU.sim $ `13CLino` $ aKG.Lv_2[2]  ; aKG.Lv_2
aKG.Lv_3 <- l.EMU.sim $ `13CLino` $ aKG.Lv_3[2]  ; aKG.Lv_3
aKG.Lv_4 <- l.EMU.sim $ `13CLino` $ aKG.Lv_4[2]  ; aKG.Lv_4
aKG.Lv_5 <- l.EMU.sim $ `13CLino` $ aKG.Lv_5[2]  ; aKG.Lv_5
aKG.Lv.mean <- mean(aKG.Lv_1, aKG.Lv_2, aKG.Lv_3, aKG.Lv_4, aKG.Lv_5)


# blood glutamine labeling
Gln.Blood_1 <- l.EMU.sim $ `13CLino` $Gln.Blood_1[2] # Gln -> aKG -> Suc: 1st carbon lost to CO2
Gln.Blood_2 <- l.EMU.sim $ `13CLino` $Gln.Blood_2[2]
Gln.Blood_3 <- l.EMU.sim $ `13CLino` $Gln.Blood_3[2]
Gln.Blood_4 <- l.EMU.sim $ `13CLino` $Gln.Blood_4[2]
Gln.Blood_5 <- l.EMU.sim $ `13CLino` $Gln.Blood_5[2]
Gln.Blood_mean <- mean(Gln.Blood_1, Gln.Blood_2, Gln.Blood_3, Gln.Blood_4, Gln.Blood_5)
Gln.Blood_2345 <- mean(Gln.Blood_2, Gln.Blood_3, Gln.Blood_4, Gln.Blood_5)

# CO2 labeling (assume proportional to blood linoleate labeling)
CO2.Blood <- l.EMU.sim $ `13CLino` $CO2.Blood_1[2]

# blood linoleate labeling
Lino.Blood_1 <- l.EMU.sim $ `13CLino` $Lino.Blood_1[2]
Lino.Blood_2 <- l.EMU.sim $ `13CLino` $Lino.Blood_2[2]
Lino.Blood_mean <- mean(Lino.Blood_1, Lino.Blood_2)


# # 3HB labeling
# HB.Blood_1 <- l.EMU.sim $ `13CLino` $HB.Blood_1 [2]
# HB.Blood_2 <- l.EMU.sim $ `13CLino` $HB.Blood_2 [2]
# HB.Blood_3 <- l.EMU.sim $ `13CLino` $HB.Blood_3 [2]
# HB.Blood_4 <- l.EMU.sim $ `13CLino` $HB.Blood_4 [2]
# HB.Blood_mean <- mean(HB.Blood_1, HB.Blood_2, HB.Blood_3, HB.Blood_4)


# ---------<>---------<>---------<>---------<>---------<>---------<>---------<>---------<>---------<>---------<>---------


# # Option 1: from the full model
# # fluxes (molecules of products produced, nmol /min/ g body weight)
# PEPCK <- 447
# CS <- 391
# ME <- 2 
# PEPCK.CS.ME <- PEPCK + CS + ME
# PC <- 444 
# PDH <- 13
# GLS <- 26



# # Option 2: from the reduced model, leading to very similar results
# fluxes (molecules of products produced, nmol /min/ g body weight)
# load enzyme index
d.enz <- read_excel(path = "../../LvM_Fasted/LvM_fasted.xlsx") %>%
  filter(!is.na(enzyme)) %>%
  select(enzyme, reactions, Marker)

optimal_keyFlux <- tibble(flux = c(env$optimal_allFlux)) %>%
  mutate(enzyme = 1:nrow(.)) %>% # flux created
  full_join(d.enz, by = "enzyme") %>%  # augment with name of enzymes
  filter(!is.na(Marker))

PEPCK <- optimal_keyFlux[optimal_keyFlux$Marker == "PEPCK", ]$flux ; PEPCK
CS    <- optimal_keyFlux[optimal_keyFlux$Marker == "CS", ]   $flux ; CS
IDH   <- CS
ME    <- optimal_keyFlux[optimal_keyFlux$Marker == "ME", ]   $flux ; ME
PEPCK.CS.ME <- PEPCK + CS + ME

PC    <- optimal_keyFlux[optimal_keyFlux$Marker == "PC", ]   $flux ; PC
PDH   <- optimal_keyFlux[optimal_keyFlux$Marker == "PDH", ]  $flux ; PDH
GLS   <- optimal_keyFlux[optimal_keyFlux$Marker == "GLS", ]  $flux ; GLS

KETOout  <- optimal_keyFlux[optimal_keyFlux$Marker == "KETOout", ]  $flux ; KETOout
KETOin   <- optimal_keyFlux[optimal_keyFlux$Marker == "KETOin", ]   $flux ; KETOin
KETO  <- KETOout - KETOin # net ketone production


OGDC <- optimal_keyFlux[optimal_keyFlux$Marker == "OGDC", ]  $flux ; OGDC

# <>---------<>---------<>---------<>---------<>---------<>---------<>---------<>---------<>---------<>---------<>---------






# calculate the time point-specific partial-molecular average labeling
d.pulse.chase.atomFlux <- d.pulse.chase %>% 
  
  # here we are using the ratio of lebeling of different EMU as reference
  # extrapolate to the pulse chase condition
  mutate(
    Glutamine_2345 = Gln.Blood_2345 / Gln.Blood_mean * Glutamine,
    CO2 = CO2.Blood / Lino.Blood_mean * `C18:2`, 
    Pyr = Lactate,                           # in pulse-chase, assume pyruvate labeling's same as lactate
    Pyr_23 = Pyr.Lv_23 / Pyr.Lv.mean * Pyr,  # use MFA to get the EMU23/full carbon average as benchmark for the same pyruvate compound
    aKG_1 = aKG.Lv_1 / Mal.Lv.mean  * Malate # use MFA to get aKG1/mal average as benchmark 
  ) %>% 
  
  # # mutate(b = 2 * Malate * (PEPCK.CS.ME) - PDH * Pyr_23 - 1.5 * PC * Pyr - .5 * PC * CO2 - 2 * GLS * Glutamine_2345 )
  # mutate(b = 2 * (PEPCK * Malate + ME * Malate + KETO * `3-Hydroxybutyric acid` ) +
  #          .5 * (IDH * Cit.Lv_6 + OGDC * aKG_1) -
  #          (# 1.5 * PC * Pyr  +
  #             .5 * PC * CO2  +
  #             PDH * Pyr_23  +
  #             2.5 * GLS * Glutamine
  #          )
  # )
  
  mutate(b = 2 * (PEPCK * Malate + ME * Malate + KETO * 0 ) +
           .5 * (IDH * 0 + OGDC * 0) -
           (# 1.5 * PC * Pyr  +
             .5 * PC * 0  +
               PDH * 0  +
               2.5 * GLS * 0
           )
  )



# d.pulse.chase.atomFlux %>% select(Malate, Lactate, `C18:2`, `TAG C18:2`) %>% pairs()


lm(data = d.pulse.chase.atomFlux %>% filter(method == "pulse"), formula = b ~ `C18:2` + `TAG C18:2` + 0) %>% summary()
lm(data = d.pulse.chase.atomFlux %>% filter(method == "chase"), formula = b ~ `C18:2` + `TAG C18:2` + 0) %>% summary()

mdl <- lm(data = d.pulse.chase.atomFlux, formula = b ~ `C18:2` + `TAG C18:2` + 0); summary(mdl)

# get summary table
mdl.summary <- summary(mdl)$coefficients ; mdl.summary  # print the standard error of the coefficient
d.mdl.summary <- mdl.summary %>% as_tibble() %>%  mutate(metab = rownames(mdl.summary), .before = 1) # add row names as columns

# get flux
NEFA_lino_kin <- d.mdl.summary %>% filter(metab == "`C18:2`")     %>% pull(Estimate) ; NEFA_lino_kin
TAG_lino_kin  <- d.mdl.summary %>% filter(metab == "`TAG C18:2`") %>% pull(Estimate) ; TAG_lino_kin

# get error
NEFA_lino_kin_err <- d.mdl.summary %>% filter(metab == "`C18:2`")     %>% pull(`Std. Error`) ; NEFA_lino_kin_err
TAG_lino_kin_err  <- d.mdl.summary %>% filter(metab == "`TAG C18:2`") %>% pull(`Std. Error`) ; TAG_lino_kin_err



# Check fitting
d.pulse.chase.atomFlux.fitted <- d.pulse.chase.atomFlux %>% 
  mutate(fitted = mdl$fitted.values)

d.pulse.chase.atomFlux.fitted %>% 
  ggplot(aes(x = b, y = fitted, color = method)) + 
  geom_point(size = 3) +
  geom_abline(slope = 1, intercept = 0) +
  # coord_cartesian(xlim = c(0, NA), ylim = c(0, NA)) +
  scale_x_continuous(expand = expansion(mult = c(0.1, .1))) +
  scale_y_continuous(expand = expansion(mult = c(0.1, .1))) 





# Get more liver nutrient burning fluxes from the LvM MFA model 
# flux in molecules of products produced, nmol /min/ g body weight

# FFA oxidation directly from MFA, before integration with the kinetics tracing
NEFA_Palm <- optimal_keyFlux[optimal_keyFlux$Marker == "NEFA-Palm", ]  $flux ; NEFA_Palm
NEFA_Ole  <- optimal_keyFlux[optimal_keyFlux$Marker == "NEFA-Ole", ]   $flux ; NEFA_Ole
NEFA_Lino <- optimal_keyFlux[optimal_keyFlux$Marker == "NEFA-Lino", ]  $flux ; NEFA_Lino
carbs     <- optimal_keyFlux[optimal_keyFlux$Marker == "PDH", ]        $flux ; carbs



# integrate with the kinetics tracing: 
# use the ratio of NEFA to linoleate in MFA to calculate the expected kinetic - derived NEFA and TAG fluxes

# NEFA flux
NEFA_Palm_kin <- NEFA_Palm / NEFA_Lino * NEFA_lino_kin; NEFA_Palm_kin
NEFA_Ole_kin  <- NEFA_Ole  / NEFA_Lino * NEFA_lino_kin; NEFA_Ole_kin

# TAG flux
TAG_Palm_kin <-  NEFA_Palm / NEFA_Lino * TAG_lino_kin; TAG_Palm_kin
TAG_Ole_kin  <-  NEFA_Ole  / NEFA_Lino * TAG_lino_kin; TAG_Ole_kin

# NEFA error 
NEFA_Palm_kin_err <- NEFA_Palm / NEFA_Lino * NEFA_lino_kin_err; NEFA_Palm_kin_err
NEFA_Ole_kin_err  <- NEFA_Ole  / NEFA_Lino * NEFA_lino_kin_err; NEFA_Ole_kin_err

# TAG error 
TAG_Palm_kin_err <- NEFA_Palm / NEFA_Lino * TAG_lino_kin_err; TAG_Palm_kin_err
TAG_Ole_kin_err  <- NEFA_Ole  / NEFA_Lino * TAG_lino_kin_err; TAG_Ole_kin_err


# get PDH error from the CI data
load("../../LvM_Fasted/6_CI_fasted.RData")
d.CI.sem <- d.CI.bounds.reactions.reduced.LvM %>%  # get the C.I dataset
  select(flux.index, R, L, flux.optimal, reactions) %>% 
  mutate(sem = (R-L) / 2 / 1.98)

d.CI.sem

# PDH error
carbs_err <- d.CI.sem %>% filter(reactions == "Pyr.Lvm->AcCoA.Lv+CO2.Blood") %>% pull(sem)

# Ketogenesis net flux
keto.net       <- d.CI.sem %>% filter(reactions == "(AcAct.Lv->HB.Blood) - (HB.Blood->AcAct.Lv)") %>% pull(flux.optimal)
keto.net_error <- d.CI.sem %>% filter(reactions == "(AcAct.Lv->HB.Blood) - (HB.Blood->AcAct.Lv)") %>% pull(sem)

CS_error  <- d.CI.sem %>% filter(reactions == "OAA.Lv+AcCoA.Lv->Cit.Lv") %>% pull(sem)
PDH_error <- d.CI.sem %>% filter(reactions == "Pyr.Lvm->AcCoA.Lv+CO2.Blood") %>% pull(sem)

# calculate 'others -> AcCoA' flux and error 
others <- (keto.net + CS) - (NEFA_Palm_kin + NEFA_Ole_kin + NEFA_lino_kin) - (TAG_Palm_kin + TAG_Ole_kin + TAG_lino_kin) - PDH; others
others.err <-  (keto.net_error^2 + CS_error^2 + (NEFA_Palm_kin_err^2 + NEFA_Ole_kin_err^2 + NEFA_lino_kin_err^2) + (TAG_Palm_kin_err^2 + TAG_Ole_kin_err^2 + TAG_lino_kin_err^2) + PDH_error^2 ) %>% sqrt()



# plot contribution sources to liver TAG
d.liver.energy <- tibble(
  source = c( "NEFA-Palm",        "NEFA-Ole",          "NEFA-Lino",         "TAG-Palm",         "TAG-Ole",         "TAG-Lino",         "others",     "carbs"),
  flux   = c( NEFA_Palm_kin,       NEFA_Ole_kin,        NEFA_lino_kin,       TAG_Palm_kin,       TAG_Ole_kin,       TAG_lino_kin,       others,       PDH),
  error  = c( NEFA_Palm_kin_err,   NEFA_Ole_kin_err,    NEFA_lino_kin_err,   TAG_Palm_kin_err,   TAG_Ole_kin_err,   TAG_lino_kin_err,   others.err,   PDH_error)) %>% 
  mutate(source = factor(source, levels = .$source, ordered = T)) %>% 
  arrange(desc(source)) %>% 
  mutate(errY = cumsum(flux))


p <- d.liver.energy %>% 
  ggplot(aes(x = 1, y = flux, fill = source, alpha = source)) +
  geom_col(color = "black") +
  # coord_polar(theta = "y") +
  # scale_x_continuous(limits = c(0, 1.5)) +
  # theme_void() +
  scale_fill_manual(values = c(
    "carbs" = "tomato",
    "NEFA-Palm" = "purple", "NEFA-Ole" = "orange", "NEFA-Lino" = "skyblue2", 
    "TAG-Palm" = "purple", "TAG-Ole" = "orange", "TAG-Lino" = "skyblue2", 
    "others" = "grey")
  ) +
  geom_errorbar(aes(ymax = errY, ymin = errY - error), width = .2, alpha = 1) +
  scale_alpha_manual(values = c(1, 1, 1, .2, .2, .2, 1, 1)) +
  scale_y_continuous(expand = expansion(mult = c(0, .1)), breaks = seq(0, 1000, 50)) +
  scale_x_continuous(expand = expansion(add = .3)) +
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text = element_text(colour = "black"),
        axis.title.x = element_blank())
p

ggsave("liver energy source.pdf", height = 4, width = 2)



# add label
p +
  geom_text(aes(label = paste(source, " - ", round(flux/sum(flux) * 100, 1), "%"), 
                y = flux, 
                x = .7), 
            position = position_stack(vjust = .5), 
            color = "black", alpha = 1, hjust = 0)  





# calculate 3HB contribution from circulating FAs and TAG
d.chase.3HB <- d.chase %>% 
  filter(sample %>% str_detect("SAC")) %>% 
  filter(Compound %in% c("C18:2", "TAG C18:2", "3-Hydroxybutyric acid")) %>% ungroup() %>% 
  mutate(mouse = str_extract(sample, "\\d$")) %>% 
  select(-c(sample, tissue)) %>% 
  pivot_wider(names_from = Compound, values_from = labeling) %>% 
  mutate(method = "chase")


d.pulse.3HB <- d.pulse %>% 
  filter(Compound %in% c("C18:2", "TAG C18:2", "3-Hydroxybutyric acid")) %>% ungroup() %>% 
  select(-tissue) %>% 
  pivot_wider(names_from = Compound, values_from = labeling) %>% 
  select(-batch) %>% relocate('time (h)') %>% 
  mutate(method = "pulse")
d.pulse.3HB <- d.pulse.3HB[complete.cases(d.pulse.3HB)  , ]

# stack chase - pulse data
d.pulse.chase.3HB <- rbind(d.pulse.3HB, d.chase.3HB)

d.pulse.chase.3HB

mdl.3HB <- lm(data = d.pulse.chase.3HB, `3-Hydroxybutyric acid` ~ `TAG C18:2` + `C18:2` + 0)

d.pulse.chase.3HB %>% 
  mutate(fitted = mdl.3HB$fitted.values) %>% 
  ggplot(aes(x = `3-Hydroxybutyric acid`, y = fitted, color = method)) + geom_point(size = 3) +
  geom_abline(intercept = 0, slope = 1) +
  coord_equal(xlim = c(0, NA), ylim = c(0, NA)) + theme(legend.position = "right")



# others vs top 3 fatty acids
(20+41+54+46+91+120) / 74

