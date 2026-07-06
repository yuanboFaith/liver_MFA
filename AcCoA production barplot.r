library(tidyverse)
library(readxl)


theme_set(
  theme_classic(base_size = 24) +
    theme(legend.position = "none", 
          # panel.border = element_rect(colour = "black", fill = NA, linewidth = .3),
          strip.background = element_blank(),
          axis.title.y = element_text(margin = margin(r = 10, unit = "pt")),
          axis.title = element_text(face = "bold", colour = "black"),
          axis.text = element_text(colour = "black"),
          axis.text.y = element_text(size = 28),
          strip.text = element_text(size = 12, face = "bold"),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),
          axis.title.x = element_blank()) 
)


(rstudioapi::getActiveDocumentContext())$path %>% dirname() %>% setwd(); getwd()


# Fasted TAG kinetics
path.compiledCI <- "CI compiled.xlsx"
d.Fast.TAGkinetics <- read_excel(path.compiledCI, sheet = "Fasted TAGkinetics", range = "P1:R11") 

# arrange in order
d.Fast.TAGkinetics <- d.Fast.TAGkinetics %>% 
  mutate(fuels = factor(fuels, levels = .$fuels, ordered = T)) %>% 
  arrange(desc(fuels)) %>% 
  mutate(err.Y = cumsum(flux.center)) %>% 
  
  # calculate percentage
  mutate(contri         = flux.center / sum(flux.center) * 100) %>% 
  mutate(contri.SEM     = SEM         / sum(flux.center) * 100) %>% 
  mutate(contri.error.Y = cumsum (contri))
d.Fast.TAGkinetics



p.Fast.TAGkinetics <- d.Fast.TAGkinetics %>% 
  ggplot(aes(x = 1, y = contri, fill = fuels, alpha = fuels)) +
  geom_col(color = "black") +
  scale_fill_manual(values = c(
    "NEFA-Palm" = "purple", "NEFA-Ole" = "orange", "NEFA-Lino" = "skyblue2",
    "TAG-Palm" = "purple", "TAG-Ole" = "orange", "TAG-Lino" = "skyblue2",
    "other NEFAs" = "grey80", 
    "other TAGs" = "grey80",
    "carbs" = "tomato",
    "others" = "green3")) +
  geom_errorbar(aes(ymax = contri.error.Y, ymin = contri.error.Y - contri.SEM), width = .2, alpha = 1) +
  scale_alpha_manual(values = c(1, 1, 1, 1, .2, .2, .2, .2, 1, .2)) +
  scale_y_continuous(expand = expansion(mult = c(0, .1)), breaks = seq(0, 100, 20)) +
  scale_x_continuous(expand = expansion(add = .3)) 


p.Fast.TAGkinetics

ggsave("figures/barplot_Fasted_TAGkinetics.pdf", height = 5*1.6, width = 3*1.3)

# With label
p.Fast.TAGkinetics + 
  geom_text(aes(label = fuels), position = position_stack(vjust = .5), alpha = 1, size = 2) 

ggsave("figures/barplot_Fasted_TAGkinetics_labeled.pdf", height = 5, width = 3)







# -$--$--$--$--$--$--$--$--$--$--$--$--$--$--$--$--$--$--$--$--$--$--$--$--$--$--$--$--$--$--$-







# Refed + TAG kinetics

d.Refed.TAGkinetics <- read_excel(path.compiledCI, sheet = "Refed TAGkinetics", range = "P1:R11") 

# arrange in order
d.Refed.TAGkinetics <- d.Refed.TAGkinetics %>% 
  
  mutate(fuels = factor(fuels, levels = .$fuels, ordered = T)) %>% 
  arrange(desc(fuels)) %>% 
  mutate(err.Y = cumsum(flux.center)) %>% 
  
  # calculate percentage
  mutate(contri         = flux.center / sum(flux.center) * 100) %>% 
  mutate(contri.SEM     = SEM         / sum(flux.center) * 100) %>% 
  mutate(contri.error.Y = cumsum (contri))

d.Refed.TAGkinetics



p.Refed.TAGkinetics <- d.Refed.TAGkinetics %>% 
  ggplot(aes(x = 1, y = contri, fill = fuels, alpha = fuels)) +
  geom_col(color = "black") +
  scale_fill_manual(values = c(
    "NEFA-Palm" = "purple", "NEFA-Ole" = "orange", "NEFA-Lino" = "skyblue2",
    "TAG-Palm" = "purple", "TAG-Ole" = "orange", "TAG-Lino" = "skyblue2",
    "other NEFAs" = "grey80", 
    "other TAGs" = "grey80",
    "carbs" = "tomato",
    "others" = "green3")) +
  geom_errorbar(aes(ymax = contri.error.Y, ymin = contri.error.Y - contri.SEM), width = .2, alpha = 1) +
  scale_alpha_manual(values = c(1, 1, 1, 1, .2, .2, .2, .2, 1, .2)) +
  scale_y_continuous(expand = expansion(mult = c(0, .1)), breaks = seq(0, 100, 20)) +
  scale_x_continuous(expand = expansion(add = .3)) +
  theme(axis.text.y = element_text(size = 24))


p.Refed.TAGkinetics

ggsave("figures/barplot_Refed_TAGkinetics.pdf",  height = 8, width = 3.5)

# With label
p.Refed.TAGkinetics + 
  geom_text(aes(label = fuels), position = position_stack(vjust = .5), alpha = 1, size = 2) 

ggsave("figures/barplot_Refed_TAGkinetics_labeled.pdf",  height = 8, width = 4)





# -$--$--$--$--$--$--$--$--$--$--$--$--$--$--$--$--$--$--$--$--$--$--$--$--$--$--$--$--$--$--$-






# Refed + TAG kinetics + hepatic portal AAs
d.Refed.TAGkinetics_hpAAs <- read_excel(
  path.compiledCI, sheet = "Refed TAGkinetics_portalAAs", range = "T1:V10") 

# arrange in order
d.Refed.TAGkinetics_hpAAs <- d.Refed.TAGkinetics_hpAAs %>% 
  mutate(fuels = factor(fuels, levels = .$fuels, ordered = T)) %>% 
  arrange(desc(fuels)) %>% 
  mutate(err.Y = cumsum(flux.center)) %>% 
  
  # calculate percentage
  mutate(contri         = flux.center / sum(flux.center) * 100) %>% 
  mutate(contri.SEM     = SEM         / sum(flux.center) * 100) %>% 
  mutate(contri.error.Y = cumsum (contri))

d.Refed.TAGkinetics_hpAAs



p.Refed.TAGkinetics.hpAAs <- d.Refed.TAGkinetics_hpAAs %>% 
  ggplot(aes(x = 1, y = contri, fill = fuels, alpha = fuels)) +
  geom_col(color = "black") +
  scale_fill_manual(values = c(
    "NEFA-Palm" = "purple", "NEFA-Ole" = "orange", "NEFA-Lino" = "skyblue2",
    "TAG-Palm" = "purple", "TAG-Ole" = "orange", "TAG-Lino" = "skyblue2",
    # "other NEFAs" = "grey80", 
    # "other TAGs" = "grey80",
    "Amino Acids" = "firebrick4",
    "carbs" = "tomato",
    "others" = "green3")) +
  geom_errorbar(aes(ymax = contri.error.Y, ymin = contri.error.Y - contri.SEM), width = .2, alpha = 1) +
  scale_alpha_manual(values = c(1,  1,  1, # NEFA
                                .2, .2, .2, # TAG 
                                .4,  1, .2  # amino acids, carbs, others   
  )) +
  scale_y_continuous(expand = expansion(mult = c(0, .1)), 
                     breaks = seq(0, 100, 20)) +
  scale_x_continuous(expand = expansion(add = .3)) +
  coord_cartesian(ylim = c(0, NA)) 


p.Refed.TAGkinetics.hpAAs

ggsave("figures/barplot_Refed_TAGkinetics_hpAAs.pdf", height = 7, width = 3.5)

# With label
p.Refed.TAGkinetics.hpAAs + 
  geom_text(aes(label = fuels), position = position_stack(vjust = .5), alpha = 1, size = 2) 

ggsave("figures/barplot_Refed_TAGkinetics_hpAAs_labeled.pdf", height = 7, width = 3.5)


