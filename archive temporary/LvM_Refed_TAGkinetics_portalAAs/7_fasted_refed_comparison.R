rm(list = ls())
library(tidyverse)


myClassic <- theme_classic(base_size = 14) +
  theme(axis.title.y = element_text(margin = margin(r = 5, unit = "pt")),
        axis.title.x = element_text(margin = margin(t = 5, unit = "pt")),
        axis.title = element_text(face = "bold"),
        strip.background = element_blank(),
        strip.text = element_text(size = 15, face = "bold"))



rstudioapi::getActiveDocumentContext()$path %>% dirname() %>% setwd(); getwd()
d.fasted <- load("6_CI_fasted.RData")             %>% get() %>% mutate(state = "fasted", .before = 1)
d.refed  <- load("../LvM_Refed/6_CI_refed.RData") %>% get() %>% mutate(state = "refed",  .before = 1)

d.CI.all <- bind_rows(d.fasted, d.refed)


# normalize fluxes to CS flux
d.norm.basis <- d.CI.all %>% 
  filter(reactions == "OAA.Lv+AcCoA.Lv->Cit.Lv") %>% 
  select(state, flux.optimal) %>% 
  rename(norm.basis = flux.optimal)



d.CI.norm <- d.CI.all %>% 
  left_join(d.norm.basis, by = "state") %>% 
  group_by(state) %>% 
  mutate(R.norm = R / norm.basis * 100, 
         L.norm = L / norm.basis * 100, 
         flux.optimal.norm = flux.optimal / norm.basis * 100)


# plot relative fluxes, normalized to CS as 100 
d.CI.norm %>% 
  filter(flux.index %in% c(12:19, 24:28, 30, 32:35, 66, 68)) %>% 
  ggplot(aes(x = reactions, y = flux.optimal.norm, color = state)) +
  # geom_point(position = position_dodge(.3),) +
  geom_errorbar(aes(ymin = L.norm, ymax = R.norm),
                position = position_dodge(.3),
                width = .2) +
  geom_text(
    aes(label = paste(round(L.norm), "-", round(R.norm)),
        y = R.norm + 2),
    position = position_dodge(.3),
    size = 2.5, fontface = "bold") +
  myClassic +
  theme(axis.text.x = element_text(angle = 60, hjust = 1),
        axis.text = element_text(size = 14)) +
  scale_y_continuous(breaks = seq(0, 200, 10),
                     expand = expansion(mult = c(0, .1))) +
  labs(y = "relative fluxes normalized to CS") +
  scale_color_manual(values = c("fasted" = "steelblue", "refed" = "tomato")) 
  # geom_vline(xintercept = seq(1, 20, 2),
  #            linewidth = 13, alpha = .05)





# plot absolute fluxes
d.CI.norm %>% 
  filter(flux.index %in% c(12:19, 24:28, 30, 32:35, 66, 68)) %>% 
  ggplot(aes(x = reactions, y = flux.optimal, color = state)) +
  # geom_point(position = position_dodge(.3),) +
  geom_errorbar(aes(ymin = L, ymax = R),
                position = position_dodge(.3),
                width = .2) +
  geom_text(
    aes(label = paste(round(L), "-", round(R)),
        y = R + 2),
    position = position_dodge(.3),
    size = 2.5, fontface = "bold") +
  myClassic +
  theme(axis.text.x = element_text(angle = 60, hjust = 1),
        axis.text = element_text(size = 14)) +
  scale_y_continuous(breaks = seq(0, 2000, 100),
                     expand = expansion(mult = c(0, .1))) +
  labs(y = "nmol molecules / g / BW") +
  scale_color_manual(values = c("fasted" = "steelblue", "refed" = "tomato")) 
# geom_vline(xintercept = seq(1, 20, 2),
#            linewidth = 13, alpha = .05)









# bar plot comparing PC vs PDH
r.ordered <- c("OAA.Lv+AcCoA.Lv->Cit.Lv", "Pyr.Lvm+CO2.Blood->OAA.Lv", "Pyr.Lvm->AcCoA.Lv+CO2.Blood")

# plot fluxes, normalized to CS as 100
d.CI.norm %>% 
  filter(flux.index %in% c(12, 17, 25)) %>% 
  mutate(reactions = factor(reactions, levels = r.ordered)) %>% 
  
  ggplot(aes(x = reactions, y = flux.optimal.norm, fill = reactions)) +
  geom_col(position = position_dodge(.9), color = "black", alpha = .6, size = 1) +
  geom_errorbar(aes(ymin = L.norm, ymax = R.norm),
                position = position_dodge(.9),
                width = .2, size = 1) +
  # geom_text(
  #   aes(label = paste(round(L.norm), "-", round(R.norm)),
  #       y = R.norm + 2),
  #   position = position_dodge(.3),
  #   size = 2.5, fontface = "bold") +
  myClassic +
  theme(axis.text.x = element_text(angle = 60, hjust = 1),
        axis.text = element_text(size = 17)) +
  scale_y_continuous(breaks = seq(0, 200, 20),
                     expand = expansion(mult = c(0, .1))) +
  labs(y = "nmol molecules / g / BW") +
  scale_fill_manual(values = c("#1FB050", "#25B0F0", "red")) +
  facet_wrap(~state, scales = "free")
# geom_vline(xintercept = seq(1, 20, 2),
#            linewidth = 13, alpha = .05)




