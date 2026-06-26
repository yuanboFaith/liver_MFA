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


rstudioapi::getActiveDocumentContext()$path %>% dirname() %>% setwd()


myPath <- "pulse data.xlsx"

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
d.labeling <- d.labeling %>% filter(! (tissue == "serum" & batch == 2))




# calculate averaged carbon labeling
d.average.C.labeling <- d.labeling %>%
  group_by(batch, mouse, `time (h)`, tissue, Compound) %>%
  mutate(labeling.weighted = C_Label / max(C_Label) * labeling) %>%
  summarise(labeling = sum(labeling.weighted)) %>%
  
  # remove failed animal
  filter(mouse != "G") # serious leaking from the red pinport found at the end of infusion




# plot labeling pulse kinetics
func.plt <- function(md = d.average.C.labeling){
  md %>% 
    ggplot(aes(x = `time (h)`, y = labeling, color = Compound)) +
    geom_point(size = 2) +
    # geom_text(aes(label = mouse), size = 5) +
    facet_wrap(~Compound, scales = "free") +
    scale_x_continuous(expand = expansion(mult = c(0, .03)), breaks = seq(0, 12, 2))  +
    scale_y_continuous(expand = expansion(mult = c(0, .03)), n.breaks = 7,
                       labels = ~.x * 100,
                       name = "labeling (%)")  +
    coord_cartesian(xlim = c(0, NA), ylim = c(0, NA)) +
    stat_smooth(
      method = "nls", 
      formula = y ~ a * (1 - exp(-b * x)),
      method.args = list(start = list(a = .01, b = .1)),
      se = FALSE,
      fullrange = T
      # color = "blue"
    ) 
}




# source labeling
set.seed(1234)
# myColors <- (c(RColorBrewer::brewer.pal(6, "Paired")) %>% colorRampPalette())(7) %>% sample() %>% rev()

myColors <- c("TAG C18:2" = "firebrick", "C18:2" = "steelblue", 
              "3-Hydroxybutyric acid" = "purple", 
              "Glucose" = "orange", "Lactate" = "green4", 
              "Glutamine" = "pink3", "Alanine" = "grey50")

d.average.C.labeling$Compound %>% unique()

sources <- c("TAG C18:2", "C18:2", "Glucose", "Lactate", "Glutamine", "Alanine", "3-Hydroxybutyric acid")


p1 <- d.average.C.labeling %>% 
  filter(Compound %in% sources) %>% 
  func.plt() +
  facet_wrap(~"fuel sources") +
  theme(legend.position = "right") +
  # scale_color_brewer(palette = "Spectral") +
  scale_color_manual(values = myColors) 
p1


targets <- c("Succinic acid", "Malate")


p2 <- d.average.C.labeling %>% 
  filter(Compound %in% targets) %>% 
  func.plt() +
  facet_wrap(~"liver TCA") +
  theme(legend.position = "right")
p2

cowplot::plot_grid(p1 + theme(legend.position = "none"), 
                   ggplot() + theme_void(),
                   p2 + theme(legend.position = "none"), 
                   nrow = 1, 
                   rel_widths = c(1, .1, 1))

ggsave(filename = "pulse.pdf", height = 3.6, width = 7)



# # calculate the contribution fraction 
# 
# # Note!!! This calculation still UNDERESTIMATE the fatty acids/TAG contribution, 
# # because PC flux dilutes the labeling, but contributes with no carbon fuel supply
# # but it's good to take a first look using this method, before proceeding to MFA 
# 
# d.Q <- d.average.C.labeling %>% 
#   filter( (tissue == "liver" & Compound %in% c("TAG C18:2", "Malate", "Succinic acid")) |
#             tissue == "serum" & Compound %in% c("C18:2")) %>% 
#   ungroup() %>% 
#   filter(mouse %in% c("A", "B", "C", "D", "E", "F", "G", "H", "I"))
# 
# # source labeling
# d.source <- d.Q %>% filter( (tissue == "liver" & Compound %in% c("TAG C18:2")) |
#                               tissue == "serum" & Compound %in% c("C18:2")) 
# d.source.spread <- d.source %>% select(-tissue) %>% spread(Compound, labeling)
# 
# 
# # TCA labeling
# d.TCA <- d.Q %>% filter( (tissue == "liver" & Compound %in% c("Malate", "Succinic acid")) ) %>% 
#   # take the average labeling 
#   group_by(mouse, `time (h)`) %>%  
#   summarise(TCA.labeling = mean(labeling))
# 
# # combine the source and TCA label, and choose complete data set
# d.label <- d.source.spread %>% left_join(d.TCA, by = c("mouse", "time (h)"))
# d.label.complete <- d.label[complete.cases(d.label), ] %>% arrange(`time (h)`)
# 
# # fit a linear regression
# mdl <- lm(data = d.label.complete, formula = TCA.labeling ~ `C18:2` + `TAG C18:2` )
# mdl
# 
# plot(d.label.complete$TCA.labeling,  mdl$fitted.values, xlab = "observed", ylab = "fitted", main = "fitted vs. observed")
# abline(a = 0, b = 1, col = "red", lty = 2)  # Red dashed reference line
# 
# 
# # force contribution to be non-zero
# A <- as.matrix(d.label.complete[, c("C18:2", "TAG C18:2")])
# b <- d.label.complete$TCA.labeling
# 
# # Non-negative least squares fit (without intercept)
# fit <- nnls::nnls(A, b)
# coef(fit)
# plot(fit$fitted, b); abline(a = 0, b = 1, col = "red", lty = 2)  # Red dashed reference line
# 
# 
# # bootstrap to estimate the error
# set.seed(42)
# B <- 1000 # number of bootstrap samples
# 
# coefs.Boot <- matrix(NA, nrow = B, nco = ncol(A))
# 
# for (i in 1:B) {
#   rows <- sample(1:nrow(A), nrow(A), replace=TRUE)
#   A.Boot <- A[rows, ]
#   b.Boot <- b[rows]
#   fit.Boot <- nnls::nnls(A.Boot, b.Boot)
#   coefs.Boot[i, ] <- coef(fit.Boot)
# }
# 
# # Standard error estimates for each coefficient
# error.Boot <- apply(coefs.Boot, 2, sd)
# 
# 
# 
# # plot the contribution
# x1 <- tibble(names = c("C18:2", "C18:2"),
#              frac = coef(fit),
#              error = error.Boot) %>% 
#   mutate(where = c("serum", "TAG"))
# x1  
# 
# # fractions relative to linoleate is from 8_3_confidence_analysis, based on U-13C fatty acid tracing result
# x2 <- x1 %>% mutate(names = "C18:1", frac = frac * .8, error = error * .8);  x2 # oleate C18:1; 
# x3 <- x1 %>% mutate(names = "C16:0", frac = frac * .4,  error = error * .4); x3 # plamitate C16:0
# x4 <- x1 %>% mutate(names = "other", frac = frac * (1 + .8 + .4) * 25/75,
#                     error = error * (1 + .8 + .4) * 25/75); x3 # other fatty acid, 30% of TOTAL fatty acids, or 
# 
# d.contri.f <- bind_rows(x1, x2, x3, x4) 
# 
# # calculate error bar position
# d.contri.f <- d.contri.f %>% 
#   mutate(names = fct_reorder(names, -frac)) %>% 
#   arrange(names) %>% 
#   group_by(where) %>% 
#   mutate(y.err = cumsum(frac))
# 
# d.contri.f %>% 
#   ggplot(aes(x = where, y = frac, fill = names)) +
#   geom_col(color = "black", position = position_stack(reverse = T),
#            alpha = c(1, 1, .5, .5, .5, .5, .3, .3)) +
#   geom_errorbar(aes(ymin = y.err  - error,
#                     ymax = y.err), 
#                 width = .15) +
#   scale_y_continuous(expand = expansion(mult = c(0, .1)),
#                      labels = ~. * 100,
#                      breaks = seq(0, 1, .05)) +
#   theme(legend.position = "right") +
#   labs(y = "contribution % to TCA", x = NULL) +
#   scale_fill_brewer(palette = "Pastel1")
# 
# ggsave("contribution Pulsing.pdf", height = 4/1.2, width = 4)
# 
# 
# d.contri.f$frac %>% sum()
# 
# 
# 
# # # save project
d.pulse <- d.average.C.labeling

save(myColors, sources, targets, d.pulse, file =  "liver_kinetics_pulse.RData")

