# Analyse repeated run results

rm(list = ls())
library(tidyverse)

# set path to the folder of the current active script
rstudioapi::getActiveDocumentContext()$path %>% dirname() %>% setwd(); getwd()


myClassic <- theme_classic(base_size = 14) +
  theme(axis.title.y = element_text(margin = margin(r = 5, unit = "pt")),
        axis.title.x = element_text(margin = margin(t = 5, unit = "pt")),
        axis.title = element_text(face = "bold"),
        strip.background = element_blank(),
        strip.text = element_text(size = 15, face = "bold"))




# Compile the cost from all repeated runs
files.repeated.run <- list.files("./repeat_results",  full.names = TRUE, pattern = "\\.RData$")

d.repeat.cost <- tibble(.rows = 0)
d.repeat.u.initial <- tibble(.rows = 0)
d.repeat.u.initial.converging <- tibble(.rows = 0)

for (file.repeat.i in files.repeated.run) {
  # file.repeat.i=files.repeated.run[1]
  message("loading file: ", file.repeat.i, "...")
  
  load(file.repeat.i) # care was taken not to overwrite the 'file.repeat.i' when loading 'file.repeat.i'
  
  # 1. cost
  d.repeat.cost.i <- tibble(
    cost.iterations = cost.iterations, 
    file = file.repeat.i %>% basename() %>% str_remove(".RData")) %>% 
    mutate(ite = 1:nrow(.))
  
  # combine repeats
  d.repeat.cost <- bind_rows(d.repeat.cost, d.repeat.cost.i)
  
  
  
  # 2. initial values
  d.repeat.u.initial.i <- tibble(
    flux = names(u.initial),
    u.initial = u.initial,
    file = file.repeat.i %>% basename() %>% str_remove(".RData"))
  # combine repeats
  d.repeat.u.initial <- bind_rows(d.repeat.u.initial, d.repeat.u.initial.i)
  
  
  
  # 3. fluxes during iteration
  d.repeat.u.initial.converging.i <-  d.u.iterations %>% 
    mutate(file = file.repeat.i %>% basename() %>% str_remove(".RData"))
  
  # combine repeats
  d.repeat.u.initial.converging <- bind_rows(d.repeat.u.initial.converging, d.repeat.u.initial.converging.i)
}

# make file in order
d.repeat.cost <- d.repeat.cost %>% 
  mutate(repeat.index = str_extract(file, "\\d{1,2}$") %>% as.integer()) %>% 
  arrange(repeat.index) %>% 
  mutate(file = factor(file, levels = .$file %>% unique()))

files.ordered <- d.repeat.cost$file %>% unique()


# color of repeats
set.seed(1)
color.repeats <- colorRampPalette(c(RColorBrewer::brewer.pal(8, "Set1"), "black"))(n_distinct(d.repeat.cost$file)) %>% sample()
names(color.repeats) <- unique(d.repeat.cost$file)


# cost plot 1
d.repeat.cost %>% 
  ggplot(aes(x = ite, y = cost.iterations, color = file, fill = file)) +
  # geom_point(size = 1) + 
  geom_line() +
  scale_x_continuous(breaks = seq(1, 21, 2)) +
  # scale_y_log10(breaks = c(seq(.1, 1, .1), 1:10, seq(10, 100, 10)), expand = expansion(add = 0.1)) + 
  scale_y_log10(breaks = c(.5, 1, 5, 10, 50, 100, 500, 1000)) + 
  annotation_logticks(sides = "l") +
  myClassic +
  # annoate with the terminal point
  geom_point(data = d.repeat.cost %>% group_by(file) %>%  filter(ite == max(ite)),
             shape = 23, size = 2, color = "black",
             position = position_jitter(.2, 0.000005, seed = 1234)) +
  # # annoate with the repeat index
  ggrepel::geom_text_repel(
    data = d.repeat.cost %>% group_by(file) %>%  filter(ite == max(ite)),
    aes(label = file %>% str_extract("\\d{1,2}")),
    color = "black", size = 1.5, box.padding = unit(2, "pt"), max.overlaps = 7) +
  scale_color_manual(values = color.repeats) +
  scale_fill_manual(values = color.repeats) +
  theme(legend.position = "none")


ggsave("./plots/MFA repeat convergence.pdf", height = 3.5*1.1, width = 4.5*1.1)




# initial flux values
d.repeat.u.initial %>% 
  mutate(file = factor(file, levels = files.ordered)) %>% 
  ggplot(aes(x = file, y = u.initial, color = file, fill = file)) +
  ggbeeswarm::geom_quasirandom(width = .2, shape = 21, size = .1, alpha = .7) +
  # geom_line(aes(group = flux), linewidth = .2) +
  scale_y_log10(breaks = 10^c(-5:5),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  annotation_logticks(sides = "l") +
  labs(y = "initial fluxes", x = NULL) +
  myClassic +
  theme(legend.position = "none") +
  scale_x_discrete(labels = function(x)str_extract(x, "\\d{1,2}$")) +
  scale_color_manual(values = color.repeats) +
  scale_fill_manual(values = color.repeats) 

ggsave("./plots/MFA repeat initial fluxes.pdf", height = 4/1.3, width = 5/1.3)




# select fluxes with the minimum cost
d.repeat.best <- d.repeat.cost %>% 
  group_by(repeat.index) %>% 
  filter(ite == max(ite)) %>% # for each repeat, select the last data points
  ungroup() %>% 
  filter(cost.iterations == min(cost.iterations)) %>% # select the repeat with minimal cost
  slice_head(n = 1) # in case there are two repeats with identical cost, select the first one

file.repeat.best <- str_c("./repeat_results/", d.repeat.best$file, ".RData")
load(file.repeat.best) # overwrite the RData loaded in the last loop with the best fit


save(
  optimal_cost, optimal_free_flux, optimal_allFlux, l.EMU.sim,
  file = "5_2_optimal_solution.RData"
)


# Log up the best repeat index number
# sink("./plots/converging_history/best converge.txt", append = T)
# sink()
mmmp <- paste0("Congrats! Selected the best repeat with the minimal cost: ", d.repeat.best$file)
message("Congrats! Selected the best repeat with the minimal cost: ", d.repeat.best$file)

cat("\n\n", mmmp, "\n", file = "./plots/converging_history/log converge.txt", append = TRUE)





d.obs.simu$metabolite %>% unique()

ordered.metabolite <- c("Glc.Blood", "Lac.Blood", "Ala.Blood", "Gln.Blood", "Glycerol.Blood", "HB.Blood",
                        "Palm.Blood", "Ole.Blood", "Lino.Blood", "Mal.Lv", "Suc.Lv", "aKG.Lv")


d.obs.simu <- d.obs.simu %>% mutate(metabolite = factor(metabolite, levels = ordered.metabolite))



# plot simulated vs. observed
# sim vs. obs plot2
p <- d.obs.simu %>% 
  mutate(metabolite = str_replace(metabolite, "\\.Blood", "\\.Bld")) %>% 
  ggplot(aes(x = obs, y = sim, color = tracer)) +
  # geom_text(aes(label = label)) +
  geom_point(size = 3, shape = 21, stroke = .3) +
  # geom_text(aes(label = label), size = 3 , fontface = "bold") +
  scale_x_continuous(transform = "sqrt") +
  scale_y_continuous(transform = "sqrt") +
  # facet_grid(metabolite ~ tracer) +
  # facet_wrap(~tracer) +
  geom_abline(slope = 1, intercept = 0, color = "black", linewidth = .3) +
  theme_bw(base_size = 16) +
  theme(axis.text.x = element_text(angle = 40, hjust = 1),
        # legend.position = "none",
        panel.spacing = unit(1, "pt"),
        panel.grid = element_line(size = .3)) +
  guides(color = guide_legend(override.aes = list(stroke = 1, shape = 19)))  +
  
  geom_point(data = d.obs.simu %>%
               group_by(tracer, metabolite, label) %>%
               summarise(sim = mean(sim),
                         obs = mean(obs)),
             aes(fill = tracer), color = "black", fill = "black", size = .5,
             shape = 21)

p02 <- p + coord_equal(xlim = c(0, .2), ylim = c(0, .2));   p02
ggsave(filename = "./plots/sim vs obs LOW range.pdf", width = 5, height = 5)

p81 <- p + coord_equal(xlim = c(0.9, 1), ylim = c(0.9, 1)); p81
ggsave(filename = "./plots/sim vs obs HIGH range.pdf", width = 5, height = 5)



set.seed(123)
myColors <- colorRampPalette(c(RColorBrewer::brewer.pal(8, "Dark2") %>% sample(), "black")) (d.obs.simu$metabolite %>% n_distinct())

p <- d.obs.simu %>% 
  mutate(metabolite = str_replace(metabolite, "\\.Blood", "\\.Bld")) %>% 
  ggplot(aes(x = obs, y = sim, color = metabolite)) +
  scale_x_continuous(transform = "sqrt") +
  scale_y_continuous(transform = "sqrt") +
  facet_wrap(~tracer, nrow = 2) +
  geom_abline(slope = 1, intercept = 0, color = "black", linewidth = .3) +
  theme_bw(base_size = 13) +
  theme(axis.text.x = element_text(angle = 40, hjust = 1),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        legend.position = "right",
        panel.spacing = unit(5, "pt"),
        panel.grid = element_line(size = .3)) +
  geom_text(data = d.obs.simu %>%
              group_by(tracer, metabolite, label) %>%
              summarise(sim = mean(sim),
                        obs = mean(obs)),
            aes(color = metabolite, label = label %>% str_remove("#")),
            size = 3, key_glyph = draw_key_point) +
  guides(color = guide_legend(override.aes = list(size = 5, fontface = "bold"))) +
  labs(color = NULL) +
  scale_color_manual(values = myColors)

p02.facet <- p + coord_equal(xlim = c( 0, .2), ylim = c(0, .2));   p02.facet
ggsave(filename = "./plots/sim vs obs LOW range faceted.pdf", width = 10, height = 5)

p81.facet <- p + coord_equal(xlim = c(.8, 1),  ylim = c(.8, 1));   p81.facet
ggsave(filename = "./plots/sim vs obs HIGH range faceted.pdf", width = 10, height = 5)
