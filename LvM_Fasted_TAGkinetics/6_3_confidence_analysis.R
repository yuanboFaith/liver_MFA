library(readxl)
library(tidyverse)
rm(list = ls())

rstudioapi::getActiveDocumentContext()$path %>% dirname() %>% setwd()



myClassic <- theme_classic(base_size = 14) +
  theme(axis.title.y = element_text(margin = margin(r = 5, unit = "pt")),
        axis.title.x = element_text(margin = margin(t = 5, unit = "pt")),
        axis.title = element_text(face = "bold"),
        strip.background = element_blank(),
        strip.text = element_text(size = 15, face = "bold"))




# Read all the files of CI for each flux --------------------------------------------
# print existing files of CI
list.files("CI", full.names = TRUE) 

# Import simulated result
d.CI <- tibble(.rows = 1)

for (file in list.files("CI")) { # files are ordered by 1, 10, 11, 12...19, 2, 20, 21...
  print(file)
  # check
  # file <- list.files("CI")[1]
  # Load the file into the temporary environment
  load(file = paste0("CI/", file))
  
  # Retrieve objects from the environment and sequentially compile them together
  d.CI <- rbind(d.CI, d.flux.cost)
  
}


# optimal_cost
# F.critical.95 <- qf(p = .95,  df1 = 1, df2 = n.minus.p, lower.tail = T)
# cost.threshold.95 <- F.critical.95 / n.minus.p * optimal_cost + optimal_cost; cost.threshold.95




d.CI <- d.CI[-1, ] %>% # remove the first empty line
  # filter(cost >= optimal_cost) %>% 
  mutate(flux.index = factor(flux.index)) %>% 
  
  mutate(CI95 = cost.threshold.95) %>% # upper boundary limit indicating 95% CI
  # mark lowest point for each flux
  group_by(flux.index) %>% 
  mutate(lowestPoint = ifelse(cost == min(cost), T, F))

d.CI_side <- d.CI %>% left_join(
  # find optimal flux value
  d.CI %>% filter(lowestPoint == T) %>% rename(flux.optimal = flux.value) %>% select(flux.index, flux.optimal) 
) %>% 
  # find left and right side of the C.I.
  mutate(side = ifelse(flux.value > flux.optimal, "R", "L")) %>% 
  mutate(side = ifelse(lowestPoint == T, "C", side)) %>% # right, left, or center (optimal point)
  
  # mark the index of steps
  group_by(flux.index, side) %>% 
  mutate(binary.index = 1:length(side))




# Get the optimal flux acquired from the repeated run
d.optimal.repeatRun <- optimal_allFlux %>% as_tibble() %>% 
  mutate(flux.index = 1:nrow(.) %>% factor()) %>%  
  rename(flux.optimal.repeatRun = V1)

# add the Net differences from repeat runs
# d.CI_side %>% filter(flux.index %>% str_detect("-")) %>% # select net diffference reactions
#   separate(flux.index, into = c("a", "b"), sep = "-", remove = F)

flux.index <- d.CI_side$flux.index %>% unique()
flux.index.net.diff <- flux.index[flux.index %>% str_detect("-")] %>% as.character()

for (a in flux.index.net.diff) {
  # test a = flux.index.net.diff[1]
  print(a)
  index.r1_r2 <- str_split_1(a, pattern = "-")  # reaction index of the two reactions
  
  f.r1 <- filter(d.optimal.repeatRun, flux.index == index.r1_r2[1]) $ flux.optimal.repeatRun # flux of the first reaction
  f.r2 <- filter(d.optimal.repeatRun, flux.index == index.r1_r2[2]) $ flux.optimal.repeatRun # flux of the second reaction
  d.new.f.i <- tibble("flux.optimal.repeatRun" = f.r1 - f.r2, "flux.index" = a) # a data frame of new net flux to add to the optimal repeat run
  d.optimal.repeatRun <- bind_rows(d.optimal.repeatRun, d.new.f.i)
}





# add summation of fluxes
# d.CI_side %>% filter(flux.index %>% str_detect("-")) %>% # select net diffference reactions
#   separate(flux.index, into = c("a", "b"), sep = "-", remove = F)

flux.index.sum <- flux.index[flux.index %>% str_detect("\\+")] %>% as.character()

for (a in flux.index.sum) {
  # test a = flux.index.sum[1]
  print(a)
  index.rs <- str_split_1(a, pattern = "-") %>% str_split_1(pattern = "\\+")# reaction index of the fluxes summed up
  
  # extract the flux for each reaction, and then sum them up
  f.sum <- d.optimal.repeatRun %>%
    filter(flux.index %in% index.rs) %>%
    summarise(total = sum(flux.optimal.repeatRun)) %>% pull(total) 
  
  d.new.f.i <- tibble("flux.optimal.repeatRun" = f.sum, "flux.index" = a) # a data frame of new fluxes summed up
  d.optimal.repeatRun <- bind_rows(d.optimal.repeatRun, d.new.f.i)
}

d.CI_side <- d.CI_side %>% left_join(d.optimal.repeatRun, by = "flux.index")







# select fluxes of interest
d.CI.selected <- d.CI_side # %>%
# filter(flux.index %in% c(12, 25))


# plot -----
d.CI.selected %>%  
  ggplot(aes(x = flux.value, y = cost)) + 
  
  # highlight the C.I. bowl shape
  geom_ribbon(data = d.CI.selected %>% filter(cost <= cost.threshold.95),
              aes(ymax = CI95, ymin = cost, xmin = flux.value, xmax = flux.value),
              fill = "snow2", alpha = .7) +
  
  # CI boundary 95
  geom_hline(yintercept = c(cost.threshold.95), color = "snow4",linewidth = .5, linetype = "dashed") +
  # mark fluxes determined by binary search
  geom_vline(data = d.CI.selected %>% filter(side != "C"),
             aes(xintercept = flux.value, color = side), linewidth = .2, alpha = .7) +
  
  # CI curve
  geom_line(aes(group = flux.index), linewidth = .6) +
  
  geom_point(size = 1) +
  
  # mark lowest point acquired from the CI determination procedure
  geom_segment(data = filter(d.CI.selected, lowestPoint == T), 
               aes(x=flux.value, xend = flux.value, yend = optimal_cost, y = Inf),
               linetype = "dashed", color = "green3", linewidth = .5) +
  
  # mark lowest point acquired from repeated run (this should be consistent with the lowest point from CI analysis)
  geom_vline(aes(xintercept = flux.optimal.repeatRun), color = "purple", linetype = "dotted") +
  
  # mark the flux number
  geom_text(data = filter(d.CI.selected, lowestPoint == T), 
            # if flux is zero, then mark the flux number shifted 5 units to the right
            # otherwise mark using the flux value as is
            aes(x = ifelse( round(flux.value)==0, flux.value + 5, flux.value), 
                label = flux.index, y = (optimal_cost + cost.threshold.95)/2), 
            size = 3, fontface = "bold", color = "red4") +
  
  # mark the binary steps
  geom_text(data = filter(d.CI.selected, lowestPoint != T), 
            aes(y = cost.threshold.95 * 1.01, x = flux.value, label = binary.index, color = side),
            size = 2.5, alpha = .7) +
  
  theme_classic(base_size = 14) +
  theme(strip.background = element_blank(),
        panel.background = element_rect(colour = "black"),
        strip.text = element_blank(),
        # strip.text = element_text(face = "bold", size = 15, vjust = -3), 
        strip.placement = "outside",
        panel.spacing.y =  unit(10, "pt"),
        panel.spacing.x = unit(0, "pt"),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, margin = margin(t = 1)),
        legend.position = "none") +
  facet_wrap(~flux.index, scales = "free_x", nrow = 6) +
  scale_x_continuous(limits = c(0, NA), expand = expansion(mult = c(0, .2)), position = "bottom") +  
  scale_y_continuous(expand = c(0, 0),
                     labels = function(x) x * 1000,
                     breaks = seq(optimal_cost, cost.threshold.95, length.out = 5) %>% round(5),
                     name = expression(cost~(x~10^-3))) +
  coord_cartesian(ylim = c(optimal_cost, cost.threshold.95 + diff(c(optimal_cost, cost.threshold.95)) * 0.4))


ggsave(filename = "./plots/CI profile_fasted.pdf", height = 10, width = 15) 




# get C.I.
d.CI.bounds <- d.CI.selected %>% 
  group_by(flux.index) %>% 
  filter(cost < CI95) %>% 
  filter(side != "C") %>% 
  # find the max value on each side
  group_by(side, .add = T) %>% 
  filter(cost == max(cost)) %>%
  arrange(flux.index) %>% 
  ungroup()




x1 <- d.CI.bounds %>% select(flux.index, flux.value, side) %>%  # two side bounds
  pivot_wider(names_from = side, values_from = flux.value) %>% 
  mutate(L = replace_na(L, 0))
x2 <- d.CI.bounds %>% select(flux.index, flux.optimal) %>% distinct() # optimal

d.CI.bounds <- left_join(x1, x2, by = "flux.index") 





# combine with reaction names 
reaction_data.clean <- reaction_data %>% 
  select(enzyme, reactions) %>% # remove irrelevant columns
  rename(flux.index = enzyme) %>% mutate(flux.index = flux.index %>% as.character())

# add names for Net Difference Reactions
for (a in flux.index.net.diff) {
  # test a = flux.index.net.diff[1]
  print(a)
  index.r1_r2 <- str_split_1(a, pattern = "-")  # reaction index of the two reactions
  
  name.r1 <- filter(reaction_data.clean, flux.index == index.r1_r2[1]) $ reactions # name of the first reaction
  name.r2 <- filter(reaction_data.clean, flux.index == index.r1_r2[2]) $ reactions # name of the second reaction
  d.new.name.i <- tibble("flux.index" = a, "reactions" = paste0("(", name.r1, ") - (", name.r2, ")") ) # a data frame of new net flux to add to the optimal repeat run
  
  reaction_data.clean <- bind_rows(reaction_data.clean, d.new.name.i) # add this new reaction name to the reaction name database
}


# add names for Sum of Fluxes
for (a in flux.index.sum) {
  # test a = flux.index.sum[1]
  print(a)
  index.rs <- a %>% str_split_1(pattern = "\\+")
  
  # collect names of reactions to be summed up
  name.r <- reaction_data.clean %>%
    filter(flux.index %in% index.rs) %>%
    mutate(reactions = paste0("(",reactions, ")")) %>% 
    pull(reactions) %>% 
    paste0(collapse = "+") # %>% 
  # str_wrap(width = 10)
  
  d.new.name.i <- tibble("flux.index" = a, "reactions" = name.r) # a data frame of new fluxes to to add to the optimal repeat run
  
  reaction_data.clean <- bind_rows(reaction_data.clean, d.new.name.i) # add this new reaction name to the reaction name database
}


d.CI.bounds.reactions <- d.CI.bounds %>% left_join(reaction_data.clean) 
d.CI.bounds.reactions 





# More clean up 
d.CI.bounds.reactions <- d.CI.bounds.reactions %>% 
  # calculate the mid value of CI range
  mutate(flux.center = (R+L)/2, .after = flux.optimal) %>% 
  # add the index number before the reaction name
  mutate(index.reactions = str_c(flux.index, ". ", reactions)) %>% 
  
  # arrange by the flux index
  # for net reaction, arrange by the last reaction index
  mutate(index.plot.order = str_extract(flux.index, pattern = "\\d{1,3}$") %>% as.integer()) %>% 
  mutate(index.reactions = fct_reorder(index.reactions, index.plot.order)) %>% 
  mutate(reactions       = fct_reorder(reactions,       index.plot.order)) 


# liver reactions
d.CI.bounds.reactions.lv <- d.CI.bounds.reactions %>% 
  filter(str_detect(reactions, "\\.Lv"))




# # plot confidence interval - version 1
# d.CI.bounds.reactions.lv %>% 
#   ggplot(aes(x = reactions, xend = reactions, y = L, yend = R)) +
#   geom_segment() +
#   # best fit
#   geom_point(aes(y = flux.optimal), size = 3, shape = 21, fill = "tomato") +
#   ggrepel::geom_text_repel(aes(label = round(flux.optimal), y = flux.optimal), color = "tomato") +
#   
#   # # center of C.I.
#   # geom_point(aes(y = flux.center),  size = 3, shape = 23, fill = "tomato") +
#   # ggrepel::geom_text_repel(aes(label = round(flux.center), y = flux.center), color = "tomato") +
#   
#   # Left of C.I.
#   ggrepel::geom_text_repel(aes(label = round(L), y = L), color = "turquoise4", box.padding = unit(0, "pt")) +
#   # Right of C.I.
#   ggrepel::geom_text_repel(aes(label = round(R), y = R), color = "turquoise4", box.padding = unit(0, "pt")) +
#   
#   # # name of reactions index
#   geom_text(aes(label = flux.index), y = max(d.CI.bounds.reactions.lv$R+20), 
#             color = "steelblue2", fontface = "bold") +
#   
#   scale_y_continuous(expand = expansion(mult = c(0, .1)), breaks = seq(0, 300, 50)) +
#   myClassic +
#   theme(axis.text.x = element_text(angle = 50, hjust = 1),
#         axis.ticks.x = element_blank()) +
#   coord_cartesian(clip = "off") +
#   labs(y = "confidence interval (nmol products / g BW / min)", x = NULL) +
#   scale_x_discrete(labels = scales::label_wrap(width = 15))  
# 
# ggsave("./plots/Confidence interval_horizontal.pdf", width = 25, height = 10)




# plot confidence interval - version 2
func.plot.CI <- function(L.limit = 0, R.limit = 500, shiftName = 40) {
  
  # test L.limit = 0; R.limit = 500
  
  dd <- d.CI.bounds.reactions.lv %>% 
    filter(R < R.limit) %>% 
    filter(L > L.limit) 
  
  message(nrow(dd), " reactions found")
  
  dd %>%   
    ggplot(aes(x = reactions, xend = reactions, y = L, yend = R)) +
    geom_segment() +
    # best fit
    geom_point(aes(y = flux.optimal), size = 3, shape = 21, fill = "tomato") +
    ggrepel::geom_text_repel(aes(label = round(flux.optimal), y = flux.optimal), color = "tomato", nudge_x = .1) +
    
    # center of C.I.
    geom_point(aes(y = flux.center),  size = 3, shape = 3, color = "black") +
    # ggrepel::geom_text_repel(aes(label = round(flux.center), y = flux.center), color = "tomato", nudge_x = .1) +
    
    
    # Left of C.I.
    geom_text(aes(label = round(L), y = L - 10), color = "turquoise4", hjust = 1, size = 3) + # nudge_x = .1,
    # Right of C.I.
    geom_text(aes(label = round(R), y = R + 10), color = "turquoise4", hjust = 0, size = 3) +
    
    
    # name of reactions index
    geom_text(aes(label = flux.index, y = -30), color = "black", fontface = "bold", hjust = 1, angle = 60) +
    # name of reactions
    geom_text(
      aes(label = reactions, y = R + shiftName), 
      hjust = 0, 
      # nudge_x = -.05, # slightly shift downward
      color = "grey30", # box.padding = unit(0, "pt"),force = 0, 
    ) +
    
    scale_y_continuous(expand = expansion(mult = c(0, .3)), breaks = seq(0, 300, 50)) +
    myClassic +
    theme(axis.text.y = element_blank(),
          plot.margin = margin(l = 100, t = 30, b = 30, r = 80),
          axis.line.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          axis.line = element_blank(),
          axis.title.x = element_blank()
    ) +
    labs(y = "confidence interval (nmol products / g BW / min)", x = NULL) +
    coord_flip(expand = expansion(add = 1), clip = "off") +
    scale_x_discrete(limits = rev, expand = expansion(add = .2)) %>% 
  return()
  
  
}

p.small  <- func.plot.CI(L.limit = 0,    R.limit = 500, shiftName = 40)
p.small
# p.big   <- func.plot.CI(L.limit = 500.1, R.limit = 4000, shiftName = 200)
# cowplot::plot_grid(p.small, p.big, rel_heights = c(29, 3), ncol = 1, align = "v")
ggsave("./plots/Confidence interval_vertical.pdf", width = 10, height = 12)


# note the tissue name and extract clean reaction name
# Note the tissue name
d.CI.bounds.reactions.reduced.LvM <- d.CI.bounds.reactions.lv %>% 
  # remove blood, and extract tissue name
  mutate(tissue = str_remove(reactions, "\\.Blood") %>% str_extract("\\.[a-zA-Z]{1,6}")) %>% 
  # remove mitochondria suffix
  mutate(tissue = str_remove(tissue, "m$") %>% str_remove(".")) %>% 
  # remove tissue name from reaction
  mutate(reaction.pure = str_remove_all(reactions, "\\.[a-zA-Z]{1,6}"))

d.CI.bounds.reactions.reduced.LvM

save(d.CI.bounds.reactions.reduced.LvM, file = "6_CI_fasted.RData")
