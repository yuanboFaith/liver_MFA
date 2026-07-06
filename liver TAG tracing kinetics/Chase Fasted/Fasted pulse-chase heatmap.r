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
  
  # # keep only mice that has both serum and liver data
  # group_by(mouse) %>%
  # filter(any(tissue == "serum") & any(tissue == "LIVER")) %>%
  # ungroup() %>% 
  
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
d.labeling.Fast.Chase$mouse %>% unique()


# Clean up
d.labeling.Fast.Chase.clean <- d.labeling.Fast.Chase # %>% 

# # keep only mice that has both serum and liver data
#   group_by(mouse) %>%
#   filter(any(tissue == "serum") & any(tissue == "liver")) %>%
#   ungroup() 

d.labeling.Fast.Chase.clean %>% select(Compound, tissue) %>% table()
d.labeling.Fast.Chase.clean







# $$-$$-$$-$$-$$-$$-$$-$$-$$-$$-$$-$$-$$-$$-$$-$$-$$-$$-$$-$$-$$-$$-$$-$$-$$-$$-$$-$$-$$-$$-$$-$$-$$-$$-$$-$$-$$-$$-$$-$$-$$-$$-$$-$$-$$-$$-$$-$$-
x.Fasted.pulse <- d.labeling.Fast.Pulse.clean %>% select(-c(`time (h)`, batch)) # mouse id as big letters

x.Fasted.chase <- d.labeling.Fast.Chase.clean %>% 
  select(-c(timepoint, `time (h)`, sample)) %>% 
  relocate(Compound, C_Label, tissue, mouse, labeling) # %>% 
  # convert mouse id to small letters for chase
  # mutate(mouse = letters[as.integer(mouse)]) 

# ensure the mouse id of the two datasets are unique ()
x.Fasted.pulse$mouse %>% unique() %>% sort()
x.Fasted.chase$mouse %>% unique() %>% sort()

# x.Fasted.chase <- x.Fasted.chase %>% 
#   
#   # there is a slight M+1 labeling not removed completely from natural abundance correction
#   # due to the large n available, here this mouse is simply removed to make code faster
#   filter(mouse != "e") %>%
# 
#   # these two mice present data structure issue (error in MFA: one node produced an error: 'names' attribute [68] must be the same length as the vector [39])
#   # due to the large n we have, here these two mice data are thrown away
#   filter(! mouse %in% c("i", "h")) %>%
#   
#   mutate(mouse = str_replace(mouse, "a", "O")) %>% 
#   mutate(mouse = str_replace(mouse, "b", "P")) %>% 
#   mutate(mouse = str_replace(mouse, "c", "Q")) %>% 
#   mutate(mouse = str_replace(mouse, "d", "R")) %>% 
#   mutate(mouse = str_replace(mouse, "f", "S")) %>% 
#   mutate(mouse = str_replace(mouse, "g", "T")) 


d.Fast.kinetics <- bind_rows(x.Fasted.pulse, x.Fasted.chase)
# d.Fast.kinetics <- bind_rows(x.Fasted.pulse)
x.Fasted.pulse $mouse %>% unique()
x.Fasted.chase $mouse %>% unique()
d.Fast.kinetics$mouse %>% unique()


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

filter(!is.na(m.id))

# # each mouse need to have the Blood Lino and Liver TAGLino
# group_by(m.id) %>% 
# filter(any(Compound == "Lino") & any(Compound == "TAGLino")) %>% 
# ungroup()


d.dataForMFA.TAGkinetics.Fasted$m.id %>% unique()
d.dataForMFA.TAGkinetics.Fasted$Compound %>% unique()




# a quick visual check
d.dataForMFA.TAGkinetics.Fasted %>% 
  filter(Compound %in% c("Lino", "TAGLino" )) %>% 
  filter(C_Label != 0) %>% 
  ggplot(aes(x = Compound, y = labeling, fill = as.character(C_Label))) +
  geom_col() + facet_wrap(~m.id) + theme(legend.position = "right")




# Check mice level order
d.dataForMFA.TAGkinetics.Fasted$m.id %>% unique()


mouse_pulse_time <- tribble(
  ~mouse, ~time_h,
  "M", 2.53,
  "H", 2.57,
  "L", 2.57,
  "N", 2.65,
  "I", 2.68,
  "G", 2.73,
  "E", 5.05,
  "F", 5.05,
  "J", 7.03,
  "K", 7.65,
  "C", 9.50,
  "A", 9.92,
  "B", 11.67,
  "D", 11.75
)


mouse_chase_time <- tribble(
  ~mouse, ~time_h,
  4, 0.00,
  7, 0.00,
  6, 1.07,
  8, 1.13,
  9, 1.22,
  5, 1.57,
  1, 2.00,
  2, 2.08,
  3, 2.13
)

d.mouse.order <- rbind(mouse_pulse_time, mouse_chase_time)
ordered.mouse <- d.mouse.order$mouse; ordered.mouse


# plot heatmap vertically

x.iso  <- c("Glc", "Lac", "Ala", "Gln", "HB", "Mal", "Suc"); x.iso  # show isotopologue labeling pattersn

d.heatmap <- d.dataForMFA.TAGkinetics.Fasted %>% 
  filter(Compound %in% x.iso) %>% 
  mutate(Compound = factor(Compound, levels = x.iso, ordered = T)) %>% 
  arrange(Compound) %>% 
  filter( (tissue == "Lv" & Compound %in% c("Mal", "Suc")) | 
            (tissue == "Blood" & Compound %in% c("Glc", "Lac", "Ala", "Gln", "HB")) ) %>% 
  
  # reorder mice in order
  mutate(m.id = factor(m.id, levels = ordered.mouse, ordered = T))

d.heatmap      

d.heatmap$m.id %>% unique()

myColors <- colorRampPalette(
  c("grey90", "#d1e6f7", "#a6d3f1", "#57b7e7", "#1bab90", "#3aab70", 
    "#fdd53e", "orange", "firebrick2", "firebrick4", "black"),
  bias = 4)(100) 



library(patchwork)


func.heatmap.labeling <- function(
    bloodOrTissue = "Blood"){
  
  max.labeling <- .1 # set an upper max labeling beyond which the same deep saturated color is used
  
  # test
  # bloodOrTissue = "Blood"
  
  if (bloodOrTissue == "Blood") {
    x <- d.heatmap %>% 
      filter(tissue == "Blood")
  } else {
    x <- d.heatmap %>% 
      filter(tissue != "Blood")
  }
  
  y <- x %>% 
    # keep only the full carbon number in FAs
    filter(! (Compound %in% "palmitate" & C_Label != 16)) %>% 
    filter(! (Compound %in% "oleate" & C_Label != 18)) %>% 
    filter(! (Compound %in% "linoleate" & C_Label != 18)) %>% 
    
    # compound isotopologues
    mutate(Compound_C_Label = paste(Compound, C_Label), .after = 2) %>% 
    mutate(Compound_C_Label = factor(Compound_C_Label, levels = .$Compound_C_Label %>% unique() %>% rev() )) %>% 
    filter(C_Label != 0) %>% 
    
    # use the same deep saturated color for very big labeling
    mutate(labeling = ifelse(labeling > max.labeling, max.labeling, labeling))
  
  p.main <- y %>% 
    ggplot(aes(x = m.id, y = Compound_C_Label, fill = labeling)) +
    geom_tile(color = "white") +
    # geom_text(aes(label = round(labeling, 3)), size = .5) +
    # geom_raster() +
    # facet_grid(tissue ~ Infusate, scales = "free_x", space = "free") +
    # scale_fill_distiller(palette = "Spectral") +
    scale_fill_gradientn(colours = myColors, 
                         breaks = seq(from = 0, to = .5, by = .02), 
                         limits = c(0, max.labeling), # ensure the same color scale for the specified range of data
                         values = seq(0, 1, length.out = length(myColors))) +
    theme_classic(base_size = 19) +
    theme(
      strip.clip = "off",
      # axis.text.x = element_text(angle = 50, hjust = 1, size = 6),
      axis.text = element_blank(),
      axis.title = element_blank(),
      axis.ticks.length.x = unit(-2, "pt"),
      axis.ticks = element_blank(),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", size = 10),
      strip.text.x.top = element_text(angle = 90, hjust = 0),
      panel.border = element_rect(color = "black", fill = NA, linewidth = .5),
      panel.spacing = unit(0, "pt"),
      legend.text = element_text(angle = 60, hjust = 1),
      legend.position = "bottom"
    ) +
    scale_x_discrete(expand = expansion(add = 0)) +
    scale_y_discrete(expand = expansion(add = 0)) +
    guides(fill = guide_colorbar(barwidth = unit(200, "pt"),
                                 barheight = unit(5, "pt"), title = NULL)) 
  
  
  # add color bar denoting compounds
  p.sideBar <- y %>% select(Compound_C_Label, Compound, tissue) %>% 
    
    mutate(whichCompound = ifelse(Compound == "a-ketoglutarate", "aKG", as.character(Compound))) %>% 
    # if the isotope label is 1, or 16, 18 for fatty acids, then put the compund label there
    mutate(whichCompound = ifelse(str_extract(Compound_C_Label, "\\d{1,2}$") %in% c("1", "16", "18"), 
                                  paste(as.character(whichCompound), "   "), NA)) %>% 
    
    ggplot(aes(x = 0, y = Compound_C_Label, fill = Compound, color = Compound)) +
    # geom_tile(color = "white", size = .5) +
    
    # label isotopologue index M+1, M+2....
    geom_text(aes(label = Compound_C_Label %>% str_extract("\\d{1,2}$"), x = -1),
              size = 6, hjust = 1) +
    
    # label compound name at M+1 position
    geom_text(aes(label = whichCompound, x = -1),
              size = 6, hjust = 1) +
    
    facet_grid(tissue~.) +
    scale_x_continuous(expand = expansion(mult = 0)) +
    coord_cartesian(clip = "off") +
    theme_minimal() +
    theme(
      plot.margin = margin(l = 40, unit = "pt"),
      strip.background = element_blank(),
      strip.text = element_blank(),
      legend.position = "none", 
      panel.spacing = unit(0, "pt"), # needs to be identical as the main plot
      axis.title = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid = element_blank(),
      axis.text = element_blank()
    ) +
    scale_fill_brewer(palette =  ifelse(bloodOrTissue == "Blood", "Paired", "Set2")) +
    scale_color_brewer(palette =  ifelse(bloodOrTissue == "Blood", "Paired", "Set2")) 
  
  p.sideBar + 
    p.main + 
    plot_layout(widths = c(.3, 8))
  
}


# fasted
x1 <- func.heatmap.labeling(bloodOrTissue = "Blood")   #; x1
x2 <- func.heatmap.labeling(bloodOrTissue = "Tissues") #; x2
cowplot::plot_grid(x1 + theme(plot.margin = margin(b = 0), legend.position = "none"), 
                   
                   ggplot() + theme_void(),
                   
                   x2 + theme(strip.text.x.top = element_blank(),
                              plot.margin = margin(t = 0)), 
                   ncol = 1, 
                   align = "v",
                   rel_heights = c(4.8, .12,  3.2))


ggsave("fasted TAG tracing kinetics heatmap.pdf", height = 9.5, width = 3.5)

