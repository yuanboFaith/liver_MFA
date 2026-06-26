rm(list = ls())

library(tidyverse)
library(readxl)
library(patchwork)
# default path to the current folder
rstudioapi::getActiveDocumentContext()$path %>% dirname() %>% setwd(); getwd()

# List all files and directories in the specified directory
files_and_dirs <- list.files(path = ".", full.names = TRUE)


# load .RData files (the cleaned data from each tracer experiment) into global environment
for (folder in files_and_dirs) {
  # Check if the path is a directory
  if (dir.exists(folder)) {
    # List all .RData files in the directory
    rdata_files <- list.files(path = folder, pattern = "\\.RData$", full.names = TRUE)
    
    # Load the .RData file
    for (rdata_file in rdata_files) {
      load(rdata_file, envir = .GlobalEnv)
      message("Loaded: ", rdata_file)
    }
  }
}

# Retrieve all datasets with names starting with "d.13C_"
all_objects <- ls(envir = .GlobalEnv) # Get all object names in the global environment
matching_objects <- all_objects[str_detect(all_objects, "^d.13C_")] # Filter for object names that start with the specified prefix

# Extract these objects from the global environment, and combine into a tibble
d.13C <- mget(matching_objects, envir = .GlobalEnv) %>% bind_rows()


# update 3HB name
d.13C <- d.13C %>% mutate(Compound = str_replace(Compound, "3-hydroxybutyrate", "3-HB"))

d.13C$State %>% unique()
d.13C$Infusate %>% unique()
d.13C$tissue %>% unique()
d.13C$Compound %>% unique()






# compare labeling in TA and Quadriceps
exp.containing_TA <- filter(d.13C, tissue == "TA")$mouse.when.who  %>% unique() # experiments containing the tissue of TA

# select TA and Q, and representative metabolites
func.plot <- function(myState = "fasted", myInfusate = "glucose") {
  d.TA_Q <- d.13C %>% 
    filter(mouse.when.who %in% exp.containing_TA ) %>% 
    filter(Infusate %>% str_detect(myInfusate)) %>% 
    filter(tissue %in% c("TA", "Q")) %>% 
    filter(Compound %in% c("glucose", "lactate", "alanine", "malate", "succinate")) %>% 
    filter(State == myState) %>% 
    filter(C_Label != 0 )
  
  message("infused tracer: ",       d.TA_Q$Infusate %>% unique()             )
  message("tracer concentration: ", d.TA_Q$`Concentration (mM)` %>% unique() )
  message("tracer infuion rate: ",  d.TA_Q$`Rate (µl/min/g)` %>% unique()    )
  message("mouse, when, who: ",     d.TA_Q$mouse.when.who %>% unique()       )
  
  d.TA_Q %>% 
    mutate(C_Label = factor(C_Label)) %>% 
    ggplot(aes(x = tissue, y = labeling, color = C_Label)) +
    geom_point() +
    geom_line(aes(group = C_Label)) +
    facet_grid(Compound ~ mouse.when.who, scales = "free") +
    theme_bw(base_size = 15) +
    ggtitle(paste0("[U-13C]", myInfusate, ", ", myState))
}

# experiments containing TA and Q
d.13C %>% filter(mouse.when.who %in% exp.containing_TA ) %>% select(Infusate, State) %>% distinct()



func.plot(myState = "fasted", myInfusate = "glucose")
func.plot(myState = "refed",  myInfusate = "glucose")
func.plot(myState = "fasted", myInfusate = "glutamine")
func.plot(myState = "fasted", myInfusate = "lactate")

# When there is TA, there is Q collected from the same mouse. We'll use the Q as the more representative muscle

# remove extra tissues
d.13C.selectedTissues <- d.13C %>% 
  # in TT-alanine data, there is unclear tissue "in", along with gWat, iWat and SI. so "in" must not be small intestine. Is it intestine in general?
  filter(! tissue %in% c("TA", "t", "gast", "sol", "lymphn", "skin", "thymus", "LI", "in"))


(d.13C.selectedTissues %>% 
    filter(mouse.when.who == "M1_2016-06-30_SH"))$tissue %>% unique()


# unify the tissue and compound name
d.13C.selectedTissues$tissue %>% unique()

d.tissueNames <- readxl::read_excel("suppl_data.xlsx", sheet = "tissueNames")
d.tissueNames
# create a named vector of tissue names, with element names being the empirical names, and elements being the finalized abbreviated names
tissueNames.abbr <- d.tissueNames$tiss
names(tissueNames.abbr) <- d.tissueNames$tissue

d.tissueRenamed <- d.13C.selectedTissues %>% 
  mutate(tissue = tissueNames.abbr[tissue]) %>% 
  
  # some further non-sophisticated cleanup
  select(-c(When, Who, mouse)) %>%  # here we'll use 'mouse.when.who' as the unique mouse ID
  mutate(Infusate = str_remove(Infusate, "\\[U-13C\\]")) %>%  # remove "[U-13C]" prefix from the tracer name
  # rename fatty acids
  mutate(Compound = str_replace(Compound, "palmitic acid", "palmitate"),
         Compound = str_replace(Compound, "oleic acid",    "oleate"   ),
         Compound = str_replace(Compound, "linoleic acid", "linoleate")) %>% 
  # calculate infusion rate
  mutate(infuse.nmol.min.g = `Concentration (mM)` * `Rate (µl/min/g)`)


d.tissueRenamed$tissue          %>% unique()
d.tissueRenamed$mouse.when.who  %>% unique()
d.tissueRenamed$State           %>% unique()
d.tissueRenamed$Infusate        %>% unique() 
d.tissueRenamed$Compound        %>% unique()



# keep only relevant metabolites in the dataset
ordered.tracers <- c("glucose", "lactate", "alanine", "glutamine", "palmitate", "oleate", "linoleate", "glycerol", "citrate", "3-HB")
ordered.TCAs    <- c("malate",  "succinate", "citrate", "a-ketoglutarate", "glutamate") # , "citrate", "a-ketoglutarate")
ordered.tissues <- c("Blood", "Lv", "M", "H", "Ln", "K", "P", "S", "I", "Br", "Bat", "gW", "iW")

d.compoundSelected <- d.tissueRenamed %>% 
  # selected tracing experiment
  filter(Infusate %in% ordered.tracers) %>% 
  
  # keep only tracer metabolites in the blood
  filter(! (tissue == "Blood" &  !Compound %in% c(ordered.tracers))) %>% 
  # in tissues, keep only malate, succinate, and citrate and aKG (if any)
  filter(! (tissue != "Blood" &  !Compound %in% ordered.TCAs)) # may also include in the model: "citrate", "a-ketoglutarate"

d.compoundSelected




# arrange data in order of tracer, compound, tissue, and mouse ID
d.compoundSelected_ordered <- 
  d.compoundSelected %>% 
  mutate(Infusate = factor(Infusate, levels = ordered.tracers),
         Compound = factor(Compound, levels = c(ordered.tracers, ordered.TCAs) %>% unique()),
         tissue   = factor(tissue, levels = ordered.tissues)) %>% 
  # put mouse ID in such order
  arrange(State, Infusate, mouse.when.who, tissue, Compound) %>% 
  mutate(mouse.when.who = factor(mouse.when.who, levels = .$mouse.when.who %>% unique()))



# check data structure
# 1. Check duplication
x <- d.compoundSelected_ordered %>% 
  select(State, Infusate, mouse.when.who, tissue, Compound, C_Label) %>% duplicated()
sum(x)


# 2. plot: overall structure, per organ
d.compoundSelected_ordered %>% 
  select(State, Infusate, mouse.when.who, tissue, Compound) %>% distinct() %>% 
  ggplot(aes(x = tissue, y = Compound, color = Infusate)) +
  geom_point(position = position_jitter(.2, .2, 123), size = 3, shape = 21) +
  theme_bw(base_size = 15) +
  facet_wrap(~State)


# 3. plot: overall structure, per infusate
d.compoundSelected_ordered %>% 
  # filter(Compound %in% c("succinate", "malate")) %>% 
  filter(Compound %in% c("succinate", "malate", "citrate", "a-ketoglutarate", "glutamate")) %>% 
  select(State, Infusate, mouse.when.who, tissue, Compound) %>% distinct() %>% 
  ggplot(aes(x = Infusate, y = tissue, color = Compound)) +
  geom_point(position = position_jitter(.2, .2, 123), size = 3, shape = 21) +
  theme_bw(base_size = 15) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_wrap(~State) 



# 4. plot: mouse-wise
d.compoundSelected_ordered %>% 
  # filter(State == "refed") %>% 
  filter(State == "fasted") %>% 
  select(Infusate, mouse.when.who, tissue, Compound) %>% distinct() %>% 
  ggplot(aes(x = mouse.when.who, y = Compound, color = Infusate)) +
  geom_point(size = 2) +
  facet_wrap(~tissue, nrow = 1, scales = "free_x") +
  coord_flip() +
  theme_bw(base_size = 15)  +
  theme(axis.text.x = element_text(angle = 50, hjust = 1),
        axis.text.y = element_text(size = 10))


# 5. plot: infusion parameters
# 5.1
d.compoundSelected_ordered %>% 
  select(State, Infusate, `Concentration (mM)`, `Rate (µl/min/g)`, mouse.when.who) %>% distinct() %>% 
  mutate(`Rate (µl/min/g)`    = factor(`Rate (µl/min/g)`),
         `Concentration (mM)` = factor(`Concentration (mM)`)) %>% 
  ggplot(aes(x = `Concentration (mM)`, y = `Rate (µl/min/g)`, color = mouse.when.who)) +
  geom_point(size = 5, show.legend = F, alpha = .2) +
  facet_grid(State ~ paste("[U-13C]", Infusate), scales = "free") +
  coord_flip() +
  theme_bw(base_size = 15)  +
  theme(axis.text.x = element_text(angle = 50, hjust = 1),
        axis.text.y = element_text(size = 10))

# 5.2
d.compoundSelected_ordered %>% 
  select(State, Infusate, infuse.nmol.min.g, mouse.when.who) %>% distinct() %>% 
  ggplot(aes(x = paste("[U-13C]", Infusate), y = infuse.nmol.min.g, color = mouse.when.who)) +
  geom_point(size = 5, show.legend = F, shape = 21, position = position_jitter(.3, 0)) +
  facet_grid(State ~ ., scales = "free") +
  coord_flip() +
  theme_bw(base_size = 15)  +
  theme(axis.text.x = element_text(angle = 50, hjust = 1),
        axis.text.y = element_text(size = 10))


# summarizing tracer-dose combination
f <- function(d){d %>% select(State, Infusate, infuse.nmol.min.g) %>% distinct()}
f(d = d.compoundSelected_ordered %>% filter(State == "fasted"))
f(d = d.compoundSelected_ordered %>% filter(State == "refed"))

# here we normalize the labeling to the average infusion rate for each tracer
# set up the standard benchmark infusion rate: the max of infusion rate for reach tracer at each state
d.referenceInfusionRate <- d.compoundSelected_ordered %>% 
  select(State, Infusate, infuse.nmol.min.g) %>% distinct() %>% 
  group_by(State, Infusate) %>% 
  filter(infuse.nmol.min.g == max(infuse.nmol.min.g)) %>% 
  rename(infuse.benchmark = infuse.nmol.min.g)

# normalize M+1, M+2, ....using the benchmark infusion rate
d.infusionRate.benchmark <- d.compoundSelected_ordered %>% 
  left_join(d.referenceInfusionRate) %>% 
  # overwrite the original labeling with the newly scaled labeling adjusted to the benchmark infusion rate
  # most labeling is not changed
  mutate(labeling = ifelse(C_Label !=0, infuse.benchmark / infuse.nmol.min.g * labeling, NA)) %>% 
  
  # adjust the M+0 by subtracting from 1 with the updated labeling of M+1, M+2...
  group_by(mouse.when.who, tissue, Compound) %>% 
  mutate(labeling = ifelse(C_Label == 0, 1 - sum(labeling, na.rm = T), labeling)) %>% 
  # keep the benchmark infusion data, and discard the original infusion rate, concentration, and infusion rate
  select(-c(infuse.nmol.min.g, `Concentration (mM)`, `Rate (µl/min/g)`)) %>% 
  rename(infuse.nmol.min.g = infuse.benchmark) %>% 
  ungroup()

d.infusionRate.benchmark


# calculate normalized labeling (for heatmap display purpose)
# basis of normalization : fully labeled tracer enrichment
d.labeling.blood.Tracer_Cmax <- d.infusionRate.benchmark %>% 
  mutate(tracee = as.character(Compound) == as.character(Infusate)) %>% 
  filter(tissue == "Blood") %>% 
  filter(tracee == T) %>% 
  group_by(Compound) %>% 
  filter(C_Label == max(C_Label)) %>% 
  group_by(State, Infusate, mouse.when.who) %>% 
  summarise(labeling.tracerMax = unique(labeling))

d.labeling.blood.Tracer_Cmax

d.infusionRate.benchmark.normalized <- 
  d.infusionRate.benchmark %>% 
  left_join(d.labeling.blood.Tracer_Cmax, by = c("mouse.when.who", "State", "Infusate")) %>% 
  mutate(labeling.norm = labeling / labeling.tracerMax) %>% 
  
  # not include in modeling
  filter(Infusate != "citrate") %>% 
  filter(! Compound %in%  c("citrate", "a-ketoglutarate")) # %>% 
# filter(! tissue %in% c("gW", "iW", "H", "Ln" , "S"))
# filter(! tissue %in% c("gW", "iW", "H", "Ln" , "S"))



# p5: heatmap of labeling
myColors <- colorRampPalette(
  c("grey90", "#d1e6f7", "#a6d3f1", "#57b7e7", "#1bab90", "#3aab70", 
    "#fdd53e", "orange", "firebrick2", "firebrick4", "black"),
  bias = 3)(100) 

d.infusionRate.benchmark.normalized$Infusate %>% unique()





func.heatmap.labeling <- function(
    bloodOrTissue = "Blood", 
    myState = "fasted"){
  
  max.labeling <- .2 # set an upper max labeling beyond which the same deep saturated color is used
  
  if (bloodOrTissue == "Blood") {
    x <- d.infusionRate.benchmark.normalized %>% 
      filter(tissue == "Blood")
  } else {
    x <- d.infusionRate.benchmark.normalized %>% 
      filter(tissue != "Blood")
  }
  
  y <- x %>% filter(State == myState) %>% 
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
    ggplot(aes(x = mouse.when.who, y = Compound_C_Label, fill = labeling)) +
    geom_tile(color = "white") +
    # geom_text(aes(label = round(labeling, 3)), size = .5) +
    # geom_raster() +
    facet_grid(tissue ~ Infusate, scales = "free_x", space = "free") +
    # scale_fill_distiller(palette = "Spectral") +
    scale_fill_gradientn(colours = myColors, 
                         breaks = seq(from = 0, to = .5, by = .025), 
                         limits = c(0, max.labeling), # ensure the same color scale for the specified range of data
                         values = seq(0, 1, length.out = length(myColors))) +
    theme_classic() +
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
      panel.spacing = unit(0, "pt")
    ) +
    scale_x_discrete(expand = expansion(add = 0)) +
    scale_y_discrete(expand = expansion(add = 0)) +
    guides(fill = guide_colorbar(barheight = unit(200, "pt"))) 
  
  
  # add color bar denoting compounds
  p.sideBar <- y %>% select(Compound_C_Label, Compound, tissue) %>% 
    
    mutate(whichCompound = ifelse(Compound == "a-ketoglutarate", "aKG", as.character(Compound))) %>% 
    # if the isotope label is 1, or 16, 18 for fatty acids, then put the compund label there
    mutate(whichCompound = ifelse(str_extract(Compound_C_Label, "\\d{1,2}$") %in% c("1", "16", "18"), 
                                  paste(as.character(whichCompound), "   "), NA)) %>% 
    
    ggplot(aes(x = 1, y = Compound_C_Label, fill = Compound, color = Compound)) +
    geom_tile(color = "white", size = .5) +
    
    # label isotopologue index M+1, M+2....
    geom_text(aes(label = Compound_C_Label %>% str_extract("\\d{1,2}$"), x = 0),
              size = 1.6, hjust = 1) +
    
    # label compound name at M+1 position
    geom_text(aes(label = whichCompound, x = -.1),
              size = 2.1, hjust = 1) +
    
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
x1 <- func.heatmap.labeling(myState = "fasted", bloodOrTissue = "Blood")   #; x1
x2 <- func.heatmap.labeling(myState = "fasted", bloodOrTissue = "Tissues") #; x2
cowplot::plot_grid(x1 + theme(plot.margin = margin(b = 0)), 
                   ggplot() + theme_void(),
                   x2 + theme(strip.text.x.top = element_blank(),
                              plot.margin = margin(t = 0)), 
                   ncol = 1, 
                   align = "v",
                   rel_heights = c(1, .12,  2.3*1.8))

ggsave("fasted labeling.pdf", height = 10*1.3, width = 5)
 


# refed
y1 <- func.heatmap.labeling(myState = "refed", bloodOrTissue = "Blood")   # ; y1
y2 <- func.heatmap.labeling(myState = "refed", bloodOrTissue = "Tissues") # ; y2
cowplot::plot_grid(y1 + theme(plot.margin = margin(b = 0), axis.ticks.x = element_blank()), 
                   ggplot() + theme_void(),
                   y2 + theme(strip.text.x.top = element_blank(),
                              plot.margin = margin(t = 0)), 
                   ncol = 1, 
                   align = "v", 
                   rel_heights = c(1, .12,  2.3*1.8))
ggsave("refed labeling.pdf", height = 10*1.3, width = 5)



# 
# d.infusionRate.benchmark.normalized %>% 
#   filter(tissue == "Blood") %>% 
#   filter(Compound == "3-HB") %>% view()



# plot heatmap horizontally
func.heatmap.labeling <- function(
    bloodOrTissue = "Blood", 
    myState = "fasted"){
  
  max.labeling <- .2 # set an upper max labeling beyond which the same deep saturated color is used
  
  if (bloodOrTissue == "Blood") {
    x <- d.infusionRate.benchmark.normalized %>% 
      filter(tissue == "Blood")
  } else {
    x <- d.infusionRate.benchmark.normalized %>% 
      filter(tissue != "Blood")
  }
  
  y <- x %>% filter(State == myState) %>% 
    # keep only the full carbon number in FAs
    filter(! (Compound %in% "palmitate" & C_Label != 16)) %>% 
    filter(! (Compound %in% "oleate" & C_Label != 18)) %>% 
    filter(! (Compound %in% "linoleate" & C_Label != 18)) %>% 
    
    
    # compound isotopologues
    mutate(Compound_C_Label = paste(Compound, C_Label), .after = 2) %>% 
    mutate(Compound_C_Label = factor(Compound_C_Label, levels = .$Compound_C_Label %>% unique() )) %>% 
    filter(C_Label != 0) %>% 
    
    # use the same deep saturated color for very big labeling
    mutate(labeling = ifelse(labeling > max.labeling, max.labeling, labeling))
  
  p.main <- y %>% 
    ggplot(aes(y = mouse.when.who, x = Compound_C_Label, fill = labeling)) +
    geom_tile(color = "white") +
    # geom_text(aes(label = round(labeling, 3)), size = .5) +
    # geom_raster() +
    facet_grid(Infusate ~ tissue, scales = "free_y", space = "free", switch = "both") +
    # scale_fill_distiller(palette = "Spectral") +
    scale_fill_gradientn(colours = myColors, 
                         breaks = seq(from = 0, to = .5, by = .025), 
                         limits = c(0, max.labeling), # ensure the same color scale for the specified range of data
                         values = seq(0, 1, length.out = length(myColors))) +
    theme_classic() +
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
      panel.spacing = unit(0, "pt")
    ) +
    scale_x_discrete(expand = expansion(add = 0)) +
    scale_y_discrete(expand = expansion(add = 0)) +
    guides(fill = guide_colorbar(barheight = unit(200, "pt"))) 
  
  
  # add color bar denoting compounds
  p.sideBar <- y %>% select(Compound_C_Label, Compound, tissue) %>% 
    # filter(Compound == "3-HB") %>% 
    
    mutate(whichCompound = ifelse(Compound == "a-ketoglutarate", "aKG", as.character(Compound))) %>% 
    # if the isotope label is 1, or 16, 18 for fatty acids, then put the compound label there
    mutate(whichCompound = ifelse(str_extract(Compound_C_Label, "\\d{1,2}$") %in% c("1", "16", "18"), 
                                  paste(as.character(whichCompound), "   "), NA)) %>% 
    
    ggplot(aes(y = 1, x = Compound_C_Label, fill = Compound, color = Compound)) +
    geom_tile(color = "white", size = .5) +
    
    # label isotopologue index M+1, M+2....at the END of the isotopologue name 
    geom_text(aes(label = Compound_C_Label %>% str_extract("\\d{1,2}$"), y = 2.5),
              size = 1.7) +
    
    # label compound name at M+1 position
    geom_text(aes(label = whichCompound, y = 3.5),
              size = 2.2, hjust = 0, angle = 90) +
    
    facet_grid(.~tissue, switch = "both") +
    scale_y_continuous(expand = expansion(mult = 0)) +
    coord_cartesian(clip = "off") +
    theme_minimal() +
    theme(
      plot.margin = margin(t = 30, unit = "pt"),
      strip.background = element_blank(),
      strip.text = element_blank(),
      legend.position = "none", 
      panel.spacing = unit(0, "pt"), # needs to be identical as the main plot
      axis.title = element_blank(),
      axis.ticks.x = element_blank(),
      panel.grid = element_blank(),
      axis.text = element_blank(),
      # axis.text.y = element_text(size = 6, margin = margin(r = 40, unit = "pt")
    ) +
    scale_fill_brewer(palette =  ifelse(bloodOrTissue == "Blood", "Paired", "Set2")) +
    scale_color_brewer(palette =  ifelse(bloodOrTissue == "Blood", "Paired", "Set2")) 
  
  (p.sideBar / p.main) + plot_layout(heights = c(.4, 8))
  
}


# fasted
x1 <- func.heatmap.labeling(myState = "fasted", bloodOrTissue = "Blood")   #; x1
x2 <- func.heatmap.labeling(myState = "fasted", bloodOrTissue = "Tissues") #; x2
cowplot::plot_grid(x1 + theme(plot.margin = margin(r = 35), legend.position = "none"), 
                   x2 + theme(# strip.text.y.left = element_blank(),
                     plot.margin = margin(l = 0)), 
                   nrow = 1, align = "h", rel_widths = c(.7, 2.3))

ggsave("fasted labeling horizontal.pdf", height = 4, width = 12)








# labeling heatmap of blood and liver alone
func.heatmap.labeling <- function(
    bloodOrTissue = "Blood", 
    myState = "fasted"){
  
  max.labeling <- .2 # set an upper max labeling beyond which the same deep saturated color is used
  
  if (bloodOrTissue == "Blood") {
    x <- d.infusionRate.benchmark.normalized %>% 
      filter(tissue == "Blood")
  } else {
    x <- d.infusionRate.benchmark.normalized %>% 
      filter(tissue != "Blood") %>% 
      filter(tissue == "Lv")
  }
  
  y <- x %>% filter(State == myState) %>% 
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
    ggplot(aes(x = mouse.when.who, y = Compound_C_Label, fill = labeling)) +
    geom_tile(color = "white") +
    # geom_text(aes(label = round(labeling, 3)), size = .5) +
    # geom_raster() +
    facet_grid(tissue ~ Infusate, scales = "free_x", space = "free") +
    # scale_fill_distiller(palette = "Spectral") +
    scale_fill_gradientn(colours = myColors, 
                         breaks = seq(from = 0, to = .5, by = .025), 
                         limits = c(0, max.labeling), # ensure the same color scale for the specified range of data
                         values = seq(0, 1, length.out = length(myColors))) +
    theme_classic() +
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
      panel.spacing = unit(0, "pt")
    ) +
    scale_x_discrete(expand = expansion(add = 0)) +
    scale_y_discrete(expand = expansion(add = 0)) +
    guides(fill = guide_colorbar(barheight = unit(200, "pt"))) 
  
  
  # add color bar denoting compounds
  p.sideBar <- y %>% select(Compound_C_Label, Compound, tissue) %>% 
    
    mutate(whichCompound = ifelse(Compound == "a-ketoglutarate", "aKG", as.character(Compound))) %>% 
    # if the isotope label is 1, or 16, 18 for fatty acids, then put the compund label there
    mutate(whichCompound = ifelse(str_extract(Compound_C_Label, "\\d{1,2}$") %in% c("1", "16", "18"), 
                                  paste(as.character(whichCompound), "   "), NA)) %>% 
    
    ggplot(aes(x = 1, y = Compound_C_Label, fill = Compound, color = Compound)) +
    geom_tile(color = "white", size = .5) +
    
    # label isotopologue index M+1, M+2....
    geom_text(aes(label = Compound_C_Label %>% str_extract("\\d{1,2}$"), x = 0),
              size = 2.5, hjust = 1) +
    
    # label compound name at M+1 position
    geom_text(aes(label = whichCompound, x = -.1),
              size = 3.3, hjust = 1) +
    
    facet_grid(tissue~.) +
    scale_x_continuous(expand = expansion(mult = 0)) +
    coord_cartesian(clip = "off") +
    theme_minimal() +
    theme(
      plot.margin = margin(l = 70, r = 5, unit = "pt"),
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
x1 <- func.heatmap.labeling(myState = "fasted", bloodOrTissue = "Blood")   #; x1
x2 <- func.heatmap.labeling(myState = "fasted", bloodOrTissue = "Tissues") #; x2
cowplot::plot_grid(x1 + theme(plot.margin = margin(b = 0)), 
                   ggplot() + theme_void(),
                   x2 + theme(strip.text.x.top = element_blank(),
                              plot.margin = margin(t = 0)), 
                   ncol = 1, 
                   align = "v",
                   rel_heights = c(2, .2,  .9))

ggsave("fasted labeling_Lv.pdf", height = 6.2, width = 5)








 # Fatty acids are barely labeled unless infused; only labeling data of FAs in blood is used
# turn fatty acids into C2 equivalents: consider M+16 / m+18 as M+2, M+1 as 0, and M+0 is unchanged
# plot fatty acids
f.plt_FAs <- function(myData, myState){
  myData %>% 
    filter(State == myState) %>% 
    filter(Infusate %in% c("palmitate", "oleate", "linoleate")) %>% 
    filter(Compound %in% c("palmitate", "oleate", "linoleate")) %>% 
    filter(tissue == "Blood") %>% 
    ggplot(aes(x = mouse.when.who , y = labeling, fill = factor(C_Label))) +
    geom_col(position = "stack", color = "black") +
    geom_text(aes(label = C_Label), position = position_stack(vjust = .5), size = 3, color = "white") +
    facet_grid(Compound~Infusate, scales = "free") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
}
f.plt_FAs(myData = d.infusionRate.benchmark, myState = "fasted")
f.plt_FAs(myData = d.infusionRate.benchmark, myState = "refed")

f.FAs_2C <- function(whichFA, C){
  d.infusionRate.benchmark %>% 
    filter(Compound == c(whichFA)) %>% 
    filter(C_Label %in% c(0, 1, C)) %>% 
    mutate(C_Label = ifelse(C_Label == C, 2, C_Label)) # turn (M+full label) to M+2
}

x.pal <-  f.FAs_2C(whichFA = "palmitate", C = 16) 
x.ole <-  f.FAs_2C(whichFA = "oleate",    C = 18) 
x.lino <- f.FAs_2C(whichFA = "linoleate", C = 18) 

d.clean <- d.infusionRate.benchmark %>% 
  filter(! Compound %in% c("palmitate", "oleate", "linoleate")) %>% 
  bind_rows(x.pal, x.ole, x.lino) %>% 
  
  # adjust the FAs infusion rate to C2 equivalent
  mutate(infuse.nmol.min.g    = ifelse(Infusate == "palmitate",  infuse.nmol.min.g   * 16/2, infuse.nmol.min.g),
         infuse.nmol.min.g    = ifelse(Infusate == "oleate",     infuse.nmol.min.g   * 18/2, infuse.nmol.min.g),
         infuse.nmol.min.g    = ifelse(Infusate == "linoleate",  infuse.nmol.min.g   * 18/2, infuse.nmol.min.g)) %>% 
  # arrange in order
  arrange(State, Infusate, mouse.when.who, tissue, Compound) 

d.clean



# check the updated profile
f.plt_FAs(myData = d.clean, myState = "fasted")
f.plt_FAs(myData = d.clean, myState = "refed")

d.clean %>% select(State, Infusate, infuse.nmol.min.g) %>% distinct()
d.clean  %>% select(Infusate, infuse.nmol.min.g) %>% distinct()


# update the naming of compounds and infusate to match that used in atom transition list
# Compound
d.clean$Compound %>% unique()
x        <- c("Glc",     "Lac",     "Ala",     "Gln",       "Palm",      "Ole",    "Lino",      "Cit",     "Mal",    "Suc",       "Glu",       "Glycerol",  "aKG",             "HB")
names(x) <- c("glucose", "lactate", "alanine", "glutamine", "palmitate", "oleate", "linoleate", "citrate", "malate", "succinate", "glutamate", "glycerol",  "a-ketoglutarate", "3-HB")

# tracer
d.clean$Infusate %>% unique()
y        <- c("13CGlc",  "13CLac",  "13CAla",  "13CGln",    "13CPalm",   "13COle", "13CLino",   "13CCit",  "13CGlycerol", "13CHB")
names(y) <- c("glucose", "lactate", "alanine", "glutamine", "palmitate", "oleate", "linoleate", "citrate", "glycerol",    "3-HB")

d.final <- d.clean %>% 
  # !!! Important to first convert to character. As factor, the underlying integer value is otherwise used for subsetting!!! 
  mutate(Compound = x[as.character(Compound)]) %>% 
  mutate(Infusate = y[as.character(Infusate)]) %>% 
  mutate(Compound.tissue = str_c(Compound, ".", tissue))


d.final$Compound %>% unique()
d.final$Infusate %>% unique()

# d.clean %>% filter(Infusate == "citrate") %>% duplicated() %>% sum()


# add CO2 labeling data
d.fox <- read_excel("suppl_data.xlsx", sheet = "fox") %>% 
  select(State, Infusate, fox.mean) 
d.fox



colnames(d.final)
d.CO2.M1 <- d.final %>% 
  # calculate 13C carbon atom infusion rate
  
  # here add tissue == "Blood" so that 13Ccitrate tracer will not contain duplicated rows
  # different from other tracers, citrate tracer is the only one that appears in both blood and tissues
  filter(str_remove(Infusate, "13C") == Compound & tissue == "Blood") %>%  # tracer == tracee
  group_by(Infusate) %>% 
  filter(C_Label == max(C_Label)) %>% 
  mutate(infuse.13C.atom.nmol.g = infuse.nmol.min.g * C_Label) %>% 
  
  # mouse wise infusion rate 
  select(State, Infusate, mouse.when.who, infuse.nmol.min.g, infuse.13C.atom.nmol.g) %>% 
  
  # combine with fox data
  left_join(d.fox, by = c("State", "Infusate")) %>% 
  
  # calculate expected CO2 labeling
  mutate(
    Compound = "CO2",
    tissue = "Blood",
    Compound.tissue = "CO2.Blood",
    # fasted - 1800 nmol / min / g BW total CO2 production flux; refed, 20% higher 
    labeling = infuse.13C.atom.nmol.g * fox.mean / ifelse(State == "fasted", 1800, 1800 * 1.2) 
  ) %>%    
  mutate(C_Label = 1)

# d.final  %>% filter(Infusate == "13CCit") %>% duplicated() %>% sum()
# d.CO2.M1 %>% filter(Infusate == "13CCit") %>% duplicated() %>% sum()

# create M+0 labeling
d.CO2.M0 <- d.CO2.M1 %>% mutate(C_Label = 0, labeling = 1 - labeling) 

# combine M+0 and M+1
d.CO2 <- bind_rows(d.CO2.M0, d.CO2.M1) %>% 
  select("Compound", "C_Label", "tissue", "labeling", "State", "Infusate", "mouse.when.who", "infuse.nmol.min.g", "Compound.tissue") %>% 
  arrange(State, Infusate, mouse.when.who)

d.CO2

d.ready <- d.final %>% bind_rows(d.CO2) %>% 
  arrange(State, Infusate, mouse.when.who, tissue, C_Label)



# Add a simplified mouse ID
d.baked <- d.ready %>% group_by(mouse.when.who) %>% 
  mutate(m.id = mouse.when.who %>% as.character() %>% as.factor()) %>% 
  mutate(m.id = str_c("m", as.integer(m.id)))
# mutate(m.id = str_c("M", mouse.when.who %>% as.character() %>% factor() %>% as.integer()) )



# add carbon sequence (full EMU)
d.fresh <- d.baked %>% 
  group_by(Compound) %>% 
  mutate(C.max = max(C_Label)) %>% 
  rowwise() %>% 
  mutate(Compound.tissue.seq      = str_c(Compound.tissue, "_", str_c(1:C.max, collapse = ""))) %>% 
  mutate(Compound.tissue.seq_m.id = str_c(Compound.tissue.seq, "|", m.id)) %>% 
  ungroup()



# remove FAs and glycerol from the regression fit, as they are barely labeled 
# (e.g. DNL not included in model, negligible), unless labeled by its own tracer

d.meaningful <- d.fresh %>% 
  filter(! (Infusate %in% c("13CGlc", "13CLac", "13CAla", "13CGln") & Compound %in% c("Glycerol", "Palm", "Ole", "Lino"))) %>% 
  filter(! (Infusate %in% c("13CPalm")                              & Compound %in% c("Glycerol", "Ole", "Lino"))) %>% 
  filter(! (Infusate %in% c("13COle")                               & Compound %in% c("Glycerol", "Palm", "Lino"))) %>% 
  filter(! (Infusate %in% c("13CLino")                              & Compound %in% c("Glycerol", "Palm", "Ole"))) %>% 
  filter(! (Infusate %in% c("13CGlycerol")                          & Compound %in% c("Palm", "Ole", "Lino"))) %>% 
  
  arrange(Infusate, mouse.when.who, tissue, Compound, C_Label) # this is the final determinant of row arrangement order


d.meaningful$Compound %>% unique()

# # testing if glucose and lactate labeling is two fold higher
# x1 <- d.meaningful %>% filter(Infusate %in% c("13CGlc", "13CLac") ) %>% 
#   mutate(labeling = ifelse(C_Label != 0, labeling * 2, NA)) %>% 
#   group_by(State, Infusate, mouse.when.who, tissue, Compound) %>% 
#   mutate(labeling = ifelse(C_Label == 0, 1 - sum(labeling, na.rm = T), labeling))
# 
# x2 <- d.meaningful %>% filter(! Infusate %in% c("13CGlc", "13CLac") )
# d.meaningful <- bind_rows(x1, x2)



# Calculate standard deviation for each isotopologue
d.good <- d.meaningful %>% group_by(State, Infusate, tissue, Compound, C_Label) %>% 
  mutate(labeling.sd = sd(labeling) %>% round(7), .after = labeling) %>% 
  ungroup() %>% 
  # if there is a single replicate, there is no SD; set it to the mean SD of all labeling
  mutate(labeling.sd = ifelse(is.na(labeling.sd), mean(labeling.sd, na.rm = T), labeling.sd)) %>% 
  # replace zero standard deviation with some arbitrary number (mean of the standard deviations of the rest of the data)  
  mutate(labeling.sd = ifelse(labeling.sd == 0, mean(labeling.sd, na.rm = T), labeling.sd))








# as the glutamate is used as labeling proxy for aKG, 
# remove the actual aKG labeling (usually of low intensity)from the dataset, 
# and relabel glutamate data as aKG

xx1 <- d.good %>% filter(Compound != "aKG") # remove original aKG data, where labeling is typically low and measured without enough accuracy
xx2 <- xx1 %>% filter(Compound != "Glu") # most of the data
  
xx3 <- xx1 %>% filter(Compound == "Glu") %>%  # get glutamate data, rename as aKG as its proxy 
  mutate(Compound = "aKG") %>% 
  mutate(Compound.tissue          = str_replace(Compound.tissue,          "Glu", "aKG"),
         Compound.tissue.seq      = str_replace(Compound.tissue.seq,      "Glu", "aKG"),
         Compound.tissue.seq_m.id = str_replace(Compound.tissue.seq_m.id, "Glu", "aKG"))


d.better <- bind_rows(xx2, xx3)


d.good   %>% dim()
d.better %>% dim()


# fasted and fed state each as a separate file
d.13C.fasted <- d.better %>% filter(State == "fasted")
d.13C.refed  <- d.better %>% filter(State == "refed")



# d.13C.fasted %>% filter(Compound == "Lino") %>% view()
# d.13C.refed  %>% filter(Compound == "Lino") %>% view()


# Save the final cleaned datasets
save(d.13C.fasted, d.13C.refed, file = "cleaned_labeling_data.RData")

d.13C.fasted$tissue %>% unique()
d.13C.refed$Compound %>% unique()


# 
# plot glucose labeling pattern
d.13C.fasted %>% filter(tissue == "Blood" & Compound == "Glc") %>%
  filter(Infusate == "13CGlc") %>%
  group_by(C_Label) %>%
  summarise(lab = mean(labeling),
            sd = sd(labeling)) %>%
  # filter(C_Label != 0) %>%
  ggplot(aes(x = 1.3, y = lab, fill = as.character(C_Label))) +
  geom_col(position = "stack", color = "black") +
  # scale_fill_brewer(palette = "Pastel1") +
  coord_polar(theta = "y", direction = -1, start = 13.1) +
  theme_void() +
  scale_fill_manual(values = c("steelblue1", "orange", "green2", "pink", "black", "black", "red")) +
  scale_x_continuous(limits = c(.2, 2))





# calculate average C labeling for each tracer, mouse, tissue, compound combination
d.summary.average.C.labeling <- d.clean %>% 
  group_by(State, Infusate, Compound, tissue, mouse.when.who) %>% 
  mutate(labeling.weighted = C_Label/max(C_Label) * labeling) %>% 
  summarise(average.C.labeling = sum(labeling.weighted))


# average across mice replicates
d.summary.average.C.labeling.mean <- d.summary.average.C.labeling %>% 
  group_by(State, Infusate, Compound, tissue) %>% 
  summarise(average.C.labeling = mean(average.C.labeling)) %>% 
  mutate(Infusate = factor(Infusate, levels = rev(ordered.tracers), ordered = T))

# plot
func.plt.average.C <- function(myState = "fasted", isBlood = T){
  
  D <- d.summary.average.C.labeling.mean %>% 
    filter(State == myState)
  
  if (isBlood == T) { D <- D %>% filter(tissue == "Blood") }
  if (isBlood != T) { D <- D %>% filter(tissue != "Blood") }
  
  D %>% 
    ggplot(aes(x = Compound, y = Infusate, fill = average.C.labeling)) +
    geom_tile(color = "black", linewidth = .5) +
    geom_text(aes(label = round(average.C.labeling*100, 1))) +
    labs(x = "labeled metabolites", y = "13C tracer infused") +
    facet_wrap(~tissue, nrow = 2) +
    scale_fill_gradientn(colours = myColors, 
                         breaks = seq(from = 0, to = .5, by = .05), 
                         limits = c(0, .25), # ensure the same color scale for the specified range of data
                         values = seq(0, 1, length.out = length(myColors)),
                         labels = ~. * 100,
                         name = "average C labeling %") +
    guides(fill = guide_colorbar(
      barheight = unit(350, "pt"),
      barwidth = unit(10, "pt"),
      title.position = "left",
      title.theme = element_text(angle = 90, hjust = .5, face = "bold"))) +
    theme(axis.text.x = element_text(angle = 40, hjust = 1),
          strip.text = element_text(size = 18),
          axis.text = element_text(size = 15)) 
  
}


# Fasted
# plot blood for p1, and other tissues for p2
p1 <- func.plt.average.C(myState = "fasted", isBlood = T) 
p2 <- func.plt.average.C(myState = "fasted", isBlood = F) 

P <- (p1 + theme(legend.position = "none") | p2 + theme(axis.title.y = element_blank())  ) + 
  plot_annotation( title = "Fasted, average C labeling",
                   theme = theme(plot.title = element_text(size = 25,hjust = 0.5, face = "bold"))) +
  plot_layout(widths = c(1, 2.8))

print(P)

ggsave(filename = "fasted Average C labeling.pdf",width = 22, height = 8)




# Refed
# plot blood for p1, and other tissues for p2
p3 <- func.plt.average.C(myState = "refed", isBlood = T) 
p4 <- func.plt.average.C(myState = "refed", isBlood = F) 

PP <- (p3 + theme(legend.position = "none") | p4 + theme(axis.title.y = element_blank())  ) + 
  plot_annotation( title = "refed, average C labeling",
                   theme = theme(plot.title = element_text(size = 25,hjust = 0.5, face = "bold"))) +
  plot_layout(widths = c(1, 2.8))

print(PP)

ggsave(filename = "refed Average C labeling.pdf",width = 22, height = 8)





# Compare the liver and muscle labeling
# d.13C.fasted %>% 
d.better %>% 
  filter(Infusate == "13CGlc") %>% 
  filter(tissue %in% c("Lv", "M", "Br")) %>% 
  filter(Compound %in% c("Mal", "Suc")) %>% 
  group_by(State, tissue, C_Label, Compound) %>% 
  summarise(labeling.mean = mean(labeling),
            labeling.sd = sd(labeling)) %>% 
  filter(C_Label != 0) %>% 
  
  # plot
  ggplot(aes(x = C_Label, y = labeling.mean, color = tissue)) +
  geom_point(position = position_dodge(.2), size = 2) +
  geom_errorbar(aes(ymin = labeling.mean - labeling.sd,
                    ymax = labeling.mean + labeling.sd),
                width = .2, position = position_dodge(.2)) +
  geom_line(position = position_dodge(.2)) +
  facet_grid(State ~ Compound, scales = "free") +
  coord_cartesian(ylim = c(0, .2)) +
  theme(strip.text = element_text(size = 16))





# FASTED STATE
# combine Tony data - TAG kinetics 
load("data_section_Fasted_TAGkinetics.RData"); d.dataForMFA.TAGkinetics.Fasted
d.13C.fasted.TAGkinetics <- d.13C.fasted %>% rbind(d.dataForMFA.TAGkinetics.Fasted)
save(d.13C.fasted.TAGkinetics, file = "../data/cleaned_labeling_data_Fasted_TAGkinetics.RData") 

# xk <- d.13C.fasted.TAGkinetics$Compound.tissue.seq_m.id
# xk[xk %>% str_detect("1516")]


# REFED STATE
# combine Tony data - hpAA data
load("data_section_hpAAs.RData"); d.dataForMFA.hpAAs
d.13C.refed.hpAA <- d.13C.refed %>% rbind(d.dataForMFA.hpAAs)
save(d.13C.refed.hpAA, file = "../data/cleaned_labeling_data_hpAA.RData") 


xk <- d.13C.refed.hpAA$Compound.tissue.seq_m.id
xk[xk %>% str_detect("1516")]



# REFED STATE
# combine Tony data - hpAA - SCFA data 
load("data_section_hpAAs.RData"); d.dataForMFA.hpAAs
load("data_section_SCFA.RData") ; d.dataForMFA.SCFA

d.13C.refed.hpAA.SCFA <- d.13C.refed %>% rbind(d.dataForMFA.hpAAs) %>% rbind(d.dataForMFA.SCFA)
save(d.13C.refed.hpAA.SCFA, file = "../data/cleaned_labeling_data_hpAA-SCFA.RData") 

xk <- d.13C.refed.hpAA.SCFA$Compound.tissue.seq_m.id
xk[xk %>% str_detect("1516")]




# REFED STATE
# combine Tony data - TAGkinetics - hpAA - SCFA data
load("data_section_hpAAs.RData")        ; d.dataForMFA.hpAAs
load("data_section_SCFA.RData")         ; d.dataForMFA.SCFA
load("data_section_TAGkinetics.RData")  ; d.dataForMFA.TAGkinetics

d.13C.refed_TAGkinetics.hpAA.SCFA <- d.13C.refed %>% 
  rbind(d.dataForMFA.TAGkinetics) %>%  
  rbind(d.dataForMFA.hpAAs) %>% 
  rbind(d.dataForMFA.SCFA)

save(d.13C.refed_TAGkinetics.hpAA.SCFA, 
     file = "../data/cleaned_labeling_data_TAGkinetics_hpAA-SCFA.RData") 


xk <- d.13C.refed_TAGkinetics.hpAA.SCFA$Compound.tissue.seq_m.id
xk[xk %>% str_detect("1516")]


d.13C.refed.hpAA.SCFA             %>% filter(Compound %in% c("Lino"))    %>% pull(Compound.tissue.seq_m.id) %>% unique()
d.13C.refed_TAGkinetics.hpAA.SCFA %>% filter(Compound %in% c("Lino"))    %>% pull(Compound.tissue.seq_m.id) %>% unique()
d.13C.refed_TAGkinetics.hpAA.SCFA %>% filter(Compound %in% c("Ole"))    %>% pull(Compound.tissue.seq_m.id) %>% unique()
d.13C.refed_TAGkinetics.hpAA.SCFA %>% filter(Compound %in% c("Palm"))    %>% pull(Compound.tissue.seq_m.id) %>% unique()
d.13C.fasted.TAGkinetics          %>% filter(Compound %in% c("Lino"))    %>% pull(Compound.tissue.seq_m.id) %>% unique()

d.13C.refed_TAGkinetics.hpAA.SCFA %>% filter(Compound %in% c("TAGLino")) %>% pull(Compound.tissue.seq_m.id) %>% unique()

d.13C.fasted.TAGkinetics          %>% filter(Compound %in% c("TAGLino")) %>% pull(Compound.tissue.seq_m.id) %>% unique()
d.13C.fasted.TAGkinetics          %>% filter(Compound %in% c("Lino"))    %>% pull(Compound.tissue.seq_m.id) %>% unique()
d.13C.fasted                      %>% filter(Compound %in% c("Lino"))    %>% pull(Compound.tissue.seq_m.id) %>% unique()
d.13C.fasted                      %>% filter(Compound %in% c("Palm"))    %>% pull(Compound.tissue.seq_m.id) %>% unique()
d.13C.fasted                      %>% filter(Compound %in% c("Ole"))    %>% pull(Compound.tissue.seq_m.id) %>% unique()
