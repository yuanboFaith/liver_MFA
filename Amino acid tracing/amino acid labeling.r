library(tidyverse)
library(readxl)

rm(list = ls())

(rstudioapi::getActiveDocumentContext())$path %>% dirname() %>% setwd(); getwd()

d <- read_excel("labeling.xlsx")

# subtract background
d.background.subtracted <- d %>% 
  mutate(blank = (Blank_1 + Blank_2)/2, .after = 3) %>% 
  rowwise() %>% 
  mutate(across(liver_1 : last_col(), ~. - blank), .keep = "unused") %>% 
  select(-contains("Blank")) %>% 
  ungroup()

# natural abundance correction
l.corrected <- d.background.subtracted %>% 
  accucor::natural_abundance_correction(resolution = 120000)

d.labeling <- l.corrected$Normalized


# tidy up
d.tidy <- d.labeling %>% 
  pivot_longer(-c(1:2), names_to = "sample", values_to = "labeling") %>% 
  separate(sample, into = c("tissue", "mouse.id")) 

# liver
d.liver <- d.tidy %>% filter(tissue == "liver") %>% mutate(part = "liver")
d.serum <- d.tidy %>% filter(tissue == "serum")

# serum
d.serum <- d.serum %>% mutate(part = str_extract(mouse.id, "[A-Z]"),
                              mouse.id = str_extract(mouse.id, "\\d{1,2}"))

# combine liver and serum
d.clean <- bind_rows(d.liver, d.serum) %>% 
  # use alphabet letters for mouse id, so that the mouse id is incorporated as tracer names 
  mutate(mouse.id = letters[mouse.id %>% as.integer()])

d.clean$mouse.id %>% unique()


# calculate average carbon labeling
d.average.C.labeling <- d.clean %>% 
  group_by(mouse.id, tissue, part, Compound) %>% 
  mutate(C_Max = max(C_Label)) %>% 
  mutate(label.weighted = labeling * C_Label / C_Max) %>% 
  summarise(labeling = sum(label.weighted))


# put compounds in order
EAA  <- c("Histidine", "Isoleucine", "Leucine", "Lysine", "Methionine", "Phenylalanine", "Threonine", "Tryptophan", "Valine")

NEAA <- c("Alanine", "Arginine", "Asparagine", "Aspartate", "Cysteine",
          "Glutamate", "Glutamine", "Glycine", "Proline", "Serine", "Tyrosine")

other <- c("Glucose", "Lactate", "Malate" , "Succinic acid", "3-Hydroxybutyric acid")

d.average.C.labeling <- d.average.C.labeling %>% 
  mutate(Compound = factor(Compound, levels = c(EAA, NEAA, other)))

# x <- d.clean$Compound %>% unique()
# y <- d.average.C.labeling$Compound %>% unique()
# x %in% y
# y %in% x




# plot amino acid labeling 
d.average.C.labeling %>% 
  # filter(part != "L") %>% 
  
  filter(mouse.id != 1) %>% 
  ggplot(aes(x = part, y = labeling * 100, color = part)) +
  # geom_boxplot(outlier.alpha = 0) +
  # geom_point() +
  geom_text(aes(label = mouse.id), position = position_jitter(.3, 0)) +
  # geom_line(aes(group = Compound)) +
  facet_wrap(~Compound, scales = "free") +
  theme(legend.position = "none") +
  labs (y = "label %") +
  scale_x_discrete(limits = c("liver", "P", "L", "T")) +
  theme_bw(base_size = 15) +
  scale_y_continuous(limits = c(0, NA))


# plot liver TCA labeling of isotopologues
d.liver %>% filter(Compound %in% c("Malate", "Succinic acid", "Glutamate")) %>% 
  ggplot(aes(x = mouse.id, y = labeling, fill = as.character(C_Label))) +
  geom_col(color = "black", position = position_stack(reverse = T)) +
  facet_wrap(~Compound) +
  scale_fill_brewer(palette = "Set1") +
  coord_cartesian(ylim = c(.97, 1), expand = 0)


# plot serum amino acid labeling of isotopolouges
# L = blood after cut liver; P = blood after cut portal vein; T = tail bllod; liver = liver tissue
d.clean %>% 
  # filter(part == "liver") %>% 
  # filter(Compound %in% c(EAA, NEAA)) %>% 
  group_by(Compound, C_Label, part) %>% 
  summarise(labeling = mean(labeling)) %>% 
  
  ggplot(aes(x = Compound, y = labeling, fill = C_Label)) +
  geom_col(color = "black", alpha = .7) +
  geom_text(aes(label = C_Label), position = position_stack(vjust = .5)) +
  coord_cartesian(ylim = c(.96, 1), expand = 0) +
  scale_fill_viridis_c() +
  facet_wrap(~part, nrow = 3) +
  labs(x = NULL) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(face = "bold", size = 18),
        strip.background = element_rect(fill = "beige", color = "black")) 



d.clean$Compound %>% unique()
d.clean %>% select(part, tissue) %>% table()


# select only the portal amino acids, systemic blood glucose, lactate, alanine, 3HB, and liver malate/succinate/glutamate
d.useful <- d.clean %>% filter(
  (tissue == "liver" & Compound %in% c("Malate", "Succinic acid", "Glutamate")) | # liver TCA cycle indicators
    (part == "P"  & Compound %in% c(EAA, NEAA) | # portal (P) blood for amino acids
       (part == "T"  & Compound %in% c("Glucose", "Lactate", "3-Hydroxybutyric acid", "Alanine")) # tail (T) blood for systemic circulating nutrients
    ))

d.useful$Compound %>% unique()


# simplify the names of compounds
d.metab_names <- tribble(
  ~Compound,               ~abbr,  ~C.max,
  "Glucose",               "Glc",   6,
  "3-Hydroxybutyric acid", "HB",    4,
  "Lactate",               "Lac",   3,
  "Malate",                "Mal",   4,
  "Succinic acid",         "Suc",   4,
  
  "Glutamate",             "Glu", 5,
  "Alanine",               "Ala", 3,
  "Arginine",              "Arg", 5,
  "Asparagine",            "Asn", 4,
  "Aspartate",             "Asp", 4,
  "Glutamine",             "Gln", 5,
  "Glycine",               "Gly", 2,
  "Isoleucine",            "Ile", 2,
  "Leucine",               "Leu", 2,
  "Lysine",                "Lys", 2,
  "Phenylalanine",         "Phe", 4,
  "Serine",                "Ser", 3,
  "Threonine",             "Thr", 2,
  "Tryptophan",            "Trp", 3,
  "Tyrosine",              "Tyr", 4,
  "Valine",                "Val", 4)


d.simpleNames <- d.useful %>% 
  # update compound names with acronyms
  left_join(d.metab_names, by = "Compound") %>% 
  select(-Compound) %>% rename(Compound = abbr) %>% 
  arrange(mouse.id, tissue, Compound) %>% 
  # update tissue names
  mutate(tissue = str_replace(tissue, "liver", "Lv"),
         tissue = str_replace(tissue, "serum", "Blood")) %>% 
  # update mouse id
  rename(m.id = mouse.id) %>% 
  # mutate(m.id = paste0("m", m.id)) 
  
  # for glutamate in liver, change the name to aKG (use glutamate as proxy of aKG labeling)
  mutate(Compound = ifelse(tissue =="Lv" & Compound == "Glu", "aKG", Compound)) %>% 
  # for portal blood, change the tissue name to hp (hepatic portal) - "hp" follows the network EMU naming rule
  mutate(tissue = ifelse(part == "P", "hp", tissue))

d.simpleNames %>% select(part, Compound) %>% table()
d.simpleNames %>% select(tissue, Compound) %>% table()

# combine with the major tracing data
load(file = "../data/cleaned_labeling_data.RData")


d.dataForMFA.hpAAs <- d.simpleNames %>% 
  group_by(tissue, Compound, C_Label) %>%  # there are only three compartments: liver, blood (systemic), and portal vein blood
  mutate(labeling.sd = sd(labeling)) %>% ungroup() %>% 
  
  mutate(State = "refed",
         Infusate = paste0("13ChpAA", m.id)) %>% 
  mutate(mouse.when.who = paste0(m.id, "_2026-Jan_BY")) %>% 
  mutate(infuse.nmol.min.g = NA) %>% 
  mutate(Compound.tissue          = paste0(Compound, ".", tissue)) %>% 
  
  # rowwise important to add the correct C-max
  rowwise() %>% mutate(Compound.tissue.seq      = str_c(Compound.tissue, "_", str_c(1:C.max, collapse = ""))) %>% ungroup() %>% 
  
  mutate(Compound.tissue.seq_m.id = str_c(Compound.tissue.seq, "|", m.id)) %>% 
  relocate(Compound, C_Label, tissue, labeling, labeling.sd, State, Infusate, mouse.when.who, 
           infuse.nmol.min.g, Compound.tissue, m.id, C.max, Compound.tissue.seq, Compound.tissue.seq_m.id) %>% 
  select(-part) %>% 
  
  # only mouse b, c, d, e, f, j, k, l has full set of portal blood, liver sample, and systemic blood
  filter(m.id %in% c("b", "c", "d", "e", "f", "j", "k", "l"))

d.dataForMFA.hpAAs


save(d.dataForMFA.hpAAs, file = "../data/data_section_hpAAs.RData") 






# d.organized %>% filter(Infusate == "13ChpAAa") %>% 
#   filter(tissue == "Blood") %>% 
#   filter(m.id == "a") %>% 
#   filter(Compound == "Ala") %>% view()

d.dataForMFA.hpAAs$tissue %>% unique()

# test 1
d.dataForMFA.hpAAs %>% filter(tissue == "hp") %>% 
  group_by(Compound) %>% 
  filter(C_Label == max(C_Label)) %>% 
  # filter(Compound == "Ala") %>% 
  # filter(labeling  == 0) %>% view()
  
  ggplot(aes(x = m.id, y = labeling, color = Compound)) + 
  geom_line(aes(group = Compound)) +
  ggrepel::geom_text_repel(aes(label = Compound))


# test 2
d.dataForMFA.hpAAs %>% filter(tissue == "Lv") %>% 
  filter(Infusate == "13ChpAAd")






# data for RongYa
# plot amino acid labeling 
o1 <- d.average.C.labeling %>% 
  filter(part == "P") %>%
  filter(! Compound %in% c("Asparagine", "Aspartate")) %>% 
  filter(! Compound %in% c("Glucose", "Lactate", "Malate", "Succinic acid", "3-Hydroxybutyric acid")) %>% 
  filter(mouse.id != 1) %>% 
  ungroup() %>% 
  mutate(Compound = fct_reorder(Compound, -labeling, mean))

p.portalAAs <- o1 %>% 
  ggplot(aes(x = Compound, y = labeling * 100)) +
  stat_summary(geom = "bar", fun = mean, fill = "turquoise4", alpha = .5, color = "black") +
  stat_summary(geom = "errorbar", fun.data = mean_sdl, fun.args = list(mult = 1), width = .2) +
  ggbeeswarm::geom_quasirandom(width = .2, size = 2) +
  theme(legend.position = "none") +
  labs (y = "label %") +
  theme_classic(base_size = 16) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.y = element_text(face = "bold")) +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, .1))) +
  scale_x_discrete(expand = expansion(add = 1)) +
  labs(x = NULL) +
  coord_cartesian(ylim = c(0, 5)) 

p.portalAAs




o2 <- d.average.C.labeling %>% 
  filter(part == "liver") %>%
  # filter(! Compound %in% c("Asparagine", "Aspartate")) %>% 
  filter(Compound %in% c( "Malate", "Succinic acid")) %>% 
  # filter(mouse.id != 1) %>% 
  ungroup() %>% 
  mutate(Compound = fct_reorder(Compound, -labeling, mean))

p.liverTCA <- o2 %>% 
  ggplot(aes(x = Compound, y = labeling * 100)) +
  stat_summary(geom = "bar", fun = mean, fill = "turquoise4", alpha = .5, color = "black") +
  stat_summary(geom = "errorbar", fun.data = mean_sdl, fun.args = list(mult = 1), width = .2) +
  ggbeeswarm::geom_quasirandom(width = .3, size = 2) +
  theme(legend.position = "none") +
  labs (y = "label %") +
  theme_classic(base_size = 16) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.y = element_text(face = "bold")) +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, .1))) +
  scale_x_discrete(expand = expansion(add = 1)) +
  coord_cartesian(ylim = c(0, 5)) +
  labs(x = NULL)
  
p.liverTCA

cowplot::plot_grid(p.portalAAs, p.liverTCA, nrow = 1, rel_widths = c(10, 3.5))


d.Rongya <- bind_rows(o1, o2) %>% arrange(tissue, Compound, mouse.id) %>% select(-part)
write.csv(d.Rongya, file = "13C protein tracing data.csv")
