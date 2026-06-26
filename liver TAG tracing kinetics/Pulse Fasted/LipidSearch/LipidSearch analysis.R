rm(list = ls())
library(tidyverse)
library(readxl)

theme_set(theme_classic(base_size = 14) +
            theme(strip.background = element_blank(),
                  strip.text = element_text(face = "bold"),
                  axis.text.x = element_text(angle = 0, hjust = 1)))

rstudioapi::getActiveDocumentContext()$path %>% dirname() %>% setwd()

func.clean <- function(path.data, path.ID){
  
  path.data = "LipidSearch liver lipidomics"
  path.ID = "LipidSearch liver id.xlsx"
  
  d.cell <- read_delim(file = path.data)
  d.cell.id <- read_excel(path = path.ID)
  
  # select useful columns
  d.cell.selected <- d.cell %>% select(ID:TotalGrade, contains("OrgMeanArea")) 
  
  # extract sample name
  colnames(d.cell.selected) <- colnames(d.cell.selected) %>% 
    str_remove("OrgMeanArea\\[") %>% str_remove(pattern = "\\]")
  
  # tidy up
  d.cell.tidy <- d.cell.selected %>% 
    pivot_longer(-c(ID:TotalGrade), names_to = "SampleID", values_to = "area") %>% 
    # select only grades A and B and C with MS2
    filter(TotalGrade %in% c("A", "B", "C"))
  
  d.cell.tidy$MainIon %>% unique()
  
  # summarized dataset
  d.cell.summarized <- d.cell.tidy %>% 
    
    # for the same molecule with same adduct ion (with different RT), sum up the ion intensity
    group_by(BaseRt, SampleID, LipidMolec, LipidMolecGroup, ClassKey, SubClassKey, FAKey, FAGroupKey, MainIon, TotalGrade, MolFormula) %>% 
    summarise(area = sum(area, na.rm = T)) %>% 
    
    # select only the most characteristic adduct ions
    # filter(MainIon %in% c("M+H", "M-H", "M+HCOO")) %>% 
    # if detected in both modes, select the one with highest identity
    group_by(BaseRt, SampleID, LipidMolec, LipidMolecGroup, ClassKey, SubClassKey, FAKey, FAGroupKey, TotalGrade, MolFormula) %>% 
    summarise(area = sum(area, na.rm = T))
  
  d.cell.id$SampleID %>% unique()
  d.cell.id$`File Name` %>% unique()
  
  # Clean up the data id dataset
  d.cell.id.clean <- d.cell.id %>% 
    # cleanup sample name
    mutate(biosample = str_remove(`File Name`, ".raw")) %>% 
    # correct a typo
    distinct(SampleID, biosample) 
  
  # Check the SampleID matches
  sum(!unique(d.cell.summarized$SampleID) %in% unique(d.cell.id.clean$SampleID))
  sum(!unique(d.cell.id.clean$SampleID)   %in% unique(d.cell.summarized$SampleID))
  
  # combine MS data with sample ID
  d.cell.all <- d.cell.summarized %>% ungroup() %>%  
    left_join(d.cell.id.clean, by = "SampleID", multiple = "all") %>% 
    select(-SampleID)
  
  # manuall check biosample names
  d.cell.all$biosample %>% unique()
  return(d.cell.all)
}

d <- func.clean(
  path.data = "LipidSearch liver lipidomics",
  path.ID = "LipidSearch liver id.xlsx") %>% 
  filter(! LipidMolecGroup %>% str_detect("d5-")) 

# select TAG species
d.TG <- d %>% filter(ClassKey == "TG") 

# clean TAG for isotope tracing selection
d.TG.clean <- d.TG %>% 
  filter(! LipidMolec %>% str_detect("O-")) %>% 
  filter(! LipidMolec %>% str_detect("P-")) %>% 
  filter(! LipidMolec %>% str_detect("\\+O")) %>% 
  filter(TotalGrade %in% c("A"))

d.TG.clean$LipidMolec %>% unique()


# calculate the formula with NH3 added 
d.TG.library <- d.TG.clean %>% 
  
  # select strong peaks
  group_by(BaseRt, LipidMolec, MolFormula) %>% 
  summarise(area = mean(area)) %>% 
  filter(area >  10^6) %>% 
  
  select(BaseRt, LipidMolec, MolFormula) %>% distinct() %>% 
  separate(MolFormula, sep = " ", into = c("C", "H", "O")) %>% 
  mutate(H = paste0("H", str_remove(H, "H") %>% as.integer() + 3)) %>% 
  mutate(MolFormula = paste0(C, H, O, "N"), .keep = "unused") %>% 
  rename(compound = LipidMolec,
         formula = MolFormula,
         rt = BaseRt) %>% 
  select(compound, formula, rt)

d.yes <- d.TG.library %>% filter(compound %>% str_detect("18:2"))
d.no  <- d.TG.library %>% filter(!compound %>% str_detect("18:2"))

write.csv(d.TG.library, file = "TG_NH3_lib.csv")
write.csv(d.yes,        file = "TG_NH3_18:2_yes_lib.csv")
write.csv(d.no,         file = "TG_NH3_18:2_no_lib.csv")




# calculate the mean area
d.TG.clean.meanArea <- d.TG.clean %>% 
  group_by(LipidMolec, BaseRt) %>% 
  summarise(area = mean(area)) %>% 
  ungroup() %>% 
  mutate(Lipid.Rt = paste0(LipidMolec, BaseRt)) %>% 
  mutate(LipidMolec = fct_reorder(Lipid.Rt, area, mean)) %>% 
  arrange(-area)

d.TG.clean.meanArea %>% 
  ggplot(aes(x = 1, y = area, fill = LipidMolec)) +
  geom_col(position = "fill", color = "black") +
  theme(legend.position = "none")


top10 <- d.TG.clean.meanArea[1:10, ]
top10.mol <- top10$LipidMolec
