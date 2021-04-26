
library("tidyverse")
library("readxl")

HVs <- read_tsv("data/02_my_data_clean.tsv")


# Cause of Death (CoD):
# 1 = CVA/Stroke (cerebrovascular accident)
# 2 = Subarachnoidal hemorraghe/cerebral aneurysm
# 3 = Trauma
# 4 = Arrhythmia
# 5 = Myocardial infarction
# 6 = Asphyxia
# 7 = Dissection/Abdominal Aortic Aneurysm/Pulmonary Embolism
# 8 = metabolic
# 9 = Unknown


cardiovascular_CoD <- c("1", "2", "4", "5", "7")
non_cardiovascular_CoD <- c("3", "6", "8")

HVs_CoD <- HVs %>% mutate(CoD_simple = case_when(
  Cause_of_death %in% cardiovascular_CoD ~ "Cardiovascular",
  Cause_of_death %in% non_cardiovascular_CoD ~ "Other",
  Cause_of_death == "9" ~ "Unknown"
)) %>%
  rename("Atherosclerosis" = 2)

ggplot(data = HVs_CoD, aes(Age, Atherosclerosis, color = CoD_simple)) +
  geom_point() +
  ylim(0, 20)

# Write data --------------------------------------------------------------
write_tsv(x = HVs_CoD, file = "data/03_dat_aug.tsv")
