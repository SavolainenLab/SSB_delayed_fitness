####### Create database for annual fecundity measure (All Births) #######
## Louis Bell-Roberts
## 14/01/2026

library(dplyr)
library(lubridate)
library(tidyverse)
library(nadiv)
library(readxl)
library(stringr)

# --- Repository paths ---
repo_root    <- "/SSB_delayed_fitness/Data_gen/"
supp_xlsx <- file.path(repo_root, "Data/Raw_inputs/Supplementary_Tables_Collated_v1.4.xlsx")
out_dir   <- file.path(repo_root, "Data/Cleaned_outputs")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# --- Load anonymised pedigree data (Sheet 20) ---
ped_raw <- as.data.frame(readxl::read_excel(supp_xlsx, sheet = 20)) %>%
  mutate(
    AnimalID = trimws(as.character(AnimalID)),
    DamGen   = trimws(as.character(DamGen)),
    SireGen  = trimws(as.character(SireGen)),
    Sex      = toupper(trimws(as.character(Sex))),
# Standardise pedigree IDs (reconciling typo alias BGSX -> IKVN)
    SireGen = ifelse(SireGen == "BGSX", "IKVN", SireGen)
  )

# ==============================================================================
# 1. PEDIGREE CLEANING (Structure Only)
# ==============================================================================
ped_clean <- ped_raw %>%
  select(AnimalID, DamGen, SireGen) %>%
  rename(animal = AnimalID, dam = DamGen, sire = SireGen) %>%
  mutate(across(everything(), ~trimws(toupper(as.character(.))))) %>%
  mutate(across(everything(), ~na_if(., "."))) %>%
  mutate(across(everything(), ~na_if(., ""))) %>%
  mutate(across(everything(), ~na_if(., "0"))) %>%
  mutate(across(everything(), ~na_if(., "NA")))

ped_clean <- ped_clean %>% distinct(animal, .keep_all = TRUE)

# Identify and Add Missing Parents (Founders)
all_parents <- unique(c(ped_clean$dam, ped_clean$sire))
all_parents <- all_parents[!is.na(all_parents)]
missing_parents <- all_parents[!(all_parents %in% ped_clean$animal)]

if(length(missing_parents) > 0) {
  founders <- data.frame(animal = missing_parents, dam = NA, sire = NA)
  ped_clean <- bind_rows(founders, ped_clean)
}

ped_clean_df <- ped_clean %>% 
  distinct(animal, .keep_all = TRUE) %>%
  filter(!is.na(animal)) %>%
  as.data.frame()

nadiv_ped_ready <- prepPed(ped_clean_df)

# ==============================================================================
# 2. CALCULATION: ANNUAL OFFSPRING (ALL BIRTHS)
# ==============================================================================
# Map offspring birth years
id_birth_map <- ped_raw %>%
  mutate(CleanID = trimws(toupper(as.character(AnimalID)))) %>%
  mutate(OffspringBirthYear = year(as.Date(substr(as.character(DOB), 1, 10)))) %>%
  filter(!is.na(OffspringBirthYear)) %>%
  select(CleanID, OffspringBirthYear)

# Join Birth Years to Offspring in Nadiv pedigree
ped_with_years <- nadiv_ped_ready %>%
  inner_join(id_birth_map, by = c("animal" = "CleanID"))

# Calculate annual offspring for Dams
dams_annual <- ped_with_years %>%
  filter(!is.na(dam)) %>%
  group_by(dam, OffspringBirthYear) %>%
  summarise(num_offspring_all_births = n(), .groups = "drop") %>%
  rename(Animal_ID = dam, Year = OffspringBirthYear) %>%
  mutate(sex_from_parent_role = "F")

# Calculate annual offspring for Sires
sires_annual <- ped_with_years %>%
  filter(!is.na(sire)) %>%
  group_by(sire, OffspringBirthYear) %>%
  summarise(num_offspring_all_births = n(), .groups = "drop") %>%
  rename(Animal_ID = sire, Year = OffspringBirthYear) %>%
  mutate(sex_from_parent_role = "M")

annual_counts <- bind_rows(dams_annual, sires_annual)

# Create Master List of all animals in the pedigree
all_animals <- data.frame(Animal_ID = as.character(unique(nadiv_ped_ready$animal)))

# Merge and add Sex
offspring_df <- all_animals %>%
  left_join(annual_counts, by = "Animal_ID") %>%
  left_join(
    ped_raw %>% 
      mutate(CleanID = trimws(toupper(as.character(AnimalID)))) %>%
      select(CleanID, Sex) %>% 
      distinct(), 
    by = c("Animal_ID" = "CleanID")
  ) %>%
  mutate(
    num_offspring_all_births = replace_na(num_offspring_all_births, 0),
    sex = case_when(
      !is.na(Sex) ~ Sex,
      Animal_ID %in% dams_annual$Animal_ID ~ "F",
      Animal_ID %in% sires_annual$Animal_ID ~ "M",
      TRUE ~ NA_character_
    )
  ) %>%
  select(Animal_ID, Year, sex_from_parent_role, sex, num_offspring_all_births) %>%
  filter(!is.na(Year))

# Write to file
write.csv(offspring_df, file = file.path(out_dir, "Annual_fecundity.csv"), row.names = FALSE)
cat("Annual fecundity (All Births) saved to:", file.path(out_dir, "Annual_fecundity.csv"), "\n")
