####### Create database for lifetime fecundity measure #######
##Louis Bell-Roberts
#09/09/2025

library(dplyr)
library(lubridate)
library(tidyverse)
library(nadiv)
library(readxl)

# --- Repository paths ---
repo_root <- "/Users/pn26026/Documents/Sandisk/SSB_delayed_fitness"
supp_xlsx <- file.path(repo_root, "Data/Raw_inputs/Supplementary_Tables_Collated_v1.4.xlsx")
out_dir   <- file.path(repo_root, "Data/Cleaned_outputs")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# --- Load anonymised pedigree data (Sheet 20) ---
pedigree <- as.data.frame(readxl::read_excel(supp_xlsx, sheet = 20)) %>%
  mutate(
    AnimalID = trimws(as.character(AnimalID)),
    DamGen   = trimws(as.character(DamGen)),
    SireGen  = trimws(as.character(SireGen)),
    Sex      = toupper(trimws(as.character(Sex)))
  ) %>%
  filter(!is.na(AnimalID), AnimalID != "", AnimalID != "NA") %>%
  distinct(AnimalID, .keep_all = TRUE)

# Identify individuals surviving to at least 1 year of age (replaces 3-digit raw tattoo check)
surviving_1yr_ids <- pedigree %>%
  mutate(
    clean_DOB = as.Date(substr(as.character(DOB), 1, 10)),
    clean_DOD = as.Date(substr(as.character(DOD), 1, 10)),
    age_death = as.numeric(difftime(clean_DOD, clean_DOB, units = "days")) / 365.25
  ) %>%
  filter(is.na(age_death) | age_death >= 1) %>%
  pull(AnimalID)

# --- Prepare pedigree for nadiv ---
ped <- pedigree %>%
  select(AnimalID, DamGen, SireGen) %>%
  rename(animal = AnimalID, dam = DamGen, sire = SireGen)

nadiv_ped <- prepPed(ped)

# --- Calculate reproductive output for each individual ---
nadiv_ped_filtered <- nadiv_ped %>%
  filter(as.character(animal) %in% surviving_1yr_ids)

# Count offspring for dams and sires
offspring_per_dam  <- table(nadiv_ped_filtered$dam)
offspring_per_sire <- table(nadiv_ped_filtered$sire)

offspring_dam_df <- as.data.frame(offspring_per_dam)
names(offspring_dam_df) <- c("Animal_ID", "num_offspring")
offspring_dam_df$sex <- "F"

offspring_sire_df <- as.data.frame(offspring_per_sire)
names(offspring_sire_df) <- c("Animal_ID", "num_offspring")
offspring_sire_df$sex <- "M"

# Combine male and female reproductive output
offspring_df <- bind_rows(offspring_dam_df, offspring_sire_df)

# Include all individuals in the pedigree with zero offspring
all_animals <- data.frame(Animal_ID = as.character(unique(nadiv_ped$animal)))

offspring_df <- all_animals %>%
  left_join(offspring_df, by = "Animal_ID") %>%
  left_join(pedigree %>% select(AnimalID, Sex), by = c("Animal_ID" = "AnimalID")) %>%
  mutate(
    num_offspring = if_else(is.na(num_offspring), 0L, as.integer(num_offspring)),
    sex = case_when(
      !is.na(Sex) ~ Sex,
      Animal_ID %in% offspring_dam_df$Animal_ID ~ "F",
      Animal_ID %in% offspring_sire_df$Animal_ID ~ "M",
      TRUE ~ NA_character_
    )
  ) %>%
  select(-Sex)

#Write to file
write.csv(offspring_df, file = file.path(out_dir, "Lifetime_fecundity.csv"), row.names = FALSE)