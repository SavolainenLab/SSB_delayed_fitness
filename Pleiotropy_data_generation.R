####### Functional Benefits of SSB: Generate All Anonymised Pleiotropy Datasets #######
## Author: Louis Bell-Roberts
## Date: October 2026
## Description: Uses only the pre-anonymised focal male dataset (Sheet 19) and 
## pre-anonymised pedigree records (Sheet 20) from Supplementary_Tables_Collated_v1.4.xlsx
## to compute lifetime fecundity, birth dates, and exit dates for mothers, fathers, 
## sisters, and brothers, exporting all four clean datasets for pleiotropy modeling.
# ==============================================================================

library(tidyverse)
library(readxl)
library(lubridate)

cat("--- INGESTING PRE-ANONYMISED SUPPLEMENTARY DATASETS ---\n")

supp_path <- "/SSB_delayed_fitness/Supplementary_Tables_Collated_v1.4.xlsx"
output_dir <- "/SSB_delayed_fitness"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# 1. Read the pre-anonymised focal male dataset (Sheet 19)
focal_df <- readxl::read_excel(supp_path, sheet = 19) %>%
  dplyr::rename(Group = assigned_group) %>%
  dplyr::mutate(Animal_ID = trimws(as.character(Animal_ID)))

# 2. Read the pre-anonymised pedigree dataset (Sheet 20)
anonymised_ped <- as.data.frame(readxl::read_excel(supp_path, sheet = 20)) %>%
  dplyr::mutate(
    AnimalID = trimws(as.character(AnimalID)),
    DamGen   = trimws(as.character(DamGen)),
    SireGen  = trimws(as.character(SireGen)),
    Sex      = trimws(as.character(Sex))
  )

# 3. Pre-calculate true lifetime offspring counts across entire pedigree
female_lifetime_offspring <- anonymised_ped %>%
  dplyr::filter(!is.na(DamGen), DamGen != "", DamGen != "NA") %>%
  dplyr::count(DamGen, name = "num_offspring")

male_lifetime_offspring <- anonymised_ped %>%
  dplyr::filter(!is.na(SireGen), SireGen != "", SireGen != "NA") %>%
  dplyr::count(SireGen, name = "num_offspring")

# 4. Standardise dates and determine population exit dates
ped_dates <- anonymised_ped %>%
  dplyr::mutate(
    clean_DOB      = as.Date(substr(as.character(DOB), 1, 10)),
    clean_DOD      = as.Date(substr(as.character(DOD), 1, 10)),
    clean_Transfer = as.Date(substr(as.character(TransferDate), 1, 10)),
    # Exit date is earliest of death or physical removal/transfer
    exit_date      = pmin(clean_DOD, clean_Transfer, na.rm = TRUE)
  )

# Base parent lookup for focal males
focal_males_unique <- tibble::tibble(Animal_ID = unique(focal_df$Animal_ID))

focal_parent_links <- focal_males_unique %>%
  dplyr::left_join(
    anonymised_ped %>%
      dplyr::select(Animal_ID = AnimalID, DamGen, SireGen) %>%
      dplyr::distinct(),
    by = "Animal_ID"
  ) %>%
  dplyr::distinct(Animal_ID, .keep_all = TRUE)


# ==============================================================================
# A. MATERNAL FECUNDITY DATASET
# ==============================================================================
cat("\n--- 1. GENERATING MATERNAL FECUNDITY DATASET ---\n")

mother_demographics <- ped_dates %>%
  dplyr::select(DamGen = AnimalID, DamGen_DOB = clean_DOB, DamGen_death_date = exit_date) %>%
  dplyr::distinct(DamGen, .keep_all = TRUE) %>%
  dplyr::left_join(female_lifetime_offspring, by = "DamGen") %>%
  dplyr::mutate(
    DamGen_num_offspring = tidyr::replace_na(num_offspring, 0),
    DamGen_DOB           = as.character(DamGen_DOB),
    DamGen_death_date    = as.character(DamGen_death_date)
  ) %>%
  dplyr::select(DamGen, DamGen_DOB, DamGen_death_date, DamGen_num_offspring)

son_mother_links <- focal_parent_links %>%
  dplyr::filter(!is.na(DamGen), DamGen != "", DamGen != "NA") %>%
  dplyr::select(Animal_ID, DamGen) %>%
  dplyr::left_join(mother_demographics, by = "DamGen")

merged_mother_df <- focal_df %>%
  dplyr::left_join(son_mother_links, by = "Animal_ID")

readr::write_csv(
  merged_mother_df,
  file.path(output_dir, "mother_pleiotropy_data.csv")
)

cat("Mother dataset saved. Rows:", nrow(merged_mother_df), 
    "| Focal males:", dplyr::n_distinct(merged_mother_df$Animal_ID), 
    "| Mapped mothers:", dplyr::n_distinct(merged_mother_df$DamGen, na.rm = TRUE), "\n")


# ==============================================================================
# B. PATERNAL FECUNDITY DATASET
# ==============================================================================
cat("\n--- 2. GENERATING PATERNAL FECUNDITY DATASET ---\n")

father_demographics <- ped_dates %>%
  dplyr::select(SireGen = AnimalID, SireGen_DOB = clean_DOB, SireGen_death_date = exit_date) %>%
  dplyr::distinct(SireGen, .keep_all = TRUE) %>%
  dplyr::left_join(male_lifetime_offspring, by = "SireGen") %>%
  dplyr::mutate(
    SireGen_num_offspring = tidyr::replace_na(num_offspring, 0),
    SireGen_DOB           = as.character(SireGen_DOB),
    SireGen_death_date    = as.character(SireGen_death_date)
  ) %>%
  dplyr::select(SireGen, SireGen_DOB, SireGen_death_date, SireGen_num_offspring)

son_father_links <- focal_parent_links %>%
  dplyr::filter(!is.na(SireGen), SireGen != "", SireGen != "NA") %>%
  dplyr::select(Animal_ID, SireGen) %>%
  dplyr::left_join(father_demographics, by = "SireGen")

merged_father_df <- focal_df %>%
  dplyr::left_join(son_father_links, by = "Animal_ID")

readr::write_csv(
  merged_father_df,
  file.path(output_dir, "father_pleiotropy_data.csv")
)

cat("Father dataset saved. Rows:", nrow(merged_father_df), 
    "| Focal males:", dplyr::n_distinct(merged_father_df$Animal_ID), 
    "| Mapped fathers:", dplyr::n_distinct(merged_father_df$SireGen, na.rm = TRUE), "\n")


# ==============================================================================
# C. SISTER FECUNDITY DATASET
# ==============================================================================
cat("\n--- 3. GENERATING SISTER FECUNDITY DATASET ---\n")

sister_demographics <- ped_dates %>%
  dplyr::filter(tolower(Sex) %in% c("f", "female")) %>%
  dplyr::select(Sister_ID = AnimalID, Sister_DOB = clean_DOB, Sister_death_date = exit_date) %>%
  dplyr::distinct(Sister_ID, .keep_all = TRUE) %>%
  dplyr::left_join(female_lifetime_offspring, by = c("Sister_ID" = "DamGen")) %>%
  dplyr::mutate(
    Sister_num_offspring = tidyr::replace_na(num_offspring, 0),
    Sister_DOB           = as.character(Sister_DOB),
    Sister_death_date    = as.character(Sister_death_date)
  ) %>%
  dplyr::select(Sister_ID, Sister_DOB, Sister_death_date, Sister_num_offspring)

all_sisters <- anonymised_ped %>%
  dplyr::filter(
    tolower(Sex) %in% c("f", "female"),
    !is.na(DamGen), DamGen != "", DamGen != "NA"
  ) %>%
  dplyr::select(Sister_ID = AnimalID, DamGen) %>%
  dplyr::distinct()

brother_sister_links <- focal_parent_links %>%
  dplyr::filter(!is.na(DamGen), DamGen != "", DamGen != "NA") %>%
  dplyr::select(Animal_ID, DamGen) %>%
  dplyr::left_join(all_sisters, by = "DamGen", relationship = "many-to-many") %>%
  dplyr::filter(is.na(Sister_ID) | Sister_ID != Animal_ID) %>%
  dplyr::left_join(sister_demographics, by = "Sister_ID") %>%
  dplyr::select(Animal_ID, DamGen, Sister_ID, Sister_DOB, Sister_death_date, Sister_num_offspring)

merged_sister_df <- focal_df %>%
  dplyr::left_join(brother_sister_links, by = "Animal_ID", relationship = "many-to-many")

readr::write_csv(
  merged_sister_df,
  file.path(output_dir, "sister_pleiotropy_data.csv")
)

cat("Sister dataset saved. Rows:", nrow(merged_sister_df), 
    "| Focal males:", dplyr::n_distinct(merged_sister_df$Animal_ID), 
    "| Mapped sisters:", dplyr::n_distinct(merged_sister_df$Sister_ID, na.rm = TRUE), "\n")


# ==============================================================================
# D. BROTHER FECUNDITY DATASET
# ==============================================================================
cat("\n--- 4. GENERATING BROTHER FECUNDITY DATASET ---\n")

brother_demographics <- ped_dates %>%
  dplyr::filter(tolower(Sex) %in% c("m", "male")) %>%
  dplyr::select(Brother_ID = AnimalID, Brother_DOB = clean_DOB, Brother_death_date = exit_date) %>%
  dplyr::distinct(Brother_ID, .keep_all = TRUE) %>%
  dplyr::left_join(male_lifetime_offspring, by = c("Brother_ID" = "SireGen")) %>%
  dplyr::mutate(
    Brother_num_offspring = tidyr::replace_na(num_offspring, 0),
    Brother_DOB           = as.character(Brother_DOB),
    Brother_death_date    = as.character(Brother_death_date)
  ) %>%
  dplyr::select(Brother_ID, Brother_DOB, Brother_death_date, Brother_num_offspring)

all_brothers <- anonymised_ped %>%
  dplyr::filter(
    tolower(Sex) %in% c("m", "male"),
    !is.na(DamGen), DamGen != "", DamGen != "NA"
  ) %>%
  dplyr::select(Brother_ID = AnimalID, DamGen) %>%
  dplyr::distinct()

focal_brother_links <- focal_parent_links %>%
  dplyr::filter(!is.na(DamGen), DamGen != "", DamGen != "NA") %>%
  dplyr::select(Animal_ID, DamGen) %>%
  dplyr::left_join(all_brothers, by = "DamGen", relationship = "many-to-many") %>%
  dplyr::filter(is.na(Brother_ID) | Brother_ID != Animal_ID) %>%
  dplyr::left_join(brother_demographics, by = "Brother_ID") %>%
  dplyr::select(Animal_ID, DamGen, Brother_ID, Brother_DOB, Brother_death_date, Brother_num_offspring)

merged_brother_df <- focal_df %>%
  dplyr::left_join(focal_brother_links, by = "Animal_ID", relationship = "many-to-many")

readr::write_csv(
  merged_brother_df,
  file.path(output_dir, "brother_pleiotropy_data.csv")
)

cat("Brother dataset saved. Rows:", nrow(merged_brother_df), 
    "| Focal males:", dplyr::n_distinct(merged_brother_df$Animal_ID), 
    "| Mapped brothers:", dplyr::n_distinct(merged_brother_df$Brother_ID, na.rm = TRUE), "\n")

cat("\nPipeline complete! All 4 datasets generated using ONLY anonymised inputs.\n")
