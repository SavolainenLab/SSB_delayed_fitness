####### Obtain dates of death for each monkey from anonymised pedigree #######
##Louis Bell-Roberts
#18/09/2025

library(lubridate)
library(readxl)
library(tidyverse)

# --- Repository paths ---
repo_root    <- "/Users/pn26026/OneDrive - University of Bristol/SSB/Code/Github_upload_v1.1/Data_gen/"
supp_xlsx <- file.path(repo_root, "Data/Raw_inputs/Supplementary_Tables_Collated_v1.4.xlsx")
out_dir   <- file.path(repo_root, "Data/Cleaned_outputs")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# -----------------------
# Read in anonymised pedigree (Sheet 20)
# -----------------------
pedigree <- as.data.frame(readxl::read_excel(supp_xlsx, sheet = 20)) %>%
  mutate(AnimalID = trimws(as.character(AnimalID)))

# Extract death dates directly from pedigree DOD column
pedigree_with_deaths_final <- pedigree %>%
  mutate(
    death_date = ifelse(is.na(DOD) | as.character(DOD) == "", NA_character_, substr(as.character(DOD), 1, 10))
  ) %>%
  select(AnimalID, death_date) %>%
  distinct(AnimalID, .keep_all = TRUE)

# -----------------------
# Save
# -----------------------
out_path <- file.path(out_dir, "Death_dates.csv")
write.csv(pedigree_with_deaths_final, out_path, row.names = FALSE)