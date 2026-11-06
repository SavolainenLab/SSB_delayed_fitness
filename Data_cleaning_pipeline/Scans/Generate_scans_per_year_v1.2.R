#### Count individual scan frequency per year ####
## Louis Bell-Roberts
## 19/08/2025

library(data.table)
library(tidyverse)
library(lubridate)
library(readxl)

# --- Repository paths ---
repo_root    <- "/Users/pn26026/OneDrive - University of Bristol/SSB/Code/Github_upload_v1.1/Data_gen/"
supp_xlsx <- file.path(repo_root, "Data/Raw_inputs/Supplementary_Tables_Collated_v1.4.xlsx")
out_dir   <- file.path(repo_root, "Data/Cleaned_outputs")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# --- Load anonymised scan data (Sheet 15) ---
scan_data <- as.data.frame(readxl::read_excel(supp_xlsx, sheet = 15))

# --- Read anonymised pedigree (Sheet 20) ---
pedigree <- as.data.frame(readxl::read_excel(supp_xlsx, sheet = 20)) %>%
  mutate(
    AnimalID = trimws(as.character(AnimalID)),
    Sex      = toupper(trimws(as.character(Sex)))
  )

#Identify known adult males and females from the pedigree
pedigree_females <- pedigree %>%
  filter(Sex == "F")
pedigree_males <- pedigree %>%
  filter(Sex == "M")

# Step 1: Keep date + ID columns
id_cols <- scan_data %>%
  select(DATE, starts_with("Individual"))

# Step 2: Convert DATE to proper date format and extract year
id_cols <- id_cols %>%
  mutate(DATE = parse_date_time(as.character(DATE), 
                                orders = c("ymd", "ymd HMS"), 
                                tz = "UTC"),
         Year = year(DATE))

# Keep only males
id_cols <- id_cols %>%
  filter(Individual_ID %in% pedigree_males$AnimalID)

# Step 5: Count scans per individual per year
scan_counts_by_year <- id_cols %>%
  group_by(Individual_ID, Year) %>%
  summarise(Scan_Count = n(), .groups = "drop") %>%
  arrange(Individual_ID, Year)

#Write to file
write.csv(x = scan_counts_by_year, file = file.path(out_dir, "Scans_per_year.csv"), row.names = FALSE)