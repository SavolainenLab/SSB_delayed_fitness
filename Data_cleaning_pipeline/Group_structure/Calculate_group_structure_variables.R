####### Group structure: compute per-year group summaries over behaviour-sampled months #######
## Author: Louis Bell-Roberts
## Date: 30/09/2025

library(tidyverse)
library(lubridate)
library(janitor)
library(readxl)
library(stringr)

# -------------------- repository paths --------------------
repo_root           <- "/Users/pn26026/Documents/Sandisk/SSB_delayed_fitness"
raw_monthly_summary <- file.path(repo_root, "Data/Raw_inputs/group_structure_summary_LBR.csv")
input_dir           <- file.path(repo_root, "Data/Raw_inputs/Census_cleaned")
out_dir             <- file.path(repo_root, "Data/Cleaned_outputs")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

output_file_monthly <- file.path(out_dir, "group_structure_summary_LBR.csv")
output_file_annual  <- file.path(out_dir, "group_structure_annual_means_filtered_LBR.csv")

# -------------------- load monthly summary or compute from census files --------------------
if (file.exists(raw_monthly_summary)) {
  message("Reading pre-summarised monthly group structure table (no raw census files required)...")
  summary_tbl <- read.csv(raw_monthly_summary, stringsAsFactors = FALSE)
} else {
  stop("Please place 'group_structure_summary_LBR.csv' inside Data/Raw_inputs/ so raw census files are not required.")
}

# -------------------- behaviour-sampled months filter --------------------
# Only include these (year, month) when computing annual averages.
valid_months <- tibble(
  year = c(rep(2017,4), rep(2019,4), rep(2020,2), rep(2023,3), rep(2024,3), rep(2025,3)),
  month_num = c(3,4,5,6,  2,3,4,5,  2,3,  3,4,5,  2,3,4,  1,2,3)
)

summary_tbl_filtered <- summary_tbl %>%
  semi_join(valid_months, by = c("year", "month_num"))

# -------------------- per-year (per group) averages over behaviour months --------------------
annual_means <- summary_tbl_filtered %>%
  group_by(year, group) %>%
  summarise(
    across(
      c(
        group_size,
        n_males, n_females,
        n_male_adult, n_male_subadult, n_male_juvenile,
        n_female_adult, n_female_juvenile, n_female_prejuvenile,
        prop_males, prop_females,
        prop_male_adult, prop_male_subadult, prop_male_juvenile,
        prop_female_adult, prop_female_juvenile, prop_female_prejuvenile
      ),
      ~ mean(.x, na.rm = TRUE)
    ),
    n_months_used = dplyr::n(),
    .groups = "drop"
  ) %>%
  mutate(
    # 1) Adult : Subadult males
    ratio_male_adult_to_subadult =
      ifelse(n_male_subadult > 0, n_male_adult / n_male_subadult, NA_real_),
    
    # 2) (Adult + Subadult) males : (Adult + Juvenile) females
    ratio_male_adult_subadult_to_female_adult_juvenile = {
      num <- n_male_adult + n_male_subadult
      den <- n_female_adult + n_female_juvenile
      ifelse(den > 0, num / den, NA_real_)
    }
  ) %>%
  arrange(year, group)

# -------------------- write outputs --------------------
write.csv(summary_tbl,  output_file_monthly, row.names = FALSE)
write.csv(annual_means, output_file_annual,  row.names = FALSE)
message(sprintf("Wrote monthly: %s", output_file_monthly))
message(sprintf("Wrote annual means (filtered): %s", output_file_annual))