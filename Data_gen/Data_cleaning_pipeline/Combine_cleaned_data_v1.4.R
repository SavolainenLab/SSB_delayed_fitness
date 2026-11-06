####### Create analysis datasets #######
## Author: Louis Bell-Roberts
## Description: Executes modular data-wrangling scripts in isolated environments,
## collates all variables, calculates
## future fitness and future rank variables, and exports the final analysis cohort.
# ==============================================================================

library(lubridate)
library(tidyverse)
library(readxl)

# -----------------------
# 0. Repository Paths & Master Execution of Modular Wrangling Scripts
# -----------------------
repo_root    <- "/SSB_delayed_fitness/Data_gen"
supp_xlsx    <- file.path(repo_root, "Data/Raw_inputs/Supplementary_Tables_Collated_v1.4.xlsx")
out_dir      <- file.path(repo_root, "Data/Cleaned_outputs")
pipeline_dir <- file.path(repo_root, "Data_cleaning_pipeline")

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

wrangling_scripts <- c(
  "Scans/Generate_scans_per_year_v1.2.R",
  "Mounts_and_Age/Generate_behav_data_from_master_repeated_measures.R",
  "Rank/Generate_rank_David_v1.2_group_by_observer.R",
  "Fecundity/Calculate_annual_fecundity_v1.3.R",
  "Death_date/Death_dates.R",
  "Dyadic/Generate_dyadic_data.R"
)

cat("==================================================================\n")
cat(" STEP 1: RUNNING MODULAR DATA WRANGLING SCRIPTS\n")
cat("==================================================================\n")

for (rel_path in wrangling_scripts) {
  full_script_path <- file.path(pipeline_dir, rel_path)
  if (!file.exists(full_script_path)) {
    stop(sprintf("Missing pipeline script: %s", full_script_path))
  }
  cat(sprintf("\n--> Executing: %s ...\n", rel_path))
  start_time <- Sys.time()
  source(full_script_path, local = new.env(), echo = FALSE)
  elapsed <- round(as.numeric(difftime(Sys.time(), start_time, units = "secs")), 1)
  cat(sprintf("    Completed in %s seconds.\n", elapsed))
}

cat("\n==================================================================\n")
cat(" STEP 2: COLLATING CLEANED DATASETS & ASSIGNING OBSERVERS\n")
cat("==================================================================\n\n")

# -----------------------
# Read in modular datasets
# -----------------------
pedigree <- as.data.frame(readxl::read_excel(supp_xlsx, sheet = 20, col_types = "text")) %>%
  mutate(
    AnimalID = trimws(as.character(AnimalID)),
    Sex      = toupper(trimws(as.character(Sex)))
  )

mount_data_collapsed      <- read.csv(file.path(out_dir, "output_SSB_data_juv_incl_collapsed.csv"))
mount_data_uncollapsed    <- read.csv(file.path(out_dir, "output_SSB_data_juv_incl.csv"))
scan_data                 <- read.csv(file.path(out_dir, "Scans_per_year.csv")) %>% rename(Animal_ID = Individual_ID)
fecundity_annual_data_pre <- read.csv(file.path(out_dir, "Annual_fecundity.csv"))
rank_data                 <- read.csv(file.path(out_dir, "Rank_per_year.csv")) %>% rename(Animal_ID = ID)

# Demographics from Sheet 20
full_ped_pre <- pedigree
if (!("CSStatus" %in% names(full_ped_pre))) {
  full_ped_pre <- full_ped_pre %>%
    mutate(CSStatus = case_when(
      !is.na(DOD) & as.character(DOD) != ""                   ~ "DEAD",
      !is.na(TransferDate) & as.character(TransferDate) != "" ~ "REMOVE",
      TRUE                                                    ~ "IN CS"
    ))
}

full_ped <- full_ped_pre %>%
  rename(Animal_ID     = AnimalID,
         death_date    = DOD,
         cs_status     = CSStatus,
         transfer_date = TransferDate) %>%
  mutate(
    DOB = as.character(as.Date(suppressWarnings(as.numeric(DOB)), origin = "1899-12-30")),
    DOB = ifelse(is.na(DOB), substr(as.character(full_ped_pre$DOB), 1, 10), DOB),
    death_date = as.character(as.Date(suppressWarnings(as.numeric(death_date)), origin = "1899-12-30")),
    death_date = ifelse(is.na(death_date), substr(as.character(full_ped_pre$DOD), 1, 10), death_date),
    death_date = ifelse(is.na(death_date) | death_date == "", NA_character_, death_date),
    transfer_date = as.character(as.Date(suppressWarnings(as.numeric(transfer_date)), origin = "1899-12-30")),
    transfer_date = ifelse(is.na(transfer_date), substr(as.character(full_ped_pre$TransferDate), 1, 10), transfer_date),
    transfer_date = ifelse(is.na(transfer_date) | transfer_date == "", NA_character_, transfer_date)
  ) %>%
  select(Animal_ID, DOB, cs_status, death_date, transfer_date) %>%
  distinct(Animal_ID, .keep_all = TRUE)

fecundity_annual_data <- fecundity_annual_data_pre %>% 
  select(Animal_ID, Year, num_offspring_all_births)

combine_base <- function(mount_df) {
  mount_df %>%
    full_join(scan_data, by = c("Animal_ID", "Year")) %>%
    full_join(fecundity_annual_data, by = c("Animal_ID", "Year")) %>%
    full_join(rank_data, by = c("Animal_ID", "Year")) %>%
    select(-DOB) %>%
    left_join(full_ped, by = "Animal_ID") %>%
    mutate(num_offspring_all_births = replace_na(num_offspring_all_births, 0))
}

combined_data_collapsed   <- combine_base(mount_data_collapsed)
combined_data_uncollapsed <- combine_base(mount_data_uncollapsed)

# -----------------------
# Group & Observer Assignment
# -----------------------
# 1. Social Group Membership for 2023-2025 (Sheet 23 / Table 23)
sheet23_groups <- as.data.frame(readxl::read_excel(supp_xlsx, sheet = 23)) %>%
  select(Year, Animal_ID = ID, Group_Sheet23 = Group) %>%
  mutate(Animal_ID = trimws(as.character(Animal_ID)), Year = as.integer(Year)) %>%
  distinct()

# 2. Ad-lib records for 2017, 2019, 2020 (Sheet 13)
adlib_df <- as.data.frame(readxl::read_excel(supp_xlsx, sheet = 13, col_types = "text")) %>%
  mutate(Year = year(parse_date_time(as.character(Date), orders = c("ymd", "ymd HMS"), tz = "UTC")))

obs_to_group <- c(
  "JC"="R", "SR"="V", "SE"="V", "AK"="V",
  "CC"="R", "Holly"="R", "KP"="V", "WR"="R", "JF"="V"
)

obs_years <- tibble::tribble(
  ~OBS,   ~Year,
  "JC",   2017,  "JC",   2019,  "JC",   2020,
  "SR",   2017,
  "SE",   2019,  "SE",   2020,
  "CC",   2023,
  "AK",   2023,
  "JF",   2024,  "KP",   2024,
  "WR",   2024,
  "KP",   2025,
  "Holly", 2025
)

obs_assign_table <- obs_years %>%
  mutate(Group = dplyr::recode(OBS, !!!obs_to_group, .default = NA_character_)) %>%
  filter(!is.na(Group)) %>%
  group_by(Year, Group) %>%
  summarise(observers = list(sort(unique(OBS))), .groups = "drop") %>%
  mutate(
    Observer_assigned = vapply(observers, function(x) paste(x, collapse = ";"), FUN.VALUE = character(1)),
    Observer_assigned = ifelse(Year == 2024 & Group == "V", "JF+KP", Observer_assigned)
  ) %>%
  select(Year, Group, Observer_assigned)

id_cols <- names(adlib_df)[str_detect(names(adlib_df), regex("^(actor|recipient|recip)", ignore_case = TRUE))]
ids_long <- adlib_df %>%
  select(Year, Group, all_of(id_cols)) %>%
  pivot_longer(cols = all_of(id_cols), names_to = "role", values_to = "Animal_ID") %>%
  mutate(
    Animal_ID = trimws(as.character(Animal_ID)),
    Group     = trimws(as.character(Group))
  ) %>%
  filter(!is.na(Year), !is.na(Group), Animal_ID != "", !is.na(Animal_ID), !str_detect(Animal_ID, "^[0-9]+$"))

# Annual adlib majority group (with deterministic tie-breaker)
adlib_grp_annual <- ids_long %>%
  group_by(Animal_ID, Year, Group) %>%
  summarise(n = n(), .groups = "drop") %>%
  arrange(Animal_ID, Year, desc(n), Group) %>%
  group_by(Animal_ID, Year) %>%
  slice_head(n = 1) %>%
  select(Animal_ID, Year, Group_adlib_annual = Group)

# Overall lifetime adlib majority group (for 2020 scan-only males)
adlib_grp_overall <- ids_long %>%
  group_by(Animal_ID, Group) %>%
  summarise(n = n(), .groups = "drop") %>%
  arrange(Animal_ID, desc(n), Group) %>%
  group_by(Animal_ID) %>%
  slice_head(n = 1) %>%
  select(Animal_ID, Group_adlib_overall = Group)

assign_group_and_observer <- function(df_input) {
  df_input %>%
    mutate(Animal_ID = as.character(Animal_ID), Year = as.integer(Year)) %>%
    left_join(sheet23_groups, by = c("Animal_ID", "Year")) %>%
    left_join(adlib_grp_annual, by = c("Animal_ID", "Year")) %>%
    left_join(adlib_grp_overall, by = "Animal_ID") %>%
    mutate(
      assigned_group = dplyr::case_when(
        # 2023-2025: Strictly defined by Sheet 23 (excludes split-group movers MLTH, QEYC, QPIN)
        Year >= 2023                 ~ Group_Sheet23,
        # 2017, 2019, 2020: Defined by ad-lib majority, rank group, or scan-only fallback
        !is.na(Group_adlib_annual)  ~ Group_adlib_annual,
        !is.na(Group_for_DS)        ~ Group_for_DS,
        !is.na(Group_adlib_overall) ~ Group_adlib_overall,
        TRUE                        ~ NA_character_
      )
    ) %>%
    left_join(obs_assign_table, by = c("Year", "assigned_group" = "Group")) %>%
    rename(Observer = Observer_assigned) %>%
    select(-Group_Sheet23, -Group_adlib_annual, -Group_adlib_overall)
}

combined_with_group_collapsed   <- assign_group_and_observer(combined_data_collapsed)
combined_with_group_uncollapsed <- assign_group_and_observer(combined_data_uncollapsed)

# -----------------------
# Future Offspring and Future Rank Calculation Pipeline
# -----------------------
process_future_and_filter <- function(combined_df) {
  
  unique_individuals <- unique(combined_df$Animal_ID)
  skeleton_df <- crossing(Animal_ID = unique_individuals, Year = 2017:2025)
  
  data_wrangle <- combined_df %>%
    full_join(skeleton_df, by = c("Animal_ID", "Year")) %>%
    mutate(num_offspring_all_births = replace_na(num_offspring_all_births, 0)) %>%
    arrange(Animal_ID, Year)
  
  # Population exit year calculation
  data_wrangle <- data_wrangle %>%
    mutate(
      death_year    = as.numeric(substr(death_date, 1, 4)),
      transfer_year = as.numeric(substr(transfer_date, 1, 4)),
      exit_year     = pmin(death_year, transfer_year, na.rm = TRUE)
    )
  
  # 1. Compute Future Offspring (up to 6 years ahead)
  for (i in 1:6) {
    col_name        <- paste0("offspring_plus_", i)
    change_col_name <- paste0("offspring_change_plus_", i)
    
    future_offspring <- data_wrangle %>%
      select(Animal_ID, Year, num_offspring_all_births) %>%
      mutate(Year = Year - i) %>%
      rename(!!col_name := num_offspring_all_births)
    
    data_wrangle <- data_wrangle %>%
      left_join(future_offspring, by = c("Animal_ID", "Year")) %>%
      mutate(
        target_year = Year + i,
        is_eligible = is.na(exit_year) | exit_year >= target_year,
        !!col_name := if_else(target_year < 2024 & is_eligible, .data[[col_name]], NA_integer_),
        !!change_col_name := .data[[col_name]] - num_offspring_all_births
      ) %>%
      select(-target_year, -is_eligible)
  }
  
  # 2. Filter to Focal Male Analysis Cohort & Compute Scaled Variables
  analysis_df_pre <- data_wrangle %>%
    filter(
      !is.na(assigned_group),
      Scan_Count >= 10
    ) %>%
    mutate(
      Sex = "M",
      sex = "M"
    ) %>%
    group_by(Animal_ID) %>% 
    mutate(
      age_B = mean(Age, na.rm = TRUE),
      age_W = Age - age_B
    ) %>%
    ungroup() %>%
    mutate(
      Age_s                = as.numeric(scale(Age)),
      age_B_s              = as.numeric(scale(age_B)),
      age_W_s              = as.numeric(scale(age_W)),
      Group_Year           = as.factor(paste(assigned_group, Year, sep = "_")),
      adjusted_SSB         = (SSB_Raw) / log(Scan_Count),
      adjusted_DSB         = (DSB_Raw) / log(Scan_Count),
      adjusted_SSB_linear  = SSB_Raw / Scan_Count,
      adjusted_DSB_linear  = DSB_Raw / Scan_Count,
      death_year_extracted = year(parse_date_time(death_date, orders = c("ymd"), tz = "UTC")),
      died_current_year    = ifelse(!is.na(death_year_extracted) & death_year_extracted == Year, 1, 0),
      died_next_year       = ifelse(!is.na(death_year_extracted) & death_year_extracted == (Year + 1), 1, 0),
      died_current_or_next = if_else(!is.na(death_year_extracted) & (death_year_extracted == Year | death_year_extracted == Year + 1), 1, 0)
    ) %>% 
    group_by(Year, Group_for_DS) %>%
    mutate(DS_sc = as.numeric(scale(DS))) %>%
    ungroup()
  
  # 3. Compute Future Rank Variables (up to 8 years ahead)
  valid_ds_years <- c(2017, 2019, 2020, 2023, 2024, 2025)
  analysis_df    <- analysis_df_pre
  
  for (i in 1:8) {
    new_col_name    <- paste0("Rank_Future_plus_", i, "_s")
    change_col_name <- paste0("Rank_Change_plus_", i, "_s")
    
    future_data <- analysis_df %>%
      filter(Year %in% valid_ds_years) %>%
      select(Animal_ID, Year, DS_sc) %>%
      mutate(Year = Year - i) %>%
      rename(!!new_col_name := DS_sc)
    
    analysis_df <- left_join(analysis_df, future_data, by = c("Animal_ID", "Year")) %>%
      mutate(
        !!change_col_name := .data[[new_col_name]] - DS_sc
      )
  }
  
  return(analysis_df)
}

cat("==================================================================\n")
cat(" STEP 3: COMPUTING FUTURE FITNESS & RANK VARIABLES FOR FINAL COHORT\n")
cat("==================================================================\n")

final_collapsed   <- process_future_and_filter(combined_with_group_collapsed)
final_uncollapsed <- process_future_and_filter(combined_with_group_uncollapsed)

# -----------------------
# Write Final Cleaned Datasets to CSV
# -----------------------
collapsed_out   <- file.path(out_dir, "Combined_cleaned_collapsed_w_OBS.csv")
uncollapsed_out <- file.path(out_dir, "Combined_cleaned_uncollapsed_w_OBS.csv")

write.csv(final_collapsed,   file = collapsed_out,   row.names = FALSE)
write.csv(final_uncollapsed, file = uncollapsed_out, row.names = FALSE)

cat(sprintf("\nData collation complete!\nFinal analysis files saved to:\n - %s\n - %s\nRows: %d | Unique Focal Males: %d\n",
            collapsed_out, uncollapsed_out, nrow(final_collapsed), n_distinct(final_collapsed$Animal_ID)))
