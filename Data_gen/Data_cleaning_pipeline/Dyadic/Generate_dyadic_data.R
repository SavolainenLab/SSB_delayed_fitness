####### Generate Dyadic Behavioral & Attribute Dataset (2023-2025) #######
## Author: Louis Bell-Roberts
## Description: Self-contained modular script to generate the within-group dyadic
## interaction and attribute dataset for 2023-2025, matching Supplementary Table 18.
# ==============================================================================

library(lubridate)
library(tidyverse)
library(readxl)

# --- 1. Repository Paths ---
repo_root <- "/SSB_delayed_fitness/Data_gen"
supp_xlsx <- file.path(repo_root, "Data/Raw_inputs/Supplementary_Tables_Collated_v1.4.xlsx")
out_dir   <- file.path(repo_root, "Data/Cleaned_outputs")

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# --- 2. Ingest Pedigree, Sheet 23 Groups, Rank, and Age Data ---
pedigree <- as.data.frame(readxl::read_excel(supp_xlsx, sheet = 20, col_types = "text")) %>%
  mutate(
    AnimalID = trimws(as.character(AnimalID)),
    Sex      = toupper(trimws(as.character(Sex)))
  )

sheet23_groups <- as.data.frame(readxl::read_excel(supp_xlsx, sheet = 23)) %>%
  dplyr::select(Year, Animal_ID = ID, Group_Sheet23 = Group) %>%
  mutate(Animal_ID = trimws(as.character(Animal_ID)), Year = as.integer(Year)) %>%
  distinct()

rank_data <- read.csv(file.path(out_dir, "Rank_per_year.csv")) %>%
  rename(Animal_ID = ID)

mount_data_collapsed <- read.csv(file.path(out_dir, "output_SSB_data_juv_incl_collapsed.csv"))

# --- 3. Define Male & Female ID Vectors ---
non_study_males <- c("EWLM", "HYCZ", "UFOB", "XIXZ", "OYLB", "NEQG")
male_IDs        <- pedigree %>% filter(Sex == "M", !(AnimalID %in% non_study_males)) %>% pull(AnimalID) %>% unique()
female_IDs      <- pedigree %>% filter(Sex == "F") %>% pull(AnimalID) %>% unique()

placeholder_restore <- c(
  "MNQE" = "FUN",
  "LTFO" = "MUN",
  "VPBE" = "SAMUN",
  "GLXD" = "JVUN",
  "HYCZ" = "UN"
)

exclude_names   <- c("FUN", "JVF", "JVFUN", "JVM", "JVUN", "MUN", "SAMUN", "JVMUN", "INFANT", "JUV", "UN")
exclude_pattern <- "FUN|JVUN|MUN|SAMUN|JVMUN|INFANT|JUV|JVF|JVFUN|JVM|^UN$"

# --- 4. Load & Standardise Ad-Lib Data (Sheet 13) for 2023-2025 ---
master_data_pre <- as.data.frame(readxl::read_excel(supp_xlsx, sheet = 13, col_types = "text")) %>%
  mutate(
    Date     = substr(as.character(Date), 1, 10),
    Time_str = trimws(as.character(Time)),
    Time_str = ifelse(grepl(" ", Time_str), sub(".* ", "", Time_str), Time_str),
    Year     = as.integer(substr(Date, 1, 4)),
    Group    = trimws(as.character(Group)),
    Action   = trimws(as.character(Action)),
    Action   = dplyr::case_match(
      Action,
      "mounting"            ~ "M",
      "allo-grooming"       ~ "Gr",
      "aggression"          ~ "Ag",
      "submission"          ~ "Su",
      "sexual solicitation" ~ "SSo",
      "genital stimulation" ~ "GS",
      "genital inspection"  ~ "GI",
      "penile erection"     ~ "PE",
      "masturbation"        ~ "Ma",
      "copulatory plug"     ~ "CP",
      "self-grooming"       ~ "SGr",
      .default = Action
    ),
    across(
      c(ACTOR1, ACTOR2, Actor3, RECIPIENT, Recip2, Recip3),
      ~ dplyr::recode(trimws(as.character(.)), !!!placeholder_restore, .default = trimws(as.character(.)))
    )
  ) %>%
  filter(Year %in% c(2023, 2024, 2025))

# --- 5. Load Scan Data (Sheet 15) & Infer Scan Group Context ---
scan_data_pre <- as.data.frame(readxl::read_excel(supp_xlsx, sheet = 15, col_types = "text")) %>%
  mutate(
    DATE          = substr(as.character(DATE), 1, 10),
    TIME          = trimws(as.character(TIME)),
    TIME          = ifelse(grepl(" ", TIME), sub(".* ", "", TIME), TIME),
    Year          = as.integer(substr(DATE, 1, 4)),
    Individual_ID = trimws(as.character(Individual_ID)),
    Individual_ID = dplyr::recode(Individual_ID, !!!placeholder_restore, .default = Individual_ID)
  ) %>%
  filter(
    Year %in% c(2023, 2024, 2025),
    Individual_ID %in% c(male_IDs, female_IDs)
  )

indiv_observability <- scan_data_pre %>%
  group_by(Year, Individual_ID) %>%
  summarise(Total_Indiv_Scans = n(), .groups = "drop")

valid_monkeys_list <- indiv_observability %>%
  filter(Total_Indiv_Scans >= 10, Individual_ID %in% male_IDs) %>%
  dplyr::select(Year, Individual_ID)

scan_event_groups <- scan_data_pre %>%
  left_join(sheet23_groups, by = c("Year", "Individual_ID" = "Animal_ID")) %>%
  group_by(Year, DATE, TIME) %>%
  mutate(
    n_R = sum(Group_Sheet23 == "R", na.rm = TRUE),
    n_V = sum(Group_Sheet23 == "V", na.rm = TRUE),
    Scan_Group = case_when(
      !is.na(Group_Sheet23) ~ Group_Sheet23,
      n_R >= n_V & n_R > 0  ~ "R",
      n_V > n_R             ~ "V",
      TRUE                  ~ "R"
    )
  ) %>%
  ungroup()

# --- 6. Calculate Yearly Co-Occurrence (Shared Scans) ---
scan_data_filtered <- scan_event_groups %>%
  filter(Individual_ID %in% male_IDs, !(Individual_ID %in% c(exclude_names, female_IDs))) %>%
  dplyr::select(Year, Scan_Group, DATE, TIME, Individual_ID) %>%
  distinct() %>%
  mutate(present = 1L) %>%
  pivot_wider(
    names_from  = Individual_ID,
    values_from = present,
    values_fill = 0L
  )

pair_counts_yearly <- scan_data_filtered %>%
  group_by(Year) %>%
  nest() %>%
  mutate(
    pair_data = map(data, function(df) {
      id_matrix <- df %>%
        dplyr::select(-Scan_Group, -DATE, -TIME) %>%
        as.matrix()
      
      co_occurrence_matrix <- t(id_matrix) %*% id_matrix
      
      as.data.frame(as.table(co_occurrence_matrix), stringsAsFactors = FALSE) %>%
        rename(ID1 = Var1, ID2 = Var2, Shared_Scans = Freq) %>%
        filter(ID1 != ID2, Shared_Scans > 0) %>%
        mutate(
          ID1      = as.character(ID1),
          ID2      = as.character(ID2),
          true_ID1 = pmin(ID1, ID2),
          true_ID2 = pmax(ID1, ID2),
          ID1      = true_ID1,
          ID2      = true_ID2,
          dyad     = paste(ID1, ID2, sep = "-")
        ) %>%
        distinct(dyad, .keep_all = TRUE) %>%
        dplyr::select(ID1, ID2, dyad, Shared_Scans)
    })
  ) %>%
  dplyr::select(Year, pair_data) %>%
  unnest(cols = c(pair_data)) %>%
  ungroup()

# --- 7. Dyadic Grooming, Mounting (5-min Collapsed), Aggression & Submission ---
actor_cols     <- c("ACTOR1", "ACTOR2", "Actor3")
recipient_cols <- c("RECIPIENT", "Recip2", "Recip3")

extract_dyadic_action <- function(df, actions_vec) {
  long_events <- df %>%
    filter(Action %in% actions_vec, Year %in% c(2023, 2024, 2025)) %>%
    mutate(
      obs_id   = row_number(),
      Datetime = parse_date_time(paste(Date, Time_str), orders = c("ymd HMS", "ymd HM", "ymd"), tz = "UTC")
    ) %>%
    dplyr::select(obs_id, Year, Date, Datetime, Action, all_of(actor_cols), all_of(recipient_cols)) %>%
    pivot_longer(cols = all_of(actor_cols), names_to = "actor_role", values_to = "actor", values_drop_na = TRUE) %>%
    pivot_longer(cols = all_of(recipient_cols), names_to = "recip_role", values_to = "recip", values_drop_na = TRUE) %>%
    filter(
      actor != recip,
      actor != "", recip != "",
      actor != "NA", recip != "NA",
      !str_detect(actor, exclude_pattern),
      !str_detect(recip, exclude_pattern),
      actor %in% male_IDs,
      recip %in% male_IDs
    )
  
  # 5-minute (300s) bout-collapsing rule for mounting
  m_collapsed <- long_events %>%
    filter(Action == "M") %>%
    arrange(Date, actor, recip, Datetime) %>%
    group_by(Date, actor, recip) %>%
    mutate(
      diff_sec = as.numeric(difftime(Datetime, lag(Datetime), units = "secs"))
    ) %>%
    filter(is.na(diff_sec) | diff_sec > 300) %>%
    ungroup() %>%
    dplyr::select(-diff_sec)
  
  other_actions <- long_events %>%
    filter(Action != "M")
  
  bind_rows(m_collapsed, other_actions) %>%
    mutate(
      ID1  = pmin(actor, recip),
      ID2  = pmax(actor, recip),
      dyad = paste(ID1, ID2, sep = "-")
    ) %>%
    group_by(Year, ID1, ID2, dyad, Action) %>%
    summarise(
      Total_Count = n(),
      ID1_to_ID2  = sum(actor == ID1),
      ID2_to_ID1  = sum(actor == ID2),
      .groups = "drop"
    )
}

dyads_all_actions <- extract_dyadic_action(master_data_pre, c("Gr", "M", "Ag", "Su")) %>%
  pivot_wider(
    names_from  = Action,
    values_from = c(Total_Count, ID1_to_ID2, ID2_to_ID1),
    names_glue  = "{.value}_{Action}",
    values_fill = 0L
  )

# --- 8. Coalition Counts (Co-Aggressors & Co-Victims) ---
coalitions_raw <- master_data_pre %>%
  filter(Action == "Ag", Year %in% c(2023, 2024, 2025)) %>%
  mutate(
    obs_id = row_number(),
    across(
      all_of(c(actor_cols, recipient_cols)),
      ~ ifelse(. %in% female_IDs | is.na(.) | . == "NA", "", as.character(.))
    )
  )

actors_long <- coalitions_raw %>%
  dplyr::select(obs_id, Year, Date, Action, all_of(actor_cols)) %>%
  pivot_longer(cols = all_of(actor_cols), names_to = "role", values_to = "ID", values_drop_na = TRUE) %>%
  filter(ID != "", !str_detect(ID, exclude_pattern), ID %in% male_IDs)

actor_pairs <- actors_long %>%
  inner_join(actors_long, by = c("obs_id", "Year", "Date", "Action"), relationship = "many-to-many") %>%
  filter(ID.x < ID.y) %>%
  rename(ID1 = ID.x, ID2 = ID.y) %>%
  mutate(dyad = paste(ID1, ID2, sep = "-")) %>%
  dplyr::select(Year, Date, Action, ID1, ID2, dyad)

recips_long <- coalitions_raw %>%
  dplyr::select(obs_id, Year, Date, Action, all_of(recipient_cols)) %>%
  pivot_longer(cols = all_of(recipient_cols), names_to = "role", values_to = "ID", values_drop_na = TRUE) %>%
  filter(ID != "", !str_detect(ID, exclude_pattern), ID %in% male_IDs)

recip_pairs <- recips_long %>%
  inner_join(recips_long, by = c("obs_id", "Year", "Date", "Action"), relationship = "many-to-many") %>%
  filter(ID.x < ID.y) %>%
  rename(ID1 = ID.x, ID2 = ID.y) %>%
  mutate(dyad = paste(ID1, ID2, sep = "-")) %>%
  dplyr::select(Year, Date, Action, ID1, ID2, dyad)

yearly_counts_coalitions <- bind_rows(actor_pairs, recip_pairs) %>%
  group_by(Year, ID1, ID2, dyad) %>%
  summarise(Coalition_count = n(), .groups = "drop")

# --- 9. Construct Complete Qualifying Group Dyad Grid (5,944 Dyads) ---
mover_group_presence <- tibble::tribble(
  ~Year, ~Animal_ID, ~Group,
  2023L, "MLTH",     "R",
  2023L, "MLTH",     "V",
  2023L, "QPIN",     "R",
  2023L, "QPIN",     "V",
  2024L, "QEYC",     "R",
  2024L, "QEYC",     "V"
)

qualifying_group_males <- sheet23_groups %>%
  rename(Group = Group_Sheet23) %>%
  bind_rows(mover_group_presence) %>%
  inner_join(valid_monkeys_list, by = c("Year", "Animal_ID" = "Individual_ID")) %>%
  distinct(Year, Group, Animal_ID)

dyad_skeleton <- qualifying_group_males %>%
  group_by(Year, Group) %>%
  summarise(
    dyads = list({
      u <- sort(unique(Animal_ID))
      if (length(u) >= 2) {
        as.data.frame(t(combn(u, 2)), stringsAsFactors = FALSE) %>%
          rename(ID1 = V1, ID2 = V2)
      } else {
        tibble(ID1 = character(), ID2 = character())
      }
    }),
    .groups = "drop"
  ) %>%
  unnest(cols = c(dyads)) %>%
  mutate(dyad = paste(ID1, ID2, sep = "-")) %>%
  arrange(Year, Group, dyad) %>%
  distinct(Year, dyad, .keep_all = TRUE)

# --- 10. Merge Dyadic Counts & Individual Attributes (Pedigree DOB & Age) ---
rank_prep <- rank_data %>%
  dplyr::select(Year, ID = Animal_ID, DS, Group_for_DS) %>%
  distinct()

age_prep <- mount_data_collapsed %>%
  dplyr::select(Year, ID = Animal_ID, Age, DOB) %>%
  mutate(DOB = substr(as.character(DOB), 1, 10)) %>%
  distinct()

indiv_attributes <- full_join(rank_prep, age_prep, by = c("Year", "ID"))

expected_behav_cols <- c(
  "Total_Count_Gr", "Total_Count_M", "ID1_to_ID2_Gr", "ID1_to_ID2_M", "ID2_to_ID1_Gr", "ID2_to_ID1_M",
  "Total_Count_Ag", "ID1_to_ID2_Ag", "ID2_to_ID1_Ag",
  "Total_Count_Su", "ID1_to_ID2_Su", "ID2_to_ID1_Su"
)

dyadic_counts_final <- dyad_skeleton %>%
  left_join(indiv_observability %>% rename(ID1 = Individual_ID, ID1_yearly_scans = Total_Indiv_Scans), by = c("Year", "ID1")) %>%
  left_join(indiv_observability %>% rename(ID2 = Individual_ID, ID2_yearly_scans = Total_Indiv_Scans), by = c("Year", "ID2")) %>%
  left_join(pair_counts_yearly, by = c("Year", "ID1", "ID2", "dyad")) %>%
  left_join(dyads_all_actions, by = c("Year", "ID1", "ID2", "dyad")) %>%
  left_join(yearly_counts_coalitions, by = c("Year", "ID1", "ID2", "dyad"))

for (col_nm in expected_behav_cols) {
  if (!(col_nm %in% names(dyadic_counts_final))) {
    dyadic_counts_final[[col_nm]] <- 0L
  }
}

dyadic_counts_final <- dyadic_counts_final %>%
  mutate(
    across(
      all_of(c("ID1_yearly_scans", "ID2_yearly_scans", "Shared_Scans", expected_behav_cols, "Coalition_count")),
      ~ replace_na(as.integer(.), 0L)
    )
  ) %>%
  left_join(indiv_attributes, by = c("Year", "ID1" = "ID")) %>%
  rename(
    DS_ID1           = DS,
    Group_for_DS_ID1 = Group_for_DS,
    Age_ID1          = Age,
    DOB_ID1          = DOB
  ) %>%
  left_join(indiv_attributes, by = c("Year", "ID2" = "ID")) %>%
  rename(
    DS_ID2           = DS,
    Group_for_DS_ID2 = Group_for_DS,
    Age_ID2          = Age,
    DOB_ID2          = DOB
  ) %>%
  mutate(
    Age_Diff_ID1_minus_ID2  = Age_ID1 - Age_ID2,
    Age_Diff_ID2_minus_ID1  = Age_ID2 - Age_ID1,
    Rank_Diff_ID1_minus_ID2 = DS_ID1 - DS_ID2,
    Rank_Diff_ID2_minus_ID1 = DS_ID2 - DS_ID1
  ) %>%
  dplyr::select(
    Year, Group, ID1, ID2, dyad,
    ID1_yearly_scans, ID2_yearly_scans, Shared_Scans,
    Total_Count_Gr, Total_Count_M, ID1_to_ID2_Gr, ID1_to_ID2_M, ID2_to_ID1_Gr, ID2_to_ID1_M,
    Total_Count_Ag, ID1_to_ID2_Ag, ID2_to_ID1_Ag,
    Total_Count_Su, ID1_to_ID2_Su, ID2_to_ID1_Su,
    Coalition_count,
    DS_ID1, Group_for_DS_ID1, Age_ID1, DOB_ID1,
    DS_ID2, Group_for_DS_ID2, Age_ID2, DOB_ID2,
    Age_Diff_ID1_minus_ID2, Age_Diff_ID2_minus_ID1,
    Rank_Diff_ID1_minus_ID2, Rank_Diff_ID2_minus_ID1
  ) %>%
  arrange(Year, Group, dyad)

dyadic_out <- file.path(out_dir, "Dyadic_counts.csv")
write.csv(dyadic_counts_final, file = dyadic_out, row.names = FALSE)

cat(sprintf("\nDyadic dataset generation complete!\nSaved to: %s\nRows: %d | Columns: %d\n",
            dyadic_out, nrow(dyadic_counts_final), ncol(dyadic_counts_final)))
