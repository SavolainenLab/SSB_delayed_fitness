#### Calculate David scores ####
## Louis Bell-Roberts
## 17/09/2025

# Load required packages
library(tidyverse)
library(EloRating)
library(lubridate)
library(data.table)
library(janitor)
library(purrr)
library(readxl)

# --- Repository paths ---
repo_root <- "/SSB_delayed_fitness/Data_gen"
supp_xlsx <- file.path(repo_root, "Data/Raw_inputs/Supplementary_Tables_Collated_v1.4.xlsx")
out_dir   <- file.path(repo_root, "Data/Cleaned_outputs")

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

#### Read in data and standardise ####
# Load anonymised pedigree (Sheet 20)
pedigree <- as.data.frame(readxl::read_excel(supp_xlsx, sheet = 20)) %>%
  mutate(
    AnimalID = trimws(as.character(AnimalID)),
    Sex      = toupper(trimws(as.character(Sex)))
  )

# non_study_males <- c("EWLM", "HYCZ", "UFOB", "XIXZ", "OYLB", "NEQG")
non_study_males <- c()

male_ids <- pedigree %>%
  filter(Sex == "M", !(AnimalID %in% non_study_males)) %>%
  pull(AnimalID)

# Placeholder restoration vector to match original raw ad-lib coding
placeholder_restore <- c(
  "MNQE" = "FUN",
  "LTFO" = "MUN",
  "VPBE" = "SAMUN",
  "GLXD" = "JVUN",
  "HYCZ" = "UN"
)

# Load anonymised master behavioural data file (Sheet 13)
master_data_read_pre <- as.data.frame(readxl::read_excel(supp_xlsx, sheet = 13, col_types = "text")) %>%
  mutate(
    Action = trimws(as.character(Action)),
    Action = dplyr::case_match(
      Action,
      "aggression" ~ "Ag",
      "submission" ~ "Su",
      .default = Action
    ),
    Response = trimws(as.character(Response)),
    across(
      c(ACTOR1, ACTOR2, Actor3, RECIPIENT, Recip2, Recip3),
      ~ dplyr::recode(trimws(as.character(.)), !!!placeholder_restore, .default = trimws(as.character(.)))
    )
  )

# Parse Date into POSIXct (UTC, format: ymd)
master_data_pre <- master_data_read_pre %>%
  mutate(Date = parse_date_time(as.character(Date), orders = c("ymd", "ymd HMS"), tz = "UTC"))

cols_to_fix <- c("ACTOR1", "ACTOR2", "Actor3", "RECIPIENT", "Recip2", "Recip3")

# Exclude draws (properly handling 'ignore', 'IgnoreRefusal', 'IgnoreRejection')
master_data_pre <- master_data_pre %>%
  filter(
    !is.na(Date) &
      Action %in% c("Ag", "Su") & 
      (is.na(Response) | !(Response %in% c("IgnoreRefusal", "IgnoreRejection", "ignore")))
  ) %>%
  mutate(across(all_of(cols_to_fix), ~ na_if(as.character(.), "")))

# Filter participant columns strictly to recognized study males
master_data_pre <- master_data_pre %>%
  mutate(across(
    all_of(cols_to_fix),
    ~ ifelse(.x %in% male_ids, .x, NA_character_)
  ))

# Define actor and recipient columns
actor_cols <- c("ACTOR1", "ACTOR2", "Actor3")
recip_cols <- c("RECIPIENT", "Recip2", "Recip3")

# Keep only rows with at least one actor AND at least one recipient
master_data_all_years <- master_data_pre %>%
  filter(
    if_any(all_of(actor_cols), ~ !is.na(.)) &
      if_any(all_of(recip_cols), ~ !is.na(.))
  )

# Split dataset by year
master_data_all_years <- master_data_all_years %>% mutate(Year = year(Date))

df_by_year <- master_data_all_years %>% group_split(Year)
names(df_by_year) <- master_data_all_years %>% distinct(Year) %>% pull(Year)

# Function to filter interactions based on majority group for multi-group individuals
filter_by_majority_group <- function(df) {
  
  interaction_counts <- df %>%
    pivot_longer(
      cols = all_of(cols_to_fix),
      names_to = "Role",
      values_to = "AnimalID"
    ) %>%
    filter(!is.na(AnimalID)) %>%
    group_by(AnimalID, Group) %>%
    summarise(n_interactions = n(), .groups = "drop")
  
  # Deterministic tie-breaking by sorting Group alphabetically
  majority_group <- interaction_counts %>%
    arrange(AnimalID, desc(n_interactions), Group) %>%
    group_by(AnimalID) %>%
    slice_head(n = 1) %>%
    dplyr::select(AnimalID, Group_majority = Group)
  
  df_filtered <- df %>%
    rowwise() %>%
    filter(
      all(
        c(ACTOR1, ACTOR2, Actor3, RECIPIENT, Recip2, Recip3) %>% 
          na.omit() %>% 
          map_lgl(~ {
            maj_group <- majority_group$Group_majority[majority_group$AnimalID == .x]
            is.na(maj_group) || Group == maj_group
          })
      )
    ) %>%
    ungroup()
  
  return(df_filtered)
}

df_by_year_filtered <- map(df_by_year, filter_by_majority_group)

# Extract dyadic interaction from each interaction type
extract_interactions <- function(data, action_type, actor_cols, recip_cols) {
  
  temp <- data %>% filter(Action == action_type)
  interactions <- data.frame()
  
  if (nrow(temp) == 0) return(interactions)
  
  for (row_idx in 1:nrow(temp)) {
    actors <- na.omit(as.character(temp[row_idx, actor_cols]))
    recips <- na.omit(as.character(temp[row_idx, recip_cols]))
    date   <- temp$Date[row_idx]
    
    # Handle polyadic contest where MLGE was secondary co-actor with HSDF against WVTQ
    if ("MLGE" %in% actors && "HSDF" %in% actors && "WVTQ" %in% recips) {
      actors <- setdiff(actors, "MLGE")
    }
    
    if (length(actors) > 0 && length(recips) > 0) {
      row_interactions <- expand.grid(winner = actors, loser = recips, date = date, stringsAsFactors = FALSE)
      interactions <- bind_rows(interactions, row_interactions)
    }
  }
  
  return(interactions)
}

# Create function to repeat pipeline for each year
by_year_rank_function <- function(data_per_year) {
  
  ag_interactions <- extract_interactions(data = data_per_year, action_type = "Ag", actor_cols = c("ACTOR1", "ACTOR2", "Actor3"), recip_cols = c("RECIPIENT", "Recip2", "Recip3"))
  su_interactions <- extract_interactions(data = data_per_year, action_type = "Su", actor_cols = c("RECIPIENT", "Recip2", "Recip3"), recip_cols = c("ACTOR1", "ACTOR2", "Actor3"))
  
  all_interactions <- bind_rows(ag_interactions, su_interactions) %>%
    filter(!is.na(winner) & winner != "NA",
           !is.na(loser)  & loser  != "NA",
           !is.na(date))
  
  all_ids <- sort(unique(c(all_interactions$winner, all_interactions$loser)))
  
  win_matrix <- matrix(0, nrow = length(all_ids), ncol = length(all_ids),
                       dimnames = list(all_ids, all_ids))
  
  interaction_table <- table(all_interactions$winner, all_interactions$loser)
  win_matrix[rownames(interaction_table), colnames(interaction_table)] <- interaction_table
  
  diag(win_matrix) <- 0
  
  ds <- DS(win_matrix)
  return(ds)
}

# Split by Year x Group, then run DS per subset
df_by_year_group <- map(df_by_year_filtered, ~ .x %>%
                          group_by(Year, Group) %>%
                          group_split())
df_by_year_group <- flatten(df_by_year_group)

by_year_group_rank_function <- function(data_per_group) {
  year_val  <- unique(data_per_group$Year)
  group_val <- unique(data_per_group$Group)
  
  ds <- by_year_rank_function(data_per_group)
  
  tibble(
    Year         = year_val,
    Group_for_DS = group_val,
    ID           = ds$ID,
    DS           = ds$DS,
    DS_norm      = ds$normDS
  )
}

ds_results <- map_dfr(df_by_year_group, by_year_group_rank_function)

# Write output to Data/Cleaned_outputs/
write.csv(ds_results, file.path(out_dir, "Rank_per_year.csv"), row.names = FALSE)
cat("Rank data successfully generated and saved to:", file.path(out_dir, "Rank_per_year.csv"), "\n")
