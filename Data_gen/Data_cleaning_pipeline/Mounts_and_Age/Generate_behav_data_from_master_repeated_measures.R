####### Create my own behavioural records from the behaviour master database - repeated measures #######
##Louis Bell-Roberts
#09/09/2025

library(dplyr)
library(data.table)
library(tidyr)
library(stringr)
library(lubridate)
library(readxl)
library(purrr)
library(tidyverse)
library(openxlsx)

# --- Repository paths ---
repo_root <- "/SSB_delayed_fitness/Data_gen"
supp_xlsx <- file.path(repo_root, "Data/Raw_inputs/Supplementary_Tables_Collated_v1.4.xlsx")
out_dir   <- file.path(repo_root, "Data/Cleaned_outputs")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

##########
##Pedigree information (Sheet 20)####
# Read without col_types='text' to preserve Date formatting; uppercase Sex
pedigree <- as.data.frame(readxl::read_excel(supp_xlsx, sheet = 20)) %>%
  mutate(
    AnimalID = trimws(as.character(AnimalID)),
    Sex      = toupper(trimws(as.character(Sex)))
  )

#####################
#Read in master behavioural data file (Sheet 13)####
# Restore the 5 placeholder codes that were mapped in Anonymisation_Key_v1.2.csv
placeholder_restore <- c(
  "MNQE" = "FUN",
  "LTFO" = "MUN",
  "VPBE" = "SAMUN",
  "GLXD" = "JVUN",
  "HYCZ" = "UN"
)

master_data_read_pre <- as.data.frame(readxl::read_excel(supp_xlsx, sheet = 13, col_types = "text")) %>%
  mutate(
    Action = trimws(as.character(Action)),
    Action = dplyr::case_match(
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
  )

master_data_pre <- master_data_read_pre %>%
  filter(Action == "M", !is.na(Date)) #Only interested in mounts, so filter out the rest. Only want entries where the date is known

#Parse dates into chosen POSIXct format
master_data_pre <- master_data_pre %>%
  mutate(Date = parse_date_time(as.character(Date), 
                                orders = c("ymd", "ymd HMS"), 
                                tz = "UTC"))
master_data <- master_data_pre

#################
####Identify each individual's date of birth from Anonymised Pedigree (Sheet 20)####
id_dob_df <- pedigree %>%
  select(Animal_ID = AnimalID, DOB) %>%
  filter(!is.na(Animal_ID), Animal_ID != "", !is.na(DOB)) %>%
  distinct(Animal_ID, .keep_all = TRUE) %>%
  mutate(
    DOB = parse_date_time(substr(as.character(DOB), 1, 10), orders = c("ymd"), tz = "UTC"),
    DOB_conflict = FALSE
  )

#################
#Identify mounts that involve juveniles in the master database and potentially exclude them####

juve_exclude_func <- function(extract_juvenile_info) {
  
  # Local copy so running extract_juvenile_info = TRUE does not overwrite global master_data
  master_data_local <- master_data
  
  if (extract_juvenile_info) {
    master_data_local <- master_data_local %>%
      # Add DOB columns for actor and recipient
      left_join(id_dob_df %>% rename(DOB_actor = DOB), by = c("ACTOR1" = "Animal_ID")) %>%
      left_join(id_dob_df %>% rename(DOB_recipient = DOB), by = c("RECIPIENT" = "Animal_ID")) %>%
      
      # Add sex columns for actor and recipient
      left_join(pedigree %>% select(AnimalID, Sex) %>% rename(Sex_actor = Sex), by = c("ACTOR1" = "AnimalID")) %>%
      left_join(pedigree %>% select(AnimalID, Sex) %>% rename(Sex_recipient = Sex), by = c("RECIPIENT" = "AnimalID")) %>%
      
      # Calculate age, juvenile flags, and filter
      mutate(
        Age_actor = interval(DOB_actor, Date) / years(1),
        Age_recipient = interval(DOB_recipient, Date) / years(1),
        
        Actor_is_juv = case_when(
          Sex_actor == "M" & Age_actor < 3.5 ~ TRUE,
          Sex_actor == "F" & Age_actor < 3 ~ TRUE,
          TRUE ~ FALSE
        ),
        Recipient_is_juv = case_when(
          Sex_recipient == "M" & Age_recipient < 3.5 ~ TRUE,
          Sex_recipient == "F" & Age_recipient < 3 ~ TRUE,
          TRUE ~ FALSE
        ),
        
        Involves_Juvenile = Actor_is_juv | Recipient_is_juv
      ) %>%
      filter(Involves_Juvenile == FALSE | is.na(Involves_Juvenile))
  }
  
  #################
  #Identify males and females from the pedigree####
  
  # Filter pedigree to only include males (excluding 1963 'HYCZ'/'UN' artifact)
  pedigree_males <- pedigree %>%
    filter(Sex == "M", AnimalID != "HYCZ")
  all_males <- unique(pedigree_males$AnimalID)
  
  #Identify known females from the pedigree
  pedigree_females <- pedigree %>%
    filter(Sex == "F")
  
  #############
  ####Generate SSB data (Your exact original counting logic)####
  
  generate_ssb_summary <- function(df, all_males, pedigree_females, extract_juvenile_info = TRUE) {
    
    # Extract year from the Date column
    df <- df %>%
      mutate(Year = year(Date))
    
    # Create a base dataframe of all known males x all observed years
    base <- expand.grid(
      Animal_ID = all_males,
      Year = unique(df$Year),
      stringsAsFactors = FALSE
    )
    
    if (extract_juvenile_info) {
      # DSB counts: male -> female
      dsb <- df %>%
        filter(
          ACTOR1 %in% all_males,
          RECIPIENT %in% pedigree_females$AnimalID | RECIPIENT == "FUN"
        ) %>%
        count(Animal_ID = ACTOR1, Year, name = "DSB_Raw")
      
      # SSB giver counts: male -> male
      ssb_giver <- df %>%
        filter(
          ACTOR1 %in% all_males,
          RECIPIENT %in% all_males | RECIPIENT %in% c("SAMUN", "MUN")
        ) %>%
        count(Animal_ID = ACTOR1, Year, name = "SSB_Giver")
      
      # SSB receiver counts: male <- male
      ssb_receiver <- df %>%
        filter(
          ACTOR1 %in% c(all_males, "SAMUN", "MUN"),
          RECIPIENT %in% all_males
        ) %>%
        count(Animal_ID = RECIPIENT, Year, name = "SSB_Receiver")
      
    } else {
      # DSB counts (include juveniles)
      dsb <- df %>%
        filter(
          ACTOR1 %in% all_males,
          RECIPIENT %in% c(pedigree_females$AnimalID, "FUN", "JVF", "JVFUN")
        ) %>%
        count(Animal_ID = ACTOR1, Year, name = "DSB_Raw")
      
      # SSB giver counts (include juveniles)
      ssb_giver <- df %>%
        filter(
          ACTOR1 %in% all_males,
          RECIPIENT %in% c(all_males, "SAMUN", "MUN", "JVM", "JVMUN")
        ) %>%
        count(Animal_ID = ACTOR1, Year, name = "SSB_Giver")
      
      # SSB receiver counts (include juveniles)
      ssb_receiver <- df %>%
        filter(
          ACTOR1 %in% c(all_males, "SAMUN", "MUN", "JVM", "JVMUN"),
          RECIPIENT %in% all_males
        ) %>%
        count(Animal_ID = RECIPIENT, Year, name = "SSB_Receiver")
    }
    
    # Merge all counts with the base list of males x years
    summary_counts_df <- base %>%
      left_join(ssb_giver, by = c("Animal_ID", "Year")) %>%
      left_join(ssb_receiver, by = c("Animal_ID", "Year")) %>%
      left_join(dsb, by = c("Animal_ID", "Year")) %>%
      mutate(
        SSB_Giver = replace_na(SSB_Giver, 0),
        SSB_Receiver = replace_na(SSB_Receiver, 0),
        DSB_Raw = replace_na(DSB_Raw, 0),
        SSB_Raw = SSB_Giver + SSB_Receiver,
        SSB_proportion = if_else((SSB_Raw + DSB_Raw) > 0,
                                 SSB_Raw / (SSB_Raw + DSB_Raw),
                                 NA_real_)
      )
    
    return(summary_counts_df)
  }
  
  #Create summary counts
  summary_mounts <- generate_ssb_summary(df = master_data_local, all_males = all_males, pedigree_females = pedigree_females, extract_juvenile_info = extract_juvenile_info)
  
  #################
  #Adjusting for multiple mounting behaviour####
  
  mounts_df <- setDT(copy(master_data_local))
  mounts_df <- mounts_df[!is.na(Time)]
  
  mounts_df[, Time_clean := ifelse(grepl(" ", as.character(Time)), sub(".* ", "", as.character(Time)), as.character(Time))]
  mounts_df[, DateTime := parse_date_time(paste(as.character(as.Date(Date)), Time_clean), 
                                          orders = c("ymd HMS", "ymd HM"), 
                                          tz = "UTC")]
  
  mounts_df[, dyad := paste(ACTOR1, RECIPIENT, sep = "_")]
  setorder(mounts_df, dyad, DateTime)
  
  mounts_df[, time_diff := as.numeric(difftime(DateTime, shift(DateTime), units = "mins")), by = dyad]
  mounts_df[, event_id := cumsum(is.na(time_diff) | time_diff > 5), by = dyad]
  
  collapsed_mounts <- mounts_df[, .SD[1], by = .(dyad, Date, event_id)]
  
  #Create summary counts for collapsed mount dataset
  summary_collapsed_mounts <- generate_ssb_summary(df = collapsed_mounts, all_males = all_males, pedigree_females = pedigree_females, extract_juvenile_info = extract_juvenile_info)
  
  return(list(
    output_SSB_data = summary_mounts,
    output_SSB_collapsed_data = summary_collapsed_mounts
  ))
}

###############
output_SSB_data_juve_exclude <- juve_exclude_func(extract_juvenile_info = TRUE)
output_SSB_data_juve_include <- juve_exclude_func(extract_juvenile_info = FALSE)

#Calculate each individual ID's age at the beginning of each study period
start_date_df <- master_data_read_pre %>%
  filter(!is.na(Date)) %>%
  mutate(Date = parse_date_time(as.character(Date), 
                                orders = c("ymd", "ymd HMS"), 
                                tz = "UTC"))

earliest_dates <- start_date_df %>%
  mutate(Year = year(Date)) %>%
  group_by(Year) %>%
  summarise(start_date = min(Date, na.rm = TRUE), .groups = "drop")

age_by_year <- earliest_dates %>%
  mutate(dummy = 1) %>%
  inner_join(id_dob_df %>% mutate(dummy = 1), by = "dummy", relationship = "many-to-many") %>%
  select(-dummy) %>%
  mutate(
    Age = interval(DOB, start_date) / years(1)
  ) %>%
  filter(Age >= 0)

#Join to the four mount count dataframes
output_SSB_data_juve_exclude$output_SSB_data <- output_SSB_data_juve_exclude$output_SSB_data %>%
  left_join(age_by_year %>% select(Animal_ID, Year, Age, DOB), by = c("Animal_ID", "Year"))

output_SSB_data_juve_exclude$output_SSB_collapsed_data <- output_SSB_data_juve_exclude$output_SSB_collapsed_data %>%
  left_join(age_by_year %>% select(Animal_ID, Year, Age, DOB), by = c("Animal_ID", "Year"))

output_SSB_data_juve_include$output_SSB_data <- output_SSB_data_juve_include$output_SSB_data %>%
  left_join(age_by_year %>% select(Animal_ID, Year, Age, DOB), by = c("Animal_ID", "Year"))

output_SSB_data_juve_include$output_SSB_collapsed_data <- output_SSB_data_juve_include$output_SSB_collapsed_data %>%
  left_join(age_by_year %>% select(Animal_ID, Year, Age, DOB), by = c("Animal_ID", "Year"))

#Assign correct variable structure
output_SSB_data_juve_exclude$output_SSB_data$DOB <- as.character(output_SSB_data_juve_exclude$output_SSB_data$DOB)
output_SSB_data_juve_exclude$output_SSB_collapsed_data$DOB <- as.character(output_SSB_data_juve_exclude$output_SSB_collapsed_data$DOB)
output_SSB_data_juve_include$output_SSB_data$DOB <- as.character(output_SSB_data_juve_include$output_SSB_data$DOB)
output_SSB_data_juve_include$output_SSB_collapsed_data$DOB <- as.character(output_SSB_data_juve_include$output_SSB_collapsed_data$DOB)

#Write these dataframes to file
write.csv(x = output_SSB_data_juve_exclude$output_SSB_data, file = file.path(out_dir, "output_SSB_data_juv_excl.csv"), row.names = FALSE)
write.csv(x = output_SSB_data_juve_exclude$output_SSB_collapsed_data, file = file.path(out_dir, "output_SSB_data_juv_excl_collapsed.csv"), row.names = FALSE)
write.csv(x = output_SSB_data_juve_include$output_SSB_data, file = file.path(out_dir, "output_SSB_data_juv_incl.csv"), row.names = FALSE)
write.csv(x = output_SSB_data_juve_include$output_SSB_collapsed_data, file = file.path(out_dir, "output_SSB_data_juv_incl_collapsed.csv"), row.names = FALSE)
