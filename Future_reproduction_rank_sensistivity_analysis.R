#### Examine the future reproductive and rank consequences of performing SSB - sensitivity analysis to individual differences in observability
##Louis Bell-Roberts
#10/02/2026

library(tidyr)
library(dplyr)
library(tibble)
library(brms)
library(ggplot2)
library(coda)
library(tidybayes)
library(purrr)
library(stringr)
library(posterior)
library(future)
library(parallel)
library(bayestestR)
library(lubridate)

# Set base directory for data loading
base_dir <- "/Volumes/Macintosh HD/Users/pn26026/Documents/Sandisk/Imperial_HPC/Future_models/"

# 1. SET UP PARALLEL PROCESSING ACROSS DATASETS
plan(multisession, workers = 10)

# 2. Data Import ---------------------------------------------------------------
# File path (Local working path; convert to relative path prior to final archive)
supp_data_path <- "/Users/pn26026/Library/CloudStorage/OneDrive-UniversityofBristol/SSB/Data/Supplementary_Tables_Collated_v1.3.xlsx"

# Sheet 18: Combined cleaned collapsed data with observation effort
# Sheet 19: Cleaned pedigree data
data <- readxl::read_excel(path = supp_data_path, sheet = 18)
pedigree_pre <- readxl::read_excel(path = supp_data_path, sheet = 19)

# 3. Data Wrangling ------------------------------------------------------------
data <- data %>% 
  arrange(Animal_ID, Year)

# 1. Determine population exit year for each individual
# We need to know the final year they were eligible to reproduce.
data <- data %>%
  mutate(
    # Extract the 4-digit year from the death_date and transfer_date strings
    death_year = as.numeric(substr(death_date, 1, 4)),
    transfer_year = as.numeric(substr(transfer_date, 1, 4)),
    
    # The exit_year is the earliest of death_year or transfer_year.
    # 'pmin' compares row by row. 'na.rm = TRUE' ensures that if an animal 
    # only has a death date, it doesn't return NA due to a missing transfer date.
    exit_year = pmin(death_year, transfer_year, na.rm = TRUE)
  )

#Filter data and perform transformations
analysis_df <- data %>%
  filter(
    sex == "M",
    !is.na(assigned_group),
    Scan_Count >= 10
  ) %>%
  group_by(Animal_ID) %>% 
  mutate( #Calculate age variables
    age_B = mean(Age),
    age_W = Age - age_B
  ) %>%
  ungroup() %>%
  mutate( #Scale age variables
    Animal_ID = as.factor(Animal_ID),
    Age_s = as.numeric(scale(Age)),
    age_B_s = as.numeric(scale(age_B)),
    age_W_s = as.numeric(scale(age_W)),
    Group_Year = as.factor(paste(assigned_group, Year, sep = "_")),
    adjusted_SSB = (SSB_Raw)/log(Scan_Count),
    adjusted_DSB = (DSB_Raw)/log(Scan_Count),
    adjusted_SSB_linear = SSB_Raw/Scan_Count,
    adjusted_DSB_linear = DSB_Raw/Scan_Count
  ) %>%
  mutate( 
    # Extract death year for logic checks
    death_year_extracted = year(death_date),
    
    # --- Death Flags (for reference) ---
    died_current_year = ifelse(!is.na(death_year_extracted) & death_year_extracted == Year, 1, 0),
    died_next_year    = ifelse(!is.na(death_year_extracted) & death_year_extracted == (Year + 1), 1, 0),
    died_current_or_next = if_else(!is.na(death_year_extracted) & (death_year_extracted == Year | death_year_extracted == Year + 1), 1, 0),
    
    assigned_group = as.factor(assigned_group)
  ) %>% 
  group_by(Year, Group_for_DS) %>% #For scaling DS, individuals are grouped by 'Group_for_DS' as this is the group where their relative rank was calculated within.
  mutate(
    DS_sc = scale(DS)
  ) %>%
  ungroup()

#Calculate estimates for SSB and DSB####
# bivariate model specification
priors_SSB_DSB <- c(
  # --- FIXED EFFECTS ---
  prior(normal(0, 5),   class = Intercept, resp = "SSBRaw"),
  prior(normal(0, 5),   class = Intercept, resp = "DSBRaw"),
  prior(normal(0, 0.5), class = b,         resp = "SSBRaw"),
  prior(normal(0, 0.5), class = b,         resp = "DSBRaw"),
  
  # --- RANDOM EFFECTS ---
  # 1. SD for the INTERCEPTS
  prior(exponential(1), class = sd, coef = "Intercept", group = "Animal_ID", resp = "SSBRaw"),
  prior(exponential(1), class = sd, coef = "Intercept", group = "Animal_ID", resp = "DSBRaw"),
  
  # 2. SD for the SLOPES
  prior(lkj(2), class = cor, group = "Animal_ID"),
  
  # assigned_group:Year (Intercepts only)
  prior(exponential(1), class = sd, group = "assigned_group:Year", resp = "SSBRaw"),
  prior(exponential(1), class = sd, group = "assigned_group:Year", resp = "DSBRaw"),

  # --- DISTRIBUTION PARAMETERS ---
  prior(exponential(1), class = shape, resp = "SSBRaw"),
  prior(exponential(1), class = shape, resp = "DSBRaw")
)

#Model formulas
bf_SSB <- bf(SSB_Raw ~ 1 + age_W_s + age_B_s + offset(log(Scan_Count)) + (1|p|Animal_ID) + (1|assigned_group:Year)) + negbinomial()

bf_DSB <- bf(DSB_Raw ~ 1 + age_W_s + age_B_s + offset(log(Scan_Count)) + (1|p|Animal_ID) + (1|assigned_group:Year)) + negbinomial()

biv_SSB_DSB <- brm(
  bf_SSB + bf_DSB + set_rescor(F),
  prior = priors_SSB_DSB,
  data = analysis_df, seed = 1,
  warmup = 1500, iter = 4000, thin = 1, chains = 4, cores = 4,
  control = list(adapt_delta = 0.99),
  save_pars = save_pars(all = TRUE)
)

#Save model object
# saveRDS(biv_SSB_DSB, file = "/Users/louis.bell-roberts/Documents/SSB_Macaques/Rscripts_local/SEM/SSB_DSB_biv_mod")
# biv_SSB_DSB <- readRDS(file = "/Users/louis.bell-roberts/Documents/SSB_Macaques/Rscripts_local/SEM/SSB_DSB_biv_mod")

summary(biv_SSB_DSB)

# pp_check(biv_SSB_DSB, resp='SSBRaw', type = "bars") + xlim(c(-1, 60))
# pp_check(biv_SSB_DSB, resp='DSBRaw', type = "bars") + xlim(c(-1, 50))

#Extract predictions

# 1. Create the prediction grid
# We filter for the specific years first, then ensure we have the correct 
# Age/Group predictors associated with those specific years.
newdata <- analysis_df %>%
  filter(Year %in% c(2017, 2019, 2020, 2023, 2024, 2025)) %>% #Only calculating values for years in which reproduction data available
  select(Animal_ID, assigned_group, Year, age_W_s, age_B_s) %>%
  distinct() %>%
  # HANDLING THE OFFSET:
  # Using the mean Scan_Count standardises the prediction. 
  # This gives the expected count as if every animal had the "average" scan depth.
  mutate(Scan_Count = mean(analysis_df$Scan_Count))

# 2. Extract Predictions
predictions <- newdata %>%
  add_epred_draws(
    biv_SSB_DSB,
    allow_new_levels = FALSE, #do not make predictions for years that animals weren't present
    ndraws = 1500, # Adjust based on memory/needs
    re_formula = NULL # Includes Animal_ID and assigned_group:Year specific deviations
  ) %>%
  # The output will typically have:
  # .category: "SSBRaw" or "DSBRaw"
  # .epred: The predicted value on the response scale (counts)
  ungroup()

#Calculate z scores
# 1. Prepare the "Standardised" Data
me_ready_data <- predictions %>%
  
  # --- STEP A: Calculate "Local" Z-scores (Per Draw & Per Group) ---
  # We group by 'assigned_group' to force the mean to be specific to that observer.
  group_by(.draw, Year, assigned_group, .category) %>% 
  mutate(
    # The Baseline: "What is normal for THIS observer in THIS year?"
    group_mu = mean(.epred), 
    group_sd = sd(.epred),
    
    # The Z-score: "How does this animal compare to that baseline?"
    # Note: If group_sd is 0 (unlikely), this returns NaN.
    z_raw = (.epred - group_mu) / group_sd
  ) %>%
  
  # --- STEP B: Summarise the Posterior for brms ---
  # Now we collapse the 1000 draws into one Estimate and one Error per animal.
  ungroup() %>%
  group_by(Animal_ID, Year, assigned_group, .category) %>%
  summarise(
    # The Estimate (Mean of the Z-scores) -> The value used for prediction
    z_est = mean(z_raw),
    
    # The Measurement Error (SD of the Z-scores) -> The uncertainty for me()
    z_err = sd(z_raw),
    
    #Calculate raw values from the untransformed posterior distributions
    raw_mean = mean(.epred),
    raw_sd = sd(.epred),
    
    .groups = "drop"
  ) %>%
  
  # --- STEP C: Pivot Wider for the Bivariate Structure ---
  # SSBRaw_z_est, SSBRaw_z_err, SSBRaw_raw_mean, SSBRaw_raw_sd
  # DSBRaw_z_est, DSBRaw_z_err, DSBRaw_raw_mean, DSBRaw_raw_sd
  pivot_wider(
    names_from = .category,
    values_from = c(z_est, z_err, raw_mean, raw_sd),
    names_glue = "{.category}_{.value}"
  )


#Join the predicted values for SSB and DSB to the main analysis dataset
analysis_df_final <- analysis_df %>% left_join(
  me_ready_data, 
  by = c("Animal_ID", "Year", "assigned_group")
)

# MULTIPLE IMPUTATION APPROACH FOR FUTURE REPRODUCTION (LAGS 1 TO 6)

# PREPARE POSTERIOR DRAWS FOR MULTIPLE IMPUTATION

# 1. Calculate Local Z-scores but KEEP the draws intact
posterior_z_draws <- predictions %>%
  group_by(.draw, Year, assigned_group, .category) %>% 
  mutate(
    group_mu = mean(.epred), 
    group_sd = sd(.epred),
    z_raw = (.epred - group_mu) / group_sd
  ) %>%
  ungroup()

# 2. Pivot wider so SSB and DSB are columns for each draw/animal/year combo
# This creates a dataframe with 1500 versions (draws) of every observation
imputed_draws_wide <- posterior_z_draws %>%
  select(.draw, Animal_ID, Year, assigned_group, .category, z_raw) %>%
  pivot_wider(
    names_from = .category,
    values_from = z_raw,
    names_prefix = "imputed_" # This creates columns: imputed_SSBRaw, imputed_DSBRaw
  )

# Verify the maximum number of draws available
max_draws <- max(imputed_draws_wide$.draw)

# ==============================================================================
# GLOBAL SETUP: POSTERIOR MULTIPLE IMPUTATION (100 DRAWS)
# ==============================================================================

# Define the global number of imputations
num_imputations <- 100

# Select 100 specific draw indices globally. 
# We do this once here so that the exact same 100 correlated draws are used 
# across Concurrent models, Future Reproduction models, and Future Rank models.
set.seed(2026) 
selected_draws <- sample(1:max_draws, num_imputations, replace = FALSE)


# ==============================================================================
# PART 1: CONCURRENT REPRODUCTION & RANK
# ==============================================================================
cat("\n--- Preparing Concurrent Data ---\n")

# 1. SETUP: PREPARE CONCURRENT DATA
# Filter out rows missing either the reproductive outcome or the rank outcome
df_concurrent_base <- analysis_df_final %>%
  filter(!is.na(num_offspring_all_births), !is.na(DS_sc))

# 2. GENERATE IMPUTED DATASETS FOR CONCURRENT MODELS
imputed_datasets_conc <- list()

for (j in 1:num_imputations) {
  draw_id <- selected_draws[j]
  
  # Isolate the specific posterior draw for all animals
  draw_data <- imputed_draws_wide %>% filter(.draw == draw_id)
  
  # Join the true posterior draw to our concurrent base dataframe
  temp_df <- df_concurrent_base %>%
    left_join(draw_data, by = c("Animal_ID", "Year", "assigned_group")) %>%
    mutate(
      # Create the squared terms for the quadratic analysis
      imputed_SSBRaw_sq = imputed_SSBRaw^2,
      imputed_DSBRaw_sq = imputed_DSBRaw^2
    )
  
  imputed_datasets_conc[[j]] <- temp_df
}

# --- CONCURRENT REPRODUCTION MODELS ---
cat("\nRunning Multiple Imputation Models: Concurrent Reproduction\n")

# Formulas
bf_repro_conc_linear <- bf(
  num_offspring_all_births ~ 1 + 
    imputed_SSBRaw + imputed_DSBRaw + 
    Age + I(Age^2) + 
    (1 | Animal_ID) + (1 | assigned_group:Year),
  family = negbinomial()
)

bf_repro_conc_quad <- bf(
  num_offspring_all_births ~ 1 + 
    imputed_SSBRaw + imputed_SSBRaw_sq + imputed_DSBRaw + 
    Age + I(Age^2) + 
    (1 | Animal_ID) + (1 | assigned_group:Year),
  family = negbinomial()
)

# Priors
priors_repro_conc <- c(
  prior(normal(0, 2), class = Intercept),
  prior(normal(0, 1), class = b),
  prior(exponential(1), class = sd),
  prior(exponential(1), class = shape)
)

#Open sink
con <- file(file.path(base_dir, "Results/Models/model_output_summaries_test.txt"), open = "wt")
sink(con)
sink(con, type = "message")

# Run Linear Model
cat("Fitting Concurrent Reproduction (Linear)...\n")
fit_concurrent_repro_linear <- brm_multiple(
  formula = bf_repro_conc_linear,
  data    = imputed_datasets_conc,
  prior   = priors_repro_conc,
  seed    = 1,
  chains  = 1, 
  cores   = 1,
  warmup  = 1000, 
  iter    = 1500,
  control = list(adapt_delta = 0.99, max_treedepth = 12),
  save_pars = save_pars(all = TRUE)
)
summary(fit_concurrent_repro_linear)
# Save Model
saveRDS(fit_concurrent_repro_linear, file = file.path(base_dir, "Results/Models/Concurrent_repro_linear.rds"))
# Extract all fixed and random effect coefficients, format as dataframe, and save as CSV
coefs_repro_linear <- as.data.frame(posterior_summary(fit_concurrent_repro_linear)) %>%
  rownames_to_column(var = "Parameter")
write.csv(coefs_repro_linear, file = file.path(base_dir, "Results/Models/Concurrent_repro_linear_coefs.csv"), row.names = F)
# Remove model and table from memory, then force garbage collection to free RAM
rm(fit_concurrent_repro_linear, coefs_repro_linear)
gc()

# Run Quadratic Model
cat("Fitting Concurrent Reproduction (Quadratic)...\n")
fit_concurrent_repro_quad <- brm_multiple(
  formula = bf_repro_conc_quad,
  data    = imputed_datasets_conc,
  prior   = priors_repro_conc,
  seed    = 1,
  chains  = 1, 
  cores   = 1,
  warmup  = 1000, 
  iter    = 1500,
  control = list(adapt_delta = 0.99, max_treedepth = 12),
  save_pars = save_pars(all = TRUE)
)
summary(fit_concurrent_repro_quad)

# Save Model
saveRDS(fit_concurrent_repro_quad, file = file.path(base_dir, "Results/Models/Concurrent_repro_quad.rds"))
# Extract all fixed and random effect coefficients, format as dataframe, and save as CSV
coefs_repro_quad <- as.data.frame(posterior_summary(fit_concurrent_repro_quad)) %>%
  rownames_to_column(var = "Parameter")
write.csv(coefs_repro_quad, file = file.path(base_dir, "Results/Models/Concurrent_repro_quad_coefs.csv"), row.names = F)

# Clear memory and garbage collect
rm(fit_concurrent_repro_quad, coefs_repro_quad)
gc()

# --- CONCURRENT RANK MODELS ---
cat("\nRunning Multiple Imputation Models: Concurrent Rank\n")

# Formulas
bf_rank_conc_linear <- bf(
  DS_sc ~ 1 + 
    imputed_SSBRaw + imputed_DSBRaw + 
    Age + I(Age^2) +
    (1 | Animal_ID) + (1 | assigned_group:Year),
  family = gaussian()
)

bf_rank_conc_quad <- bf(
  DS_sc ~ 1 + 
    imputed_SSBRaw + imputed_SSBRaw_sq + imputed_DSBRaw + 
    Age + I(Age^2) +
    (1 | Animal_ID) + (1 | assigned_group:Year),
  family = gaussian()
)

# Priors
priors_rank_conc <- c(
  prior(normal(0, 1), class = Intercept),
  prior(normal(0, 1), class = b),
  prior(exponential(2), class = sd),
  prior(student_t(3, 0, 1), class = sigma)
)

# Run Linear Model
cat("Fitting Concurrent Rank (Linear)...\n")
fit_concurrent_rank_linear <- brm_multiple(
  formula = bf_rank_conc_linear,
  data    = imputed_datasets_conc,
  prior   = priors_rank_conc,
  seed    = 1,
  chains  = 1, 
  cores   = 1,
  warmup  = 1000, 
  iter    = 1500,
  control = list(adapt_delta = 0.99, max_treedepth = 12),
  save_pars = save_pars(all = TRUE)
)
summary(fit_concurrent_rank_linear)

# Save Model
saveRDS(fit_concurrent_rank_linear, file = file.path(base_dir, "Results/Models/Concurrent_rank_linear.rds"))
# Extract all fixed and random effect coefficients, format as dataframe, and save as CSV
coefs_rank_linear <- as.data.frame(posterior_summary(fit_concurrent_rank_linear)) %>%
  rownames_to_column(var = "Parameter")
write.csv(coefs_rank_linear, file = file.path(base_dir, "Results/Models/Concurrent_rank_linear_coefs.csv"), row.names = F)

# Clear memory and garbage collect
rm(fit_concurrent_rank_linear, coefs_rank_linear)
gc()

# Run Quadratic Model
cat("Fitting Concurrent Rank (Quadratic)...\n")
fit_concurrent_rank_quad <- brm_multiple(
  formula = bf_rank_conc_quad,
  data    = imputed_datasets_conc,
  prior   = priors_rank_conc,
  seed    = 1,
  chains  = 1, 
  cores   = 1,
  warmup  = 1000, 
  iter    = 1500,
  control = list(adapt_delta = 0.99, max_treedepth = 12),
  save_pars = save_pars(all = TRUE)
)
summary(fit_concurrent_rank_quad)

# Save Model
saveRDS(fit_concurrent_rank_quad, file = file.path(base_dir, "Results/Models/Concurrent_rank_quad.rds"))
# Extract all fixed and random effect coefficients, format as dataframe, and save as CSV
coefs_rank_quad <- as.data.frame(posterior_summary(fit_concurrent_rank_quad)) %>%
  rownames_to_column(var = "Parameter")
write.csv(coefs_rank_quad, file = file.path(base_dir, "Results/Models/Concurrent_rank_quad_coefs.csv"), row.names = F)

# Clear memory and garbage collect
rm(fit_concurrent_rank_quad, coefs_rank_quad)
gc()

# ==============================================================================
# PART 2: FUTURE REPRODUCTION (LAGS 1 TO 6)
# ==============================================================================

# Define Dataframes for each time step
repro_one_df   <- analysis_df_final %>% filter(Year %in% c(2017, 2019, 2020))
repro_two_df   <- analysis_df_final %>% filter(Year %in% c(2017, 2019, 2020))
repro_three_df <- analysis_df_final %>% filter(Year %in% c(2017, 2019, 2020))
repro_four_df  <- analysis_df_final %>% filter(Year %in% c(2017, 2019))
repro_five_df  <- analysis_df_final %>% filter(Year %in% c(2017))
repro_six_df   <- analysis_df_final %>% filter(Year %in% c(2017))

# Put them into a list
df_list_repro <- list(
  repro_one_df, repro_two_df, repro_three_df, 
  repro_four_df, repro_five_df, repro_six_df
)

# Initialize output lists
model_results_mult_linear <- list()
model_results_mult_quad <- list()

for(i in 1:6) {
  
  cat(paste0("\n--- Preparing data and models for Repro lag: +", i, " years ---\n"))
  
  current_df_base <- df_list_repro[[i]]
  col_name_offspring <- paste0("offspring_plus_", i)
  
  # Filter out missing outcomes
  current_df_base <- current_df_base %>% filter(!is.na(!!sym(col_name_offspring)))
  
  # Generate imputed datasets
  imputed_datasets <- list()
  for (j in 1:num_imputations) {
    draw_id <- selected_draws[j]
    draw_data <- imputed_draws_wide %>% filter(.draw == draw_id)
    
    temp_df <- current_df_base %>%
      left_join(draw_data, by = c("Animal_ID", "Year", "assigned_group")) %>%
      mutate(
        imputed_SSBRaw_sq = imputed_SSBRaw^2,
        imputed_DSBRaw_sq = imputed_DSBRaw^2
      )
    imputed_datasets[[j]] <- temp_df
  }
  
  # Formulas
  bf_mult_linear <- bf(
    as.formula(paste0(
      col_name_offspring, " ~ 1 + ",
      "imputed_SSBRaw + imputed_DSBRaw + Age + I(Age^2) + (1 | Animal_ID) + (1 | assigned_group:Year)"
    )),
    family = negbinomial()
  )
  
  bf_mult_quad <- bf(
    as.formula(paste0(
      col_name_offspring, " ~ 1 + ",
      "imputed_SSBRaw + imputed_SSBRaw_sq + imputed_DSBRaw + Age + I(Age^2) + (1 | Animal_ID) + (1 | assigned_group:Year)"
    )),
    family = negbinomial()
  )
  
  # Run Linear
  cat(paste0("Fitting Linear Repro Model for lag +", i, "...\n"))
  fit_linear <- brm_multiple(
    formula = bf_mult_linear,
    data = imputed_datasets, 
    prior = priors_repro_conc,
    seed = 1, chains = 1, cores = 1,
    warmup = 1000, iter = 1500, # Maintained warmup/iter here
    control = list(adapt_delta = 0.99, max_treedepth = 12),
    save_pars = save_pars(all = TRUE)
  )
  summary(fit_linear)
  model_results_mult_linear[[paste0("plus_", i)]] <- fit_linear
  
  # Run Quadratic
  cat(paste0("Fitting Quadratic Repro Model for lag +", i, "...\n"))
  fit_quad <- brm_multiple(
    formula = bf_mult_quad,
    data = imputed_datasets,
    prior = priors_repro_conc,
    seed = 1, chains = 1, cores = 1,
    warmup = 1000, iter = 1500, # Maintained warmup/iter here
    control = list(adapt_delta = 0.99, max_treedepth = 12),
    save_pars = save_pars(all = TRUE)
  )
  summary(fit_quad)
  model_results_mult_quad[[paste0("plus_", i)]] <- fit_quad
}

# Save Future Repro Models
# 1. Extract and combine coefficients for all lags in the linear list
coefs_repro_linear_all <- imap_dfr(model_results_mult_linear, ~{
  as.data.frame(posterior_summary(.x)) %>%
    rownames_to_column(var = "Parameter") %>%
    mutate(Lag = .y) # Adds a column identifying the lag (e.g., "plus_1")
})
write.csv(coefs_repro_linear_all, file = file.path(base_dir, "Results/Models/Future_repro_linear_all_coefs.csv"), row.names = F)

# 2. Extract and combine coefficients for all lags in the quadratic list
coefs_repro_quad_all <- imap_dfr(model_results_mult_quad, ~{
  as.data.frame(posterior_summary(.x)) %>%
    rownames_to_column(var = "Parameter") %>%
    mutate(Lag = .y)
})
write.csv(coefs_repro_quad_all, file = file.path(base_dir, "Results/Models/Future_repro_quad_all_coefs.csv"), row.names = F)

# 3. Save the full Model lists
saveRDS(object = model_results_mult_linear, file = file.path(base_dir, "Results/Models/Future_repro_linear_models.rds"))
saveRDS(object = model_results_mult_quad, file = file.path(base_dir, "Results/Models/Future_repro_quad_models.rds"))

# ==============================================================================
# PART 3: FUTURE RANK (LAGS 1 TO 8)
# ==============================================================================

# Define Dataframes
rank_one_df   <- analysis_df_final %>% filter(Year %in% c(2019, 2023, 2024))
rank_two_df   <- analysis_df_final %>% filter(Year %in% c(2017, 2023))
rank_three_df <- analysis_df_final %>% filter(Year %in% c(2017, 2020))
rank_four_df  <- analysis_df_final %>% filter(Year %in% c(2019, 2020))
rank_five_df  <- analysis_df_final %>% filter(Year %in% c(2019, 2020))
rank_six_df   <- analysis_df_final %>% filter(Year %in% c(2017, 2019))
rank_seven_df <- analysis_df_final %>% filter(Year %in% c(2017))
rank_eight_df <- analysis_df_final %>% filter(Year %in% c(2017))

df_list_rank <- list(
  rank_one_df, rank_two_df, rank_three_df, 
  rank_four_df, rank_five_df, rank_six_df, 
  rank_seven_df, rank_eight_df
)

# Initialize output lists
model_results_rank_mult_linear <- list()
model_results_rank_mult_quad <- list()

for(i in 1:8) { 
  
  cat(paste0("\n--- Preparing data and models for Rank lag: +", i, " years ---\n"))
  
  current_df_base <- df_list_rank[[i]]
  col_name_rank_scaled <- paste0("Rank_Future_plus_", i, "_s")
  
  # Filter missing
  current_df_base <- current_df_base %>% filter(!is.na(!!sym(col_name_rank_scaled)))
  
  # Generate imputed datasets
  imputed_datasets <- list()
  for (j in 1:num_imputations) {
    draw_id <- selected_draws[j]
    draw_data <- imputed_draws_wide %>% filter(.draw == draw_id)
    
    temp_df <- current_df_base %>%
      left_join(draw_data, by = c("Animal_ID", "Year", "assigned_group")) %>%
      mutate(
        imputed_SSBRaw_sq = imputed_SSBRaw^2,
        imputed_DSBRaw_sq = imputed_DSBRaw^2
      )
    imputed_datasets[[j]] <- temp_df
  }
  
  # Formulas
  bf_mult_linear <- bf(
    as.formula(paste0(
      col_name_rank_scaled, " ~ 1 + ",
      "imputed_SSBRaw + imputed_DSBRaw + DS_sc + Age + I(Age^2) + (1 | Animal_ID) + (1 | assigned_group:Year)"
    )),
    family = gaussian()
  )
  
  bf_mult_quad <- bf(
    as.formula(paste0(
      col_name_rank_scaled, " ~ 1 + ",
      "imputed_SSBRaw + imputed_SSBRaw_sq + imputed_DSBRaw + DS_sc + Age + I(Age^2) + (1 | Animal_ID) + (1 | assigned_group:Year)"
    )),
    family = gaussian()
  )
  
  # Run Linear
  cat(paste0("Fitting Linear Rank Model for lag +", i, "...\n"))
  fit_linear <- brm_multiple(
    formula = bf_mult_linear,
    data = imputed_datasets, 
    prior = priors_rank_conc,
    seed = 1, chains = 1, cores = 1,
    warmup = 1000, iter = 1500,
    control = list(adapt_delta = 0.99, max_treedepth = 12),
    save_pars = save_pars(all = TRUE)
  )
  summary(fit_linear)
  model_results_rank_mult_linear[[paste0("plus_", i)]] <- fit_linear
  
  # Run Quadratic
  cat(paste0("Fitting Quadratic Rank Model for lag +", i, "...\n"))
  fit_quad <- brm_multiple(
    formula = bf_mult_quad,
    data = imputed_datasets, 
    prior = priors_rank_conc,
    seed = 1, chains = 1, cores = 1,
    warmup = 1000, iter = 1500,
    control = list(adapt_delta = 0.99, max_treedepth = 12),
    save_pars = save_pars(all = TRUE)
  )
  summary(fit_quad)
  model_results_rank_mult_quad[[paste0("plus_", i)]] <- fit_quad
}

# Save Future Rank Models
# 1. Extract and combine coefficients for all lags in the rank linear list
coefs_rank_linear_all <- imap_dfr(model_results_rank_mult_linear, ~{
  as.data.frame(posterior_summary(.x)) %>%
    rownames_to_column(var = "Parameter") %>%
    mutate(Lag = .y)
})
write.csv(coefs_rank_linear_all, file = file.path(base_dir, "Results/Models/Future_rank_linear_all_coefs.csv"), row.names = F)

# 2. Extract and combine coefficients for all lags in the rank quadratic list
coefs_rank_quad_all <- imap_dfr(model_results_rank_mult_quad, ~{
  as.data.frame(posterior_summary(.x)) %>%
    rownames_to_column(var = "Parameter") %>%
    mutate(Lag = .y)
})
write.csv(coefs_rank_quad_all, file = file.path(base_dir, "Results/Models/Future_rank_quad_all_coefs.csv"), row.names = F)

# 3. Save the full Model lists
saveRDS(object = model_results_rank_mult_linear, file = file.path(base_dir, "Results/Models/Future_rank_linear_models.rds"))
saveRDS(object = model_results_rank_mult_quad, file = file.path(base_dir, "Results/Models/Future_rank_quad_models.rds"))

# 4. Wipe from RAM
rm(model_results_rank_mult_linear, model_results_rank_mult_quad, coefs_rank_linear_all, coefs_rank_quad_all)
gc()

# Revert outputs back to the console
sink(type = "message")
sink()

# ==============================================================================
# Generate Unified Model Output Tables for Imputed Models (Memory Efficient)
# ==============================================================================
cat("\n================================================================\n")
cat("Generating Unified Parameter and Diagnostic Tables...\n")
cat("================================================================\n")

# --- 1. Define the Dual-Purpose Extraction Function ---
extract_mi_model_data <- function(model_name, file_name, response_type) {
  
  filepath <- file.path(base_dir, "Results/Models", file_name)
  
  if (!file.exists(filepath)) {
    warning(paste("File missing on disk, skipping:", filepath))
    return(NULL)
  }
  
  cat("Processing model object:", model_name, "...\n")
  mdl <- readRDS(filepath)
  
  # ----------------------------------------------------------------------------
  # PART A: PARAMETER EXTRACTION
  # ----------------------------------------------------------------------------
  mdl_summ <- summary(mdl)
  
  # Extract Fixed Effects
  fixed_eff <- data.frame()
  if (!is.null(mdl_summ$fixed)) {
    fixed_eff <- as.data.frame(mdl_summ$fixed) %>%
      tibble::rownames_to_column("Parameter") %>%
      dplyr::mutate(Effect_Type = "Fixed", Group = NA_character_)
  }
  
  # Extract Random Effects
  random_eff <- data.frame()
  if (length(mdl_summ$random) > 0) {
    ran_df_list <- lapply(names(mdl_summ$random), function(grp) {
      as.data.frame(mdl_summ$random[[grp]]) %>%
        tibble::rownames_to_column("Parameter") %>%
        dplyr::mutate(Effect_Type = "Random", Group = grp)
    })
    random_eff <- dplyr::bind_rows(ran_df_list)
  }
  
  # Extract Distributional Parameters (e.g., shape, sigma)
  spec_pars_eff <- data.frame()
  if (!is.null(mdl_summ$spec_pars) && nrow(mdl_summ$spec_pars) > 0) {
    spec_pars_eff <- as.data.frame(mdl_summ$spec_pars) %>%
      tibble::rownames_to_column("Parameter") %>%
      dplyr::mutate(Effect_Type = "Distributional", Group = NA_character_)
  }
  
  # Bind parameters together and format
  param_df <- dplyr::bind_rows(fixed_eff, random_eff, spec_pars_eff) %>%
    dplyr::mutate(
      Model_Name = model_name,
      Response_Variable = response_type,
      
      # Clean Parameter Labels
      Parameter = stringr::str_replace_all(Parameter, "^cor\\((.*)\\)$", "Cor: \\1"),
      Parameter = ifelse(stringr::str_detect(Parameter, "^Cor: "), stringr::str_replace_all(Parameter, ",", " <-> "), Parameter),
      Parameter = stringr::str_replace_all(Parameter, "^sd\\((.*)\\)$", "SD: \\1"), 
      Parameter = stringr::str_replace_all(Parameter, "imputed_SSBRaw_sq", "SSB^2 (imputed)"),
      Parameter = stringr::str_replace_all(Parameter, "imputed_SSBRaw", "SSB (imputed)"),
      Parameter = stringr::str_replace_all(Parameter, "imputed_DSBRaw", "DSB (imputed)"),
      Parameter = stringr::str_replace_all(Parameter, "IAgeE2", "Age^2"),
      Parameter = stringr::str_replace_all(Parameter, "DS_sc", "Current Rank (scaled)"),
      Parameter = stringr::str_replace_all(Parameter, "assigned_group", "Group: "),
      
      # Round estimates
      Estimate    = round(Estimate, 3),
      Est.Error   = round(Est.Error, 3),
      Lower_95_CI = round(`l-95% CI`, 3),
      Upper_95_CI = round(`u-95% CI`, 3)
    ) %>%
    dplyr::select(Model_Name, Response_Variable, Effect_Type, Group, Parameter, Estimate, Est.Error, Lower_95_CI, Upper_95_CI)
  
  # ----------------------------------------------------------------------------
  # PART B: SAMPLER DIAGNOSTICS & DIVERGENCES
  # ----------------------------------------------------------------------------
  # Extract overall parameter stabilities
  rhats <- brms::rhat(mdl)
  bulks <- brms::ess_bulk(mdl)
  tails <- brms::ess_tail(mdl)
  
  # Extract Divergent transitions per imputation chain
  np <- brms::nuts_params(mdl)
  div_df <- np %>%
    dplyr::filter(Parameter == "divergent__") %>%
    dplyr::group_by(Chain) %>%
    dplyr::summarise(Divergences = sum(Value), .groups = "drop")
  
  diag_summary_df <- data.frame(
    Model_Name = model_name,
    Response_Variable = response_type,
    Total_Imputations = length(unique(div_df$Chain)),
    Max_Rhat = round(max(rhats, na.rm = TRUE), 4),
    Mean_Rhat = round(mean(rhats, na.rm = TRUE), 4),
    Min_Bulk_ESS = round(min(bulks, na.rm = TRUE), 0),
    Mean_Bulk_ESS = round(mean(bulks, na.rm = TRUE), 0),
    Min_Tail_ESS = round(min(tails, na.rm = TRUE), 0),
    Total_Divergent_Transitions = sum(div_df$Divergences),
    Max_Divergent_Single_Imputation = max(div_df$Divergences)
  )
  
  return(list(params = param_df, diagnostics = diag_summary_df))
}

# --- 2. Initialize lists ---
all_params <- list()
all_diags  <- list()

# --- 3. Process Standalone Concurrent Models ---
cat("Extracting Concurrent Models...\n")
standalone_files <- list(
  list(file = "Concurrent_repro_linear.rds", name = "Concurrent Repro (Linear)", resp = "Offspring"),
  list(file = "Concurrent_repro_quad.rds", name = "Concurrent Repro (Quadratic)", resp = "Offspring"),
  list(file = "Concurrent_rank_linear.rds", name = "Concurrent Rank (Linear)", resp = "Rank"),
  list(file = "Concurrent_rank_quad.rds", name = "Concurrent Rank (Quadratic)", resp = "Rank")
)

for (sf in standalone_files) {
  ext_data <- extract_mi_model_data(sf$name, sf$file, sf$resp)
  if (!is.null(ext_data)) {
    all_params[[sf$name]] <- ext_data$params
    all_diags[[sf$name]]  <- ext_data$diagnostics
  }
}

# --- 4. Process Future Lag Lists ---
cat("Extracting Future Lag Models...\n")
future_files <- list(
  list(file = "Future_repro_linear_models.rds", prefix = "Future Repro Linear", resp = "Offspring"),
  list(file = "Future_repro_quad_models.rds", prefix = "Future Repro Quad", resp = "Offspring"),
  list(file = "Future_rank_linear_models.rds", prefix = "Future Rank Linear", resp = "Rank"),
  list(file = "Future_rank_quad_models.rds", prefix = "Future Rank Quad", resp = "Rank")
)

for (ff in future_files) {
  path <- file.path(base_dir, "Results/Models", ff$file)
  if (file.exists(path)) {
    mod_list <- readRDS(path)
    
    for (lag_name in names(mod_list)) {
      full_name <- paste0(ff$prefix, " (", lag_name, ")")
      
      # We pass the model object directly here to avoid re-reading from disk
      # (Slight adjustment for list objects)
      cat("Processing model object:", full_name, "...\n")
      mdl <- mod_list[[lag_name]]
      
      # A: Params
      mdl_summ <- summary(mdl)
      fixed_eff <- if(!is.null(mdl_summ$fixed)) as.data.frame(mdl_summ$fixed) %>% tibble::rownames_to_column("Parameter") %>% dplyr::mutate(Effect_Type = "Fixed", Group = NA_character_) else data.frame()
      random_eff <- if(length(mdl_summ$random) > 0) dplyr::bind_rows(lapply(names(mdl_summ$random), function(grp) as.data.frame(mdl_summ$random[[grp]]) %>% tibble::rownames_to_column("Parameter") %>% dplyr::mutate(Effect_Type = "Random", Group = grp))) else data.frame()
      spec_pars_eff <- if(!is.null(mdl_summ$spec_pars) && nrow(mdl_summ$spec_pars) > 0) as.data.frame(mdl_summ$spec_pars) %>% tibble::rownames_to_column("Parameter") %>% dplyr::mutate(Effect_Type = "Distributional", Group = NA_character_) else data.frame()
      
      param_df <- dplyr::bind_rows(fixed_eff, random_eff, spec_pars_eff) %>%
        dplyr::mutate(Model_Name = full_name, Response_Variable = ff$resp, Parameter = stringr::str_replace_all(Parameter, "^cor\\((.*)\\)$", "Cor: \\1"), Parameter = ifelse(stringr::str_detect(Parameter, "^Cor: "), stringr::str_replace_all(Parameter, ",", " <-> "), Parameter), Parameter = stringr::str_replace_all(Parameter, "^sd\\((.*)\\)$", "SD: \\1"), Parameter = stringr::str_replace_all(Parameter, "imputed_SSBRaw_sq", "SSB^2 (imputed)"), Parameter = stringr::str_replace_all(Parameter, "imputed_SSBRaw", "SSB (imputed)"), Parameter = stringr::str_replace_all(Parameter, "imputed_DSBRaw", "DSB (imputed)"), Parameter = stringr::str_replace_all(Parameter, "IAgeE2", "Age^2"), Parameter = stringr::str_replace_all(Parameter, "DS_sc", "Current Rank (scaled)"), Parameter = stringr::str_replace_all(Parameter, "assigned_group", "Group: "), Estimate = round(Estimate, 3), Est.Error = round(Est.Error, 3), Lower_95_CI = round(`l-95% CI`, 3), Upper_95_CI = round(`u-95% CI`, 3)) %>%
        dplyr::select(Model_Name, Response_Variable, Effect_Type, Group, Parameter, Estimate, Est.Error, Lower_95_CI, Upper_95_CI)
      
      # B: Diagnostics
      rhats <- brms::rhat(mdl); bulks <- brms::ess_bulk(mdl); tails <- brms::ess_tail(mdl)
      np <- brms::nuts_params(mdl)
      div_df <- np %>% dplyr::filter(Parameter == "divergent__") %>% dplyr::group_by(Chain) %>% dplyr::summarise(Divergences = sum(Value), .groups = "drop")
      
      diag_summary_df <- data.frame(Model_Name = full_name, Response_Variable = ff$resp, Total_Imputations = length(unique(div_df$Chain)), Max_Rhat = round(max(rhats, na.rm = TRUE), 4), Mean_Rhat = round(mean(rhats, na.rm = TRUE), 4), Min_Bulk_ESS = round(min(bulks, na.rm = TRUE), 0), Mean_Bulk_ESS = round(mean(bulks, na.rm = TRUE), 0), Min_Tail_ESS = round(min(tails, na.rm = TRUE), 0), Total_Divergent_Transitions = sum(div_df$Divergences), Max_Divergent_Single_Imputation = max(div_df$Divergences))
      
      all_params[[full_name]] <- param_df
      all_diags[[full_name]]  <- diag_summary_df
    }
    
    rm(mod_list); gc() 
  }
}

# --- 5. Export Final Tables ---
cat("Exporting final CSVs...\n")

if (length(all_params) > 0) {
  final_params_table <- dplyr::bind_rows(all_params)
  write.csv(final_params_table, file.path(base_dir, "Results/Models/Supp_Table_1_All_Parameters.csv"), row.names = FALSE)
  cat("SUCCESS! Parameters table saved.\n")
  
  final_diags_table <- dplyr::bind_rows(all_diags)
  write.csv(final_diags_table, file.path(base_dir, "Results/Models/Supp_Table_2_MI_Diagnostics.csv"), row.names = FALSE)
  cat("SUCCESS! Diagnostics table saved.\n")
} else {
  cat("WARNING: No models were found to summarize.\n")
}

#Finish
future::plan(sequential)
gc() 