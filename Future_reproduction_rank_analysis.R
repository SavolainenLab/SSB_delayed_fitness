##### Examine the future reproductive and rank consequences of performing SSB
##Louis Bell-Roberts
#10/02/2026

# 1. Required Libraries --------------------------------------------------------
# Data import and manipulation
library(readxl)
library(tidyverse)
library(broom.mixed)

# Bayesian modeling and posterior diagnostics
library(brms)
library(coda)
library(bayesplot)
library(tidybayes)
library(posterior)
library(bayestestR)

# Visualization and causal modeling
library(dagitty)
library(ggdag)
library(ggplot2)
library(ggstance)
library(ggridges)
library(cowplot)
library(patchwork)

# 2. Data Import ---------------------------------------------------------------
# File path (Local working path; convert to relative path prior to final archive)
supp_data_path <- "/SSB_delayed_fitness/Supplementary_Tables_Collated_v1.4.xlsx"

# Sheet 18: Combined cleaned collapsed data with observation effort
# Sheet 19: Cleaned pedigree data
data <- readxl::read_excel(path = supp_data_path, sheet = 19)
pedigree_pre <- readxl::read_excel(path = supp_data_path, sheet = 20)

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

# Analysis of current reproduction and rank, applying linear and quadratic terms to SSB ####
# SETUP: PREPARE CONCURRENT DATA
# 1. Create a specific dataframe for concurrent analysis
df_concurrent <- analysis_df %>%
  # Filter out rows with missing data
  filter(!is.na(adjusted_SSB), !is.na(adjusted_DSB), !is.na(DS_sc))

# 2. Scale Variables for Concurrent Analysis
df_concurrent <- df_concurrent %>%
  mutate(
    adjusted_SSB_z = as.numeric(scale(adjusted_SSB)),
    adjusted_DSB_z = as.numeric(scale(adjusted_DSB))
  )

# MODEL 1: CONCURRENT REPRODUCTION
cat("\nRunning Joint Models: Concurrent Reproduction\n")

# Current Reproduction (QUADRATIC Outcome)
bf_repro_conc_quad <- bf(
  num_offspring_all_births ~ 1 + 
    adjusted_SSB_z + I(adjusted_SSB_z^2) + 
    adjusted_DSB_z +
    Age_s + I(Age_s^2) + 
    (1 | Animal_ID) + (1 | assigned_group:Year),
  family = negbinomial()
)

# Current Reproduction (LINEAR Outcome)
bf_repro_conc_linear <- bf(
  num_offspring_all_births ~ 1 + 
    adjusted_SSB_z + # Only the linear term here
    adjusted_DSB_z +
    Age_s + I(Age_s^2) + 
    (1 | Animal_ID) + (1 | assigned_group:Year),
  family = negbinomial()
)

# 2. DEFINE PRIORS (Shared for both linear and quadratic)
priors_repro_conc <- c(
  # Repro Priors (Negative Binomial)
  set_prior("normal(0, 2)",      class = "Intercept"),
  set_prior("normal(0, 1)",      class = "b"),
  set_prior("exponential(1)",    class = "sd"),
  set_prior("exponential(1)",    class = "shape")
)

# 3. RUN QUADRATIC MODEL
fit_concurrent_repro_quad <- brm(
  bf_repro_conc_quad,
  data    = df_concurrent,
  prior   = priors_repro_conc,
  seed    = 1,
  chains  = 4, 
  cores   = 4,
  warmup  = 1500, 
  iter    = 4000,
  control = list(adapt_delta = 0.99, max_treedepth = 15),
  save_pars = save_pars(all = TRUE)
)

summary(fit_concurrent_repro_quad)

# 4. RUN LINEAR MODEL
fit_concurrent_repro_linear <- brm(
  bf_repro_conc_linear,
  data    = df_concurrent,
  prior   = priors_repro_conc,
  seed    = 1,
  chains  = 4, 
  cores   = 4,
  warmup  = 1500, 
  iter    = 4000,
  control = list(adapt_delta = 0.99, max_treedepth = 15),
  save_pars = save_pars(all = TRUE)
)

summary(fit_concurrent_repro_linear)


# MODEL 2: CONCURRENT RANK
cat("\nRunning Joint Models: Concurrent Rank\n")

# 1. DEFINE FORMULAS
# Current Rank (QUADRATIC Outcome)
bf_rank_conc_quad <- bf(
  DS_sc ~ 1 + 
    adjusted_SSB_z + I(adjusted_SSB_z^2) + 
    adjusted_DSB_z +
    Age_s + I(Age_s^2) +
    (1 | Animal_ID) + (1 | assigned_group:Year),
  family = gaussian()
)

# Current Rank (LINEAR Outcome)
bf_rank_conc_linear <- bf(
  DS_sc ~ 1 + 
    adjusted_SSB_z + # Only the linear term here
    adjusted_DSB_z +
    Age_s + I(Age_s^2) +
    (1 | Animal_ID) + (1 | assigned_group:Year),
  family = gaussian()
)

# 2. DEFINE PRIORS (Shared for both linear and quadratic)
priors_rank_conc <- c(
  # Rank Outcome
  set_prior("normal(0, 1)",      class = "Intercept"),
  set_prior("normal(0, 1)",      class = "b"),
  set_prior("exponential(2)",    class = "sd"),
  set_prior("student_t(3, 0, 1)",class = "sigma")
)

# 3. RUN QUADRATIC MODEL
fit_concurrent_rank_quad <- brm(
  bf_rank_conc_quad,
  data    = df_concurrent,
  prior   = priors_rank_conc,
  seed    = 1,
  chains  = 4, 
  cores   = 4,
  warmup  = 1500, 
  iter    = 4000,
  control = list(adapt_delta = 0.99, max_treedepth = 15),
  save_pars = save_pars(all = TRUE)
)

summary(fit_concurrent_rank_quad)

# 4. RUN LINEAR MODEL
fit_concurrent_rank_linear <- brm(
  bf_rank_conc_linear,
  data    = df_concurrent,
  prior   = priors_rank_conc,
  seed    = 1,
  chains  = 4, 
  cores   = 4,
  warmup  = 1500, 
  iter    = 4000,
  control = list(adapt_delta = 0.99, max_treedepth = 15),
  save_pars = save_pars(all = TRUE)
)

summary(fit_concurrent_rank_linear)

# Analysis of future reproduction (LINEAR) ####
# 1. Define the Dataframes for each time step
repro_one_df   <- analysis_df %>% filter(Year %in% c(2017, 2019, 2020))
repro_two_df   <- analysis_df %>% filter(Year %in% c(2017, 2019, 2020))
repro_three_df <- analysis_df %>% filter(Year %in% c(2017, 2019, 2020))
repro_four_df  <- analysis_df %>% filter(Year %in% c(2017, 2019))
repro_five_df  <- analysis_df %>% filter(Year %in% c(2017))
repro_six_df   <- analysis_df %>% filter(Year %in% c(2017))

# Put them into a list for iteration
df_list_repro <- list(
  repro_one_df, repro_two_df, repro_three_df, 
  repro_four_df, repro_five_df, repro_six_df
)

# Scale SSB and DSB across the DFs in the list
for (i in 1:length(df_list_repro)) {
  df_list_repro[[i]]$adjusted_SSB_z <- as.numeric(scale(df_list_repro[[i]]$adjusted_SSB))
  df_list_repro[[i]]$adjusted_DSB_z <- as.numeric(scale(df_list_repro[[i]]$adjusted_DSB))
}

# 3. Run Models in a Loop
model_results_linear <- list()

for(i in 1:6) { 
  
  cat(paste0("\nRunning Joint Model (Linear) for lag: ", i, "\n"))
  
  col_name_offspring <- paste0("offspring_plus_", i)
  resp_name_ssb <- "adjustedSSBz"      
  resp_name_dsb <- "adjustedDSBz"
  resp_name_off <- paste0("offspringplus", i)
  
  # Offspring
  bf_offspring <- bf(
    as.formula(paste0(
      col_name_offspring, " ~ 1 + ",
      "adjusted_SSB_z + ",
      "adjusted_DSB_z + Age_s + I(Age_s^2) + ",
      "(1 | Animal_ID) + (1 | assigned_group:Year)"
    )),
    family = negbinomial()
  )
  
  priors_joint <- c(
    # --- Offspring Priors ---
    set_prior("normal(0, 2)",      class = "Intercept"),
    set_prior("normal(0, 1)",      class = "b"),
    set_prior("exponential(1)",    class = "sd"),
    set_prior("exponential(1)",    class = "shape")
  )
  
  fit <- brm(
    bf_offspring,
    data    = df_list_repro[[i]],
    prior   = priors_joint,
    seed    = 1,
    chains  = 4, 
    cores   = 4,
    warmup  = 1500, 
    iter    = 4000,
    control = list(adapt_delta = 0.99, max_treedepth = 15),
    save_pars = save_pars(all = TRUE)
  )
  
  model_results_linear[[paste0("plus_", i)]] <- fit
}

summary(model_results_linear[[1]])
summary(model_results_linear[[2]])

# Analysis of future rank (LINEAR) ####
rank_one_df   <- analysis_df %>% filter(Year %in% c(2019, 2023, 2024))
rank_two_df   <- analysis_df %>% filter(Year %in% c(2017, 2023))
rank_three_df <- analysis_df %>% filter(Year %in% c(2017, 2020))
rank_four_df  <- analysis_df %>% filter(Year %in% c(2019, 2020))
rank_five_df  <- analysis_df %>% filter(Year %in% c(2019, 2020))
rank_six_df   <- analysis_df %>% filter(Year %in% c(2017, 2019))
rank_seven_df   <- analysis_df %>% filter(Year %in% c(2017))
rank_eight_df   <- analysis_df %>% filter(Year %in% c(2017))

df_list_rank <- list(
  rank_one_df, rank_two_df, rank_three_df, 
  rank_four_df, rank_five_df, rank_six_df, 
  rank_seven_df, rank_eight_df
)

for (i in 1:length(df_list_rank)) {
  df_list_rank[[i]]$adjusted_SSB_z <- as.numeric(scale(df_list_rank[[i]]$adjusted_SSB))
  df_list_rank[[i]]$adjusted_DSB_z <- as.numeric(scale(df_list_rank[[i]]$adjusted_DSB))
}

model_results_rank_linear <- list()

for(i in 1:8) { 
  
  cat(paste0("\nRunning Rank Model (Linear) for lag: ", i, "\n"))
  
  col_name_rank <- paste0("Rank_Future_plus_", i, "_s")
  resp_name_ssb <- "adjustedSSBz"      
  resp_name_dsb <- "adjustedDSBz"
  resp_name_rank <- paste0("RankFutureplus", i, "s")
  
  bf_rank <- bf(
    as.formula(paste0(
      col_name_rank, " ~ 1 + ",
      "adjusted_SSB_z + ",
      "adjusted_DSB_z + DS_sc + Age_s + I(Age_s^2) +",
      "(1 | Animal_ID) + (1 | assigned_group:Year)"
    )),
    family = gaussian()
  )
  
  priors_joint <- c(
    # --- Rank Priors ---
    set_prior("normal(0, 1)",      class = "Intercept"),
    set_prior("normal(0, 1)",      class = "b"),
    set_prior("exponential(2)",    class = "sd"),
    set_prior("student_t(3, 0, 1)",class = "sigma")
  )
  
  fit <- brm(
    bf_rank,
    data    = df_list_rank[[i]],
    prior   = priors_joint,
    seed    = 1,
    chains  = 4, 
    cores   = 4,
    warmup  = 1500, 
    iter    = 4000,
    control = list(adapt_delta = 0.99, max_treedepth = 15),
    save_pars = save_pars(all = TRUE)
  )
  
  model_results_rank_linear[[paste0("plus_", i)]] <- fit
}

# Analysis of future reproduction - (QUADRATIC TERM) ####
model_results <- list()

for(i in 1:6) {
  
  cat(paste0("\nRunning Joint Model for lag: ", i, "\n"))
  
  col_name_offspring <- paste0("offspring_plus_", i)
  resp_name_ssb <- "adjustedSSBz"
  resp_name_dsb <- "adjustedDSBz"
  resp_name_off <- paste0("offspringplus", i)
  
  bf_offspring <- bf(
    as.formula(paste0(
      col_name_offspring, " ~ 1 + ",
      "adjusted_SSB_z + I(adjusted_SSB_z^2) + ",
      "adjusted_DSB_z + ",
      "Age_s + I(Age_s^2) + ",
      "(1 | Animal_ID) + (1 | assigned_group:Year)"
    )),
    family = negbinomial()
  )
  
  priors_joint <- c(
    # --- Offspring ---
    set_prior("normal(0, 2)",      class = "Intercept"),
    set_prior("normal(0, 1)",      class = "b"),
    set_prior("exponential(1)",    class = "sd"),
    set_prior("exponential(1)",    class = "shape")
  )
  
  fit <- brm(
    formula = bf_offspring,
    data    = df_list_repro[[i]],
    prior   = priors_joint,
    seed    = 1,
    chains  = 4,
    cores   = 4,
    warmup  = 1500,
    iter    = 4000,
    control = list(adapt_delta = 0.99, max_treedepth = 15),
    save_pars = save_pars(all = TRUE)
  )
  
  model_results[[paste0("plus_", i)]] <- fit
}

# Calculate the probability of direction for all parameters
p_direction(model_results[[1]])
summary(model_results[[1]])

# MULTI-PANEL PLOT: FUTURE REPRODUCTION (YEARS +1 to +3)
# 1. EXTRACT MODELS FROM THE LIST
fit_repro_1 <- model_results[["plus_1"]]
fit_repro_2 <- model_results[["plus_2"]]
fit_repro_3 <- model_results[["plus_3"]]

# 2. CREATE SSB PLOTS (Row 1 - Blue)
# PLOTTING FUNCTION: Create Spaghetti Plot (Variable Y-Limits, Linetype, and Thickness)
create_spaghetti_plot <- function(model, x_var, y_resp, x_label, y_label, line_color, y_limits, 
                                  mean_linetype = "solid", mean_linewidth = 0.7) {
  
  # A. GENERATE DATA GRID
  ce <- conditional_effects(model, effects = x_var, resp = y_resp, resolution = 100)
  dat_grid <- ce[[1]] 
  
  # B. GENERATE POSTERIOR DRAWS
  post_draws <- posterior_epred(
    model, 
    newdata = dat_grid, 
    resp = y_resp, 
    ndraws = 300,
    re_formula = NA # Use population average
  )
  
  # C. WRANGLE DATA
  draws_df <- data.frame(t(post_draws)) 
  draws_merged <- cbind(dat_grid, draws_df) 
  draws_long <- draws_merged %>%
    pivot_longer(
      cols = starts_with("X"), 
      names_to = "draw", 
      values_to = "prediction"
    )
  
  # D. PLOT
  p <- ggplot() +
    # 1. The Lines
    geom_line(data = draws_long, 
              aes(x = .data[[x_var]], y = prediction, group = draw),
              color = line_color, 
              alpha = 0.1,        
              linewidth = 0.3) +
    
    # 2. The Mean Estimate (Adjustable linetype and thickness)
    geom_line(data = dat_grid,
              aes(x = .data[[x_var]], y = estimate__),
              color = "black", 
              linewidth = mean_linewidth, # Controlled by function argument
              linetype = mean_linetype) + # Controlled by function argument
    
    # 3. Styling
    labs(x = x_label, y = y_label) +
    
    theme_cowplot() +
    theme(
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "black")
    ) +
    
    # 4. Zoom (Variable Limits)
    coord_cartesian(ylim = y_limits)
  
  return(p)
}

# Using expression() to format the t subscript for the axis labels.
p_ssb_1 <- create_spaghetti_plot(
  model      = fit_repro_1,
  x_var      = "adjusted_SSB_z",
  y_resp     = "offspringplus1",
  x_label    = expression("SSB " * t["year"] * " (std)"),
  y_label    = expression("No. offspring sired " * t["year"] * " + 1 (std)"),
  line_color = "goldenrod2",
  y_limits   = c(0, 2) # Keeping axes uniform across panels for direct comparison
) + 
  ggtitle("Lag: +1 Year") + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

p_ssb_2 <- create_spaghetti_plot(
  model      = fit_repro_2,
  x_var      = "adjusted_SSB_z",
  y_resp     = "offspringplus2",
  x_label    = expression("SSB " * t["year"] * " (std)"),
  y_label    = expression("No. offspring sired " * t["year"] * " + 2 (std)"),
  line_color = "goldenrod2",
  y_limits   = c(0, 2)
) + 
  ggtitle("Lag: +2 Years") + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

p_ssb_3 <- create_spaghetti_plot(
  model      = fit_repro_3,
  x_var      = "adjusted_SSB_z",
  y_resp     = "offspringplus3",
  x_label    = expression("SSB " * t["year"] * " (std)"),
  y_label    = expression("No. offspring sired " * t["year"] * " + 3 (std)"),
  line_color = "goldenrod2",
  y_limits   = c(0, 2),
  mean_linetype = 'dashed'
) + 
  ggtitle("Lag: +3 Years") + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

# 3. CREATE DSB PLOTS (Row 2 - Vermillion)
# Applying the same expression formatting for the DSB axis labels

p_dsb_1 <- create_spaghetti_plot(
  model      = fit_repro_1,
  x_var      = "adjusted_DSB_z",
  y_resp     = "offspringplus1",
  x_label    = expression("DSB " * t["year"] * " (std)"),
  y_label    = expression("No. offspring sired " * t["year"] * " + 1 (std)"),
  line_color = "darkorchid2",
  y_limits   = c(0, 9)
)

p_dsb_2 <- create_spaghetti_plot(
  model      = fit_repro_2,
  x_var      = "adjusted_DSB_z",
  y_resp     = "offspringplus2",
  x_label    = expression("DSB " * t["year"] * " (std)"),
  y_label    = expression("No. offspring sired " * t["year"] * " + 2 (std)"),
  line_color = "darkorchid2",
  y_limits   = c(0, 9)
)

p_dsb_3 <- create_spaghetti_plot(
  model      = fit_repro_3,
  x_var      = "adjusted_DSB_z",
  y_resp     = "offspringplus3",
  x_label    = expression("DSB " * t["year"] * " (std)"),
  y_label    = expression("No. offspring sired " * t["year"] * " + 3 (std)"),
  line_color = "darkorchid2",
  y_limits   = c(0, 9),
  mean_linetype = 'dashed'
)

# 1. Arrange the 6 plots and automatically add A-F tags
repro_grid <- plot_grid(
  p_ssb_1, p_ssb_2, p_ssb_3,
  p_dsb_1, p_dsb_2, p_dsb_3,
  nrow = 2, 
  ncol = 3,
  labels = "AUTO", # Automatically adds A, B, C, D, E, F
  label_size = 14,
  label_fontface = "bold"
)

# 2. Create the title and subtitle as their own objects
title <- ggdraw() + draw_label("Future Reproduction Predictions (Years +1 to +3)", fontface = "bold", size = 16)
subtitle <- ggdraw() + draw_label("Top Row: Quadratic effect of SSB | Bottom Row: Linear effect of DSB", size = 12)

# 3. Stack the title, subtitle, and grid together
combined_repro_plot <- plot_grid(
  # title, 
  # subtitle, 
  repro_grid, 
  ncol = 1, 
  rel_heights = c(0.1, 0.05, 1) # Keeps the titles from taking up too much vertical space
)

combined_repro_plot

# Analysis of future rank - (QUADRATIC TERM) ####
model_results_rank <- list()

for(i in 1:8) { 
  
  cat(paste0("\nRunning Model for Rank Plus ", i, "\n"))
  
  col_name_rank <- paste0("Rank_Future_plus_", i, "_s")
  resp_name_ssb <- "adjustedSSBz"      
  resp_name_dsb <- "adjustedDSBz"
  resp_name_rank <- paste0("RankFutureplus", i, "s")
  
  bf_rank <- bf(
    as.formula(paste0(
      col_name_rank, " ~ 1 + ",
      "adjusted_SSB_z + I(adjusted_SSB_z^2) + ", 
      "adjusted_DSB_z + ", 
      "DS_sc + Age_s + I(Age_s^2) +",
      "(1 | Animal_ID) + (1 | assigned_group:Year)"
    )),
    family = gaussian()
  )
  
  priors_joint <- c(
    # --- Rank ---
    set_prior("normal(0, 1)",      class = "Intercept"),
    set_prior("normal(0, 1)",      class = "b"),
    set_prior("exponential(2)",    class = "sd"),
    set_prior("student_t(3, 0, 1)",class = "sigma")
  )
  
  fit <- brm(
    bf_rank,
    data    = df_list_rank[[i]],
    prior   = priors_joint,
    seed    = 1,
    chains  = 4, 
    cores   = 4,
    warmup  = 1500, 
    iter    = 4000,
    control = list(adapt_delta = 0.99),
    save_pars = save_pars(all = TRUE)
  )
  
  model_results_rank[[paste0("plus_", i)]] <- fit
}

summary(model_results_rank[[1]])
p_direction(model_results_rank[[1]])

# MULTI-PANEL PLOT: FUTURE RANK (YEARS +1 to +3)
# 1. EXTRACT MODELS FROM THE LIST
fit_rank_1 <- model_results_rank[["plus_1"]]
fit_rank_2 <- model_results_rank[["plus_2"]]
fit_rank_3 <- model_results_rank[["plus_3"]]

# 2. CREATE SSB PLOTS (Row 1 - Blue)
p_rank_ssb_1 <- create_spaghetti_plot(
  model      = fit_rank_1,
  x_var      = "adjusted_SSB_z",
  y_resp     = "RankFutureplus1s",
  x_label    = expression("SSB " * t["year"] * " (scaled)"),
  y_label    = expression("Rank " * t["year"] * " + 1 (scaled)"),
  line_color = "goldenrod2",
  y_limits   = c(-2, 1.5) 
) + 
  ggtitle("Lag: +1 Year") + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

saveRDS(p_rank_ssb_1, file = "/SSB_delayed_fitness/future_rank.rds")

p_rank_ssb_2 <- create_spaghetti_plot(
  model      = fit_rank_2,
  x_var      = "adjusted_SSB_z",
  y_resp     = "RankFutureplus2s",
  x_label    = expression("SSB " * t["year"] * " (scaled)"),
  y_label    = expression("Rank " * t["year"] * " + 2 (scaled)"),
  line_color = "goldenrod2",
  y_limits   = c(-2, 1.5)
) + 
  ggtitle("Lag: +2 Years") + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

p_rank_ssb_3 <- create_spaghetti_plot(
  model      = fit_rank_3,
  x_var      = "adjusted_SSB_z",
  y_resp     = "RankFutureplus3s",
  x_label    = expression("SSB " * t["year"] * " (scaled)"),
  y_label    = expression("Rank " * t["year"] * " + 3 (scaled)"),
  line_color = "goldenrod2",
  y_limits   = c(-2, 1.5),
  mean_linetype = 'dashed'
) + 
  ggtitle("Lag: +3 Years") + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

# 3. CREATE DSB PLOTS (Row 2 - Vermillion)
p_rank_dsb_1 <- create_spaghetti_plot(
  model      = fit_rank_1,
  x_var      = "adjusted_DSB_z",
  y_resp     = "RankFutureplus1s",
  x_label    = expression("DSB " * t["year"] * " (scaled)"),
  y_label    = expression("Rank " * t["year"] * " + 1 (scaled)"),
  line_color = "darkorchid2",
  y_limits   = c(-2, 3)
)  + 
  ggtitle("Lag: +1 Year") + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))


saveRDS(p_rank_dsb_1, file = "/SSB_delayed_fitness/future_rank_DSB_plus1.rds")

p_rank_dsb_2 <- create_spaghetti_plot(
  model      = fit_rank_2,
  x_var      = "adjusted_DSB_z",
  y_resp     = "RankFutureplus2s",
  x_label    = expression("DSB " * t["year"] * " (scaled)"),
  y_label    = expression("Rank " * t["year"] * " + 2 (scaled)"),
  line_color = "darkorchid2",
  y_limits   = c(-2, 3)
)  + 
  ggtitle("Lag: +2 Years") + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))


saveRDS(p_rank_dsb_2, file = "/SSB_delayed_fitness/future_rank_DSB_plus2.rds")

p_rank_dsb_3 <- create_spaghetti_plot(
  model      = fit_rank_3,
  x_var      = "adjusted_DSB_z",
  y_resp     = "RankFutureplus3s",
  x_label    = expression("DSB " * t["year"] * " (scaled)"),
  y_label    = expression("Rank " * t["year"] * " + 3 (scaled)"),
  line_color = "darkorchid2",
  y_limits   = c(-2, 3),
  mean_linetype = 'dashed'
)

# 4. ASSEMBLE THE 4-PANEL PLOT USING COWPLOT
# Arrange the 4 plots and automatically add A-F tags
rank_grid <- plot_grid(
  p_rank_ssb_1, p_rank_ssb_2,
  p_rank_dsb_1, p_rank_dsb_2,
  nrow = 2, 
  ncol = 2,
  labels = "AUTO", 
  label_size = 14,
  label_fontface = "bold"
)

# Create the title and subtitle as their own objects
title <- ggdraw() + draw_label("Future Rank Predictions (Years +1 to +3)", fontface = "bold", size = 16)
subtitle <- ggdraw() + draw_label("Top Row: Quadratic effect of SSB | Bottom Row: Linear effect of DSB", size = 12)

# Stack the title, subtitle, and grid together
combined_rank_plot <- plot_grid(
  title, 
  subtitle, 
  rank_grid, 
  ncol = 1, 
  rel_heights = c(0.1, 0.05, 1) # Keeps the titles from taking up too much vertical space
)

# 5. VIEW AND SAVE
print(combined_rank_plot)

# Plot the results and save tables ####
# Generate unified model output tables for all concurrent and future models

# --- 1. Define extraction function ---
# Can be looped over lists of univariate brms models.
extract_brms_summary <- function(model, model_name, response_type) {
  
  # Extract brms summary (contains ESS and Rhat)
  mod_sum <- summary(model)
  
  # ---------------------------------------------------------
  # Extract Divergent Transitions
  # ---------------------------------------------------------
  div_trans <- tryCatch({
    np <- brms::nuts_params(model)
    sum(np$Value[np$Parameter == "divergent__"])
  }, error = function(e) { NA }) # Returns NA if it can't be calculated
  
  # Extract and Format Fixed Effects
  fixed_df <- as.data.frame(mod_sum$fixed) %>%
    tibble::rownames_to_column("Raw_Parameter") %>%
    dplyr::mutate(
      Effect_Type = "Fixed",
      Group = NA_character_
    )
  
  # Extract and Format Random Effects (Conditional Check)
  if (length(mod_sum$random) > 0) {
    ran_df_list <- lapply(names(mod_sum$random), function(grp) {
      as.data.frame(mod_sum$random[[grp]]) %>%
        tibble::rownames_to_column("Raw_Parameter") %>%
        dplyr::mutate(
          Effect_Type = "Random",
          Group = grp
        )
    })
    ran_df <- dplyr::bind_rows(ran_df_list)
  } else {
    ran_df <- data.frame() # Empty dataframe if no random effects
  }
  
  # Combine Fixed and Random
  raw_combined <- dplyr::bind_rows(fixed_df, ran_df) %>%
    dplyr::select(
      Effect_Type,
      Group,
      Raw_Parameter,
      Estimate,
      Est_Error = `Est.Error`,       # brms column name
      Lower_95_CI = `l-95% CI`,      # brms column name
      Upper_95_CI = `u-95% CI`,      # brms column name
      ESS_Bulk = Bulk_ESS,           # Extract Bulk Effective Sample Size
      ESS_Tail = Tail_ESS,           # Extract Tail Effective Sample Size
      Rhat                           # Extract Rhat convergence diagnostic
    )
  
  # Calculate Max Rhat across this specific model
  model_max_rhat <- max(raw_combined$Rhat, na.rm = TRUE)
  
  # Clean, Label, and Format the Table
  clean_table <- raw_combined %>%
    dplyr::mutate(
      Model_Name = model_name,
      Response_Variable = response_type,
      Parameter = Raw_Parameter,
      
      # Clean random effect Standard Deviation and Correlation wrappers
      Parameter = str_replace_all(Parameter, "^cor\\((.*)\\)$", "Cor: \\1"),
      Parameter = ifelse(str_detect(Parameter, "^Cor: "), str_replace_all(Parameter, ",", " <-> "), Parameter),
      Parameter = str_replace_all(Parameter, "^sd\\((.*)\\)$", "SD: \\1"), 
      
      # Clean specific predictor labels (handling brms syntax for polynomials)
      Parameter = str_replace_all(Parameter, "adjusted_SSB_z", "SSB (scaled)"),
      Parameter = str_replace_all(Parameter, "ISSB \\(scaled\\)E2", "SSB^2 (scaled)"),
      Parameter = str_replace_all(Parameter, "adjusted_DSB_z", "DSB (scaled)"),
      Parameter = str_replace_all(Parameter, "Age_s", "Age (scaled)"),
      Parameter = str_replace_all(Parameter, "IAge \\(scaled\\)E2", "Age^2 (scaled)"),
      Parameter = str_replace_all(Parameter, "DS_sc", "Current Rank (scaled)"),
      Parameter = str_replace_all(Parameter, "assigned_group", "Group: "),
      
      # ---------------------------------------------------------
      # Convergence Flags
      # ---------------------------------------------------------
      Model_Divergent_Transitions = div_trans,
      Model_Max_Rhat = model_max_rhat,
      
      # Flag individual parameter if Rhat > 1.05 or ESS < 400
      Param_Poor_Convergence = ifelse(Rhat > 1.05 | ESS_Bulk < 400 | ESS_Tail < 400, "Yes", "No"),
      
      # Flag the entire model if there are any divergences or high max Rhat
      Model_Convergence_Warning = ifelse(Model_Divergent_Transitions > 0 | Model_Max_Rhat > 1.05, "Yes", "No")
    ) %>%
    # Reorder columns logically
    dplyr::select(
      Model_Name, Response_Variable, Effect_Type, Group, Parameter, 
      Estimate, Est_Error, Lower_95_CI, Upper_95_CI, ESS_Bulk, ESS_Tail, Rhat,
      Param_Poor_Convergence, Model_Divergent_Transitions, Model_Max_Rhat, Model_Convergence_Warning
    )
  
  return(clean_table)
}

# --- 2. Extract Standalone (Concurrent) Models ---
cat("\nExtracting summaries for Concurrent models...\n")
table_conc_repro_lin  <- extract_brms_summary(fit_concurrent_repro_linear, "Concurrent Reproduction (Linear)", "Offspring")
table_conc_repro_quad <- extract_brms_summary(fit_concurrent_repro_quad, "Concurrent Reproduction (Quadratic)", "Offspring")
table_conc_rank_lin   <- extract_brms_summary(fit_concurrent_rank_linear, "Concurrent Rank (Linear)", "Rank")
table_conc_rank_quad  <- extract_brms_summary(fit_concurrent_rank_quad, "Concurrent Rank (Quadratic)", "Rank")

# --- 3. Extract Lists of Future Models ---
cat("Extracting summaries for Future model lists...\n")

# Future Repro Linear (Lags 1-6)
if(exists("model_results_linear") && length(model_results_linear) > 0) {
  table_fut_repro_lin <- purrr::map_dfr(names(model_results_linear), function(mod_name) {
    extract_brms_summary(model_results_linear[[mod_name]], paste0("Future Repro Linear (", mod_name, ")"), "Offspring")
  })
} else { table_fut_repro_lin <- data.frame() }

# Future Repro Quadratic (Lags 1-6)
if(exists("model_results") && length(model_results) > 0) {
  table_fut_repro_quad <- purrr::map_dfr(names(model_results), function(mod_name) {
    extract_brms_summary(model_results[[mod_name]], paste0("Future Repro Quad (", mod_name, ")"), "Offspring")
  })
} else { table_fut_repro_quad <- data.frame() }

# Future Rank Linear (Lags 1-8)
if(exists("model_results_rank_linear") && length(model_results_rank_linear) > 0) {
  table_fut_rank_lin <- purrr::map_dfr(names(model_results_rank_linear), function(mod_name) {
    extract_brms_summary(model_results_rank_linear[[mod_name]], paste0("Future Rank Linear (", mod_name, ")"), "Rank")
  })
} else { table_fut_rank_lin <- data.frame() }

# Future Rank Quadratic (Lags 1-8)
if(exists("model_results_rank") && length(model_results_rank) > 0) {
  table_fut_rank_quad <- purrr::map_dfr(names(model_results_rank), function(mod_name) {
    extract_brms_summary(model_results_rank[[mod_name]], paste0("Future Rank Quad (", mod_name, ")"), "Rank")
  })
} else { table_fut_rank_quad <- data.frame() }


# --- 4. Bind All Tables and Apply Final Formatting ---
cat("Combining tables and rounding numeric values...\n")

final_combined_table <- dplyr::bind_rows(
  table_conc_repro_lin,
  table_conc_repro_quad,
  table_conc_rank_lin,
  table_conc_rank_quad,
  table_fut_repro_lin,
  table_fut_repro_quad,
  table_fut_rank_lin,
  table_fut_rank_quad
) %>%
  dplyr::mutate(
    # Round estimates, errors, CIs, and Rhat to 3 decimal places
    dplyr::across(c(Estimate, Est_Error, Lower_95_CI, Upper_95_CI, Rhat, Model_Max_Rhat), ~ round(.x, 3)),
    # Round ESS to whole numbers since they represent sample sizes
    dplyr::across(c(ESS_Bulk, ESS_Tail), ~ round(.x, 0)) 
  )

# --- 5. Export the Final Unified Table ---
output_dir <- "/SSB_delayed_fitness/" 
if(!dir.exists(output_dir)) dir.create(output_dir)

csv_out_path <- file.path(output_dir, "All_Models_Unified_Summary_Table_fixed&random.csv")

# write.csv(final_combined_table, csv_out_path, row.names = FALSE)
cat("SUCCESS! Output table saved to:", csv_out_path, "\n")

#Save the output models####
# 1. Define and create the directory for saving models
model_dir <- "/SSB_delayed_fitness/"
if(!dir.exists(model_dir)) dir.create(model_dir)

cat(paste0("\nSaving model objects to: ", model_dir, "\n"))

# 2. SAVE CONCURRENT MODELS (Individual brmsfit objects)
# Save Concurrent Reproduction (Quadratic)
if(exists("fit_concurrent_repro_quad")) {
  saveRDS(fit_concurrent_repro_quad, 
          file = file.path(model_dir, "fit_concurrent_repro_quad.rds"))
  cat(" - Saved fit_concurrent_repro_quad.rds\n")
}

# Save Concurrent Reproduction (Linear)
if(exists("fit_concurrent_repro_linear")) {
  saveRDS(fit_concurrent_repro_linear, 
          file = file.path(model_dir, "fit_concurrent_repro_linear.rds"))
  cat(" - Saved fit_concurrent_repro_linear.rds\n")
}

# Save Concurrent Rank (Quadratic)
if(exists("fit_concurrent_rank_quad")) {
  saveRDS(fit_concurrent_rank_quad, 
          file = file.path(model_dir, "fit_concurrent_rank_quad.rds"))
  cat(" - Saved fit_concurrent_rank_quad.rds\n")
}

# Save Concurrent Rank (Linear)
if(exists("fit_concurrent_rank_linear")) {
  saveRDS(fit_concurrent_rank_linear, 
          file = file.path(model_dir, "fit_concurrent_rank_linear.rds"))
  cat(" - Saved fit_concurrent_rank_linear.rds\n")
}


# 3. SAVE FUTURE MODELS (Lists of brmsfit objects)
# Saving the entire list at once is much cleaner than saving each lag individually. 
# When you load it via readRDS(), it will return the fully intact list.

# Save Future Reproduction (Linear)
if(exists("model_results_linear")) {
  saveRDS(model_results_linear, 
          file = file.path(model_dir, "list_future_repro_linear.rds"))
  cat(" - Saved list_future_repro_linear.rds\n")
}

# Save Future Rank (Linear)
if(exists("model_results_rank_linear")) {
  saveRDS(model_results_rank_linear, 
          file = file.path(model_dir, "list_future_rank_linear.rds"))
  cat(" - Saved list_future_rank_linear.rds\n")
}

# Save Future Reproduction (Quadratic)
if(exists("model_results")) {
  saveRDS(model_results, 
          file = file.path(model_dir, "list_future_repro_quad.rds"))
  cat(" - Saved list_future_repro_quad.rds\n")
}

# Save Future Rank (Quadratic)
if(exists("model_results_rank")) {
  saveRDS(model_results_rank, 
          file = file.path(model_dir, "list_future_rank_quad.rds"))
  cat(" - Saved list_future_rank_quad.rds\n")
}

cat("\nAll available models successfully saved as .rds files.\n")

#Load the models ####
# 1. Define the directory where models are saved
# This must match the directory used in the saving script
model_dir <- "/SSB_delayed_fitness/"
cat(paste0("\nLoading model objects from: ", model_dir, "\n"))

# 2. LOAD CONCURRENT MODELS (Individual brmsfit objects)
# Define the file path for the Concurrent Reproduction (Quadratic) model
file_conc_repro_quad <- file.path(model_dir, "fit_concurrent_repro_quad.rds")
# Check if the file exists before loading to prevent errors
if(file.exists(file_conc_repro_quad)) {
  fit_concurrent_repro_quad <- readRDS(file_conc_repro_quad)
  cat(" - Loaded fit_concurrent_repro_quad\n")
} else {
  warning(paste("File not found:", file_conc_repro_quad))
}

# Define the file path for the Concurrent Reproduction (Linear) model
file_conc_repro_linear <- file.path(model_dir, "fit_concurrent_repro_linear.rds")
if(file.exists(file_conc_repro_linear)) {
  fit_concurrent_repro_linear <- readRDS(file_conc_repro_linear)
  cat(" - Loaded fit_concurrent_repro_linear\n")
} else {
  warning(paste("File not found:", file_conc_repro_linear))
}

# Define the file path for the Concurrent Rank (Quadratic) model
file_conc_rank_quad <- file.path(model_dir, "fit_concurrent_rank_quad.rds")
if(file.exists(file_conc_rank_quad)) {
  fit_concurrent_rank_quad <- readRDS(file_conc_rank_quad)
  cat(" - Loaded fit_concurrent_rank_quad\n")
} else {
  warning(paste("File not found:", file_conc_rank_quad))
}

# Define the file path for the Concurrent Rank (Linear) model
file_conc_rank_linear <- file.path(model_dir, "fit_concurrent_rank_linear.rds")
if(file.exists(file_conc_rank_linear)) {
  fit_concurrent_rank_linear <- readRDS(file_conc_rank_linear)
  cat(" - Loaded fit_concurrent_rank_linear\n")
} else {
  warning(paste("File not found:", file_conc_rank_linear))
}

# 3. LOAD FUTURE MODELS (Lists of brmsfit objects)
# These will be loaded back as lists, just as they were saved.
# Define the file path for the Future Reproduction (Linear) list
file_list_future_repro_linear <- file.path(model_dir, "list_future_repro_linear.rds")
if(file.exists(file_list_future_repro_linear)) {
  model_results_linear <- readRDS(file_list_future_repro_linear)
  cat(" - Loaded list_future_repro_linear into 'model_results_linear'\n")
} else {
  warning(paste("File not found:", file_list_future_repro_linear))
}

# Define the file path for the Future Rank (Linear) list
file_list_future_rank_linear <- file.path(model_dir, "list_future_rank_linear.rds")
if(file.exists(file_list_future_rank_linear)) {
  model_results_rank_linear <- readRDS(file_list_future_rank_linear)
  cat(" - Loaded list_future_rank_linear into 'model_results_rank_linear'\n")
} else {
  warning(paste("File not found:", file_list_future_rank_linear))
}

# Define the file path for the Future Reproduction (Quadratic) list
file_list_future_repro_quad <- file.path(model_dir, "list_future_repro_quad.rds")
if(file.exists(file_list_future_repro_quad)) {
  # Loading back into the original object name 'model_results'
  model_results <- readRDS(file_list_future_repro_quad) 
  cat(" - Loaded list_future_repro_quad into 'model_results'\n")
} else {
  warning(paste("File not found:", file_list_future_repro_quad))
}

# Define the file path for the Future Rank (Quadratic) list
file_list_future_rank_quad <- file.path(model_dir, "list_future_rank_quad.rds")
if(file.exists(file_list_future_rank_quad)) {
  # Loading back into the original object name 'model_results_rank'
  model_results_rank <- readRDS(file_list_future_rank_quad)
  cat(" - Loaded list_future_rank_quad into 'model_results_rank'\n")
} else {
  warning(paste("File not found:", file_list_future_rank_quad))
}

cat("\nFinished loading models.\n")
summary(model_results[[1]])
summary(model_results_rank[[1]])
p_direction(model_results_rank[[3]])

#Model comparison between quadratic model and linear model when predicting future number offspring sired and future rank ####

# SECTION 1: Compare Future Reproduction Models (Linear vs. Quadratic) - However, Pareto K is too high to use this approach. See k-fold comparison below
# Create an empty list to store the LOO comparison results for reproduction models
repro_comparisons <- list()

# Loop through each model in the reproduction lists
# seq_along ensures we iterate over the exact number of models stored in the list
for (i in 1:2) {
  
  cat("\n--- Evaluating Reproduction Model Pair", i, "---\n")
  
  # Extract the i-th linear and quadratic model from their respective lists
  lin_mod <- model_results_linear[[i]]
  quad_mod <- model_results[[i]]
  
  # Calculate Leave-One-Out (LOO) cross-validation for the linear model
  # This estimates out-of-sample predictive performance
  loo_lin <- loo(lin_mod)
  
  # Calculate LOO for the quadratic model
  loo_quad <- loo(quad_mod)
  
  # Compare the two models based on their Expected Log Predictive Density (ELPD)
  # loo_compare ranks models from best to worst.
  # A model with an elpd_diff of 0 is the best fitting model in the comparison.
  comp <- loo_compare(loo_lin, loo_quad)
  
  # Store the comparison object in our results list for later extraction if needed
  repro_comparisons[[i]] <- comp
  
  # Print the comparison output to the console
  print(comp)
}


# SECTION 2: Compare Future Rank Models (Linear vs. Quadratic) - Pareto K is acceptable when using moment matching

#Initialise the dataframe
summary_stats_rank_loo <- data.frame()

for (i in 1:1) {
  
  lag_name <- paste0("Lag_", i)
  cat("\n--- Evaluating Rank Model Pair", i, "---\n")
  
  # 1. Extract the specific models for lag i
  lin_mod_rank  <- model_results_rank_linear[[i]]
  quad_mod_rank <- model_results_rank[[i]]
  
  # 2. Re-attach active C++
  lin_mod_rank$fit@stanmodel  <- rstan::stan_model(model_code = brms::stancode(lin_mod_rank))
  quad_mod_rank$fit@stanmodel <- rstan::stan_model(model_code = brms::stancode(quad_mod_rank))
  
  # 3. Calculate LOO with moment matching
  loo_lin_rank  <- loo(lin_mod_rank, moment_match = TRUE)
  loo_quad_rank <- loo(quad_mod_rank, moment_match = TRUE)
  
  # 4. Compare models
  comp_rank <- loo_compare(list(Linear = loo_lin_rank, Quadratic = loo_quad_rank))
  print(comp_rank)
  
  # 5. Extract metrics for publication table
  winner    <- rownames(comp_rank)[1]
  elpd_diff <- comp_rank[2, "elpd_diff"]
  se_diff   <- comp_rank[2, "se_diff"]
  lin_elpd  <- loo_lin_rank$estimates["elpd_loo", "Estimate"]
  quad_elpd <- loo_quad_rank$estimates["elpd_loo", "Estimate"]
  
  summary_stats_rank_loo <- rbind(summary_stats_rank_loo, data.frame(
    Lag           = lag_name,
    Winning_Model = winner,
    ELPD_Diff     = round(elpd_diff, 3),
    SE_Diff       = round(se_diff, 3),
    Linear_ELPD   = round(lin_elpd, 3),
    Quad_ELPD     = round(quad_elpd, 3)
  ))
}

# Print the final table to check
print(summary_stats_rank_loo)

# write.csv(summary_stats_rank_loo, file = "/SSB_delayed_fitness/Future_Rank_LOO_Comparisons_Supplementary.csv", row.names = FALSE)

# ==============================================================================
# Grouped K-Fold Cross-Validation: Future Reproduction
# ==============================================================================

library(brms)
library(rstan)
library(loo)
library(future)

# Sequential execution
plan(sequential)

num_repeats <- 3
k_folds     <- 10

kfold_repro_results <- list()

for (lag in 1:3) {
  
  prefix   <- paste0("Offspring_Lag_t", lag)
  resp_var <- paste0("offspringplus", lag)
  
  cat("\n==============================================\n")
  cat("PROCESSING:", prefix, "\n")
  cat("==============================================\n")
  
  # 1. Extract and recompile models
  lin_mod  <- model_results_linear[[lag]]
  quad_mod <- model_results[[lag]]
  
  cat("Recompiling Stan models for local architecture...\n")
  lin_mod$fit@stanmodel  <- rstan::stan_model(model_code = brms::stancode(lin_mod))
  quad_mod$fit@stanmodel <- rstan::stan_model(model_code = brms::stancode(quad_mod))
  
  model_data <- lin_mod$data
  
  comparison_list_repro <- list()
  summary_stats_repro   <- data.frame(
    Run               = 1:num_repeats,
    Linear_ELPD       = NA_real_,
    Quad_ELPD         = NA_real_,
    Raw_Difference    = NA_real_,
    Winning_Model     = NA_character_,
    LooComp_elpd_diff = NA_real_,
    LooComp_se_diff   = NA_real_
  )
  
  # --- Repeated Splits Loop ---
  for (r in 1:num_repeats) {
    
    cat(paste0("\n  >>> Run ", r, " of ", num_repeats, " (Seed: ", r * 123, ") <<<\n"))
    set.seed(r * 123)
    
    explicit_folds <- kfold_split_grouped(K = k_folds, x = model_data$Animal_ID)
    
    # Fit Linear Model (chains = 2, cores = 2 to keep RAM stable)
    cat("  Fitting Linear Model K-fold...\n")
    kfold_lin <- kfold(
      lin_mod,
      folds  = explicit_folds,
      resp   = resp_var,
      chains = 2,
      cores  = 2
    )
    gc() # Flush RAM immediately after Linear fold completes
    
    # Fit Quadratic Model
    cat("  Fitting Quadratic Model K-fold...\n")
    kfold_quad <- kfold(
      quad_mod,
      folds  = explicit_folds,
      resp   = resp_var,
      chains = 2,
      cores  = 2
    )
    gc() # Flush RAM immediately after Quadratic fold completes
    
    # Model Comparison
    comp <- loo_compare(list(Linear = kfold_lin, Quadratic = kfold_quad))
    comparison_list_repro[[r]] <- comp
    
    print(comp)
    
    # Extract metrics
    lin_val  <- kfold_lin$estimates["elpd_kfold", "Estimate"]
    quad_val <- kfold_quad$estimates["elpd_kfold", "Estimate"]
    
    summary_stats_repro$Linear_ELPD[r]       <- lin_val
    summary_stats_repro$Quad_ELPD[r]         <- quad_val
    summary_stats_repro$Raw_Difference[r]    <- quad_val - lin_val
    summary_stats_repro$Winning_Model[r]     <- rownames(comp)[1]
    summary_stats_repro$LooComp_elpd_diff[r] <- comp[2, "elpd_diff"]
    summary_stats_repro$LooComp_se_diff[r]   <- comp[2, "se_diff"]
    
    # Clear large intermediate objects from memory before next run
    rm(kfold_lin, kfold_quad)
    gc()
  }
  
  # Multi-run synthesis
  mean_diff <- mean(summary_stats_repro$Raw_Difference)
  mean_se   <- mean(sapply(comparison_list_repro, function(x) x[2, "se_diff"]))
  sd_diff   <- sd(summary_stats_repro$Raw_Difference)
  
  synthesis_table <- data.frame(
    Metric = c("Mean ELPD Difference (Quad - Lin)", "Mean SE of Difference", "Stability (SD of Diffs)"),
    Value  = round(c(mean_diff, mean_se, sd_diff), 3)
  )
  
  cat(paste0("\nSynthesis Summary for ", prefix, ":\n"))
  print(synthesis_table)
  
  kfold_repro_results[[prefix]] <- list(
    repeated_summary = summary_stats_repro,
    synthesis        = synthesis_table,
    comparisons      = comparison_list_repro
  )
}