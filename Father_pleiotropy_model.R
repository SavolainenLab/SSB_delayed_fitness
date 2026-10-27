####### Functional Benefits of SSB: Sexually Antagonistic Pleiotropy Analysis #######
## Author: Louis Bell-Roberts
## Date: October 2026
## Description: Tests whether male same-sex sexual behaviour (SSB) expression is 
## associated with paternal lifetime reproductive success. Fits bivariate  
## models to estimate son behavioural propensities and right-censored mixed 
## models to evaluate paternal fecundity.
# ==============================================================================

# 0. Load Required Libraries
library(tidyverse)
library(brms)
library(cowplot)
library(lubridate)


cat("--- STEP 2: FITTING BAYESIAN PLEIOTROPY MODELS ---\n")

# ------------------------------------------------------------------------------
# 1. Ingest Clean Anonymised Dataset & Prepare Age Covariates
# ------------------------------------------------------------------------------
data_path <- "/SSB_delayed_fitness/father_pleiotropy_data.csv"

data <- readr::read_csv(data_path, show_col_types = FALSE) %>%
  dplyr::filter(Scan_Count >= 10, sex == "M", !is.na(Group)) %>%
  dplyr::mutate(Animal_ID = as.factor(Animal_ID))

# Partition age into within- and between-individual components (pure numeric scaling)
data <- data %>%
  dplyr::group_by(Animal_ID) %>%
  dplyr::mutate(
    age_B      = mean(Age, na.rm = TRUE),
    age_W      = Age - age_B,
    Group_Year = as.factor(paste(Group, Year, sep = "_"))
  ) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    age_B_s = as.numeric(scale(age_B)),
    age_W_s = as.numeric(scale(age_W))
  )

# ------------------------------------------------------------------------------
# 2. Fit Bivariate Model (Son SSB & DSB Propensities)
# ------------------------------------------------------------------------------
bf_SSB <- bf(
  SSB_Raw ~ 1 + age_B_s + age_W_s + offset(log(Scan_Count)) + 
    (1 | p | Animal_ID) + (1 | Group_Year),
  family = negbinomial() #distribution works well
)

bf_DSB <- bf(
  DSB_Raw ~ 1 + age_B_s + age_W_s + offset(log(Scan_Count)) + 
    (1 | p | Animal_ID) + (1 | Group_Year),
  family = negbinomial() #distribution works well
)

cat("Fitting Bivariate Model for Son SSB and DSB...\n")
biv_w_age_e2 <- brm(
  bf_SSB + bf_DSB + set_rescor(FALSE),
  data      = data, 
  chains    = 4, 
  cores     = 4,
  warmup    = 1000, 
  iter      = 4000, 
  thin      = 3,
  control   = list(adapt_delta = 0.95),
  seed      = 42,
  save_pars = save_pars(all = TRUE)
)

pp_check(biv_w_age_e2, resp = "DSBRaw", type = "bars")
pp_check(biv_w_age_e2, resp = "SSBRaw", type = "bars")
summary(biv_w_age_e2)

# ------------------------------------------------------------------------------
# 3. Extract Son Random Intercepts (BLUPs with Standard Errors)
# ------------------------------------------------------------------------------
ran_vals <- ranef(biv_w_age_e2)$Animal_ID

# Robust column matching (handles both 'SSBRaw_Intercept' and 'SSB_Raw_Intercept')
ssb_col <- grep("SSB.*Intercept", dimnames(ran_vals)[[3]], value = TRUE)[1]
dsb_col <- grep("DSB.*Intercept", dimnames(ran_vals)[[3]], value = TRUE)[1]

sons <- tibble::tibble(
  Animal_ID = as.character(rownames(ran_vals)),
  SSB_mu    = ran_vals[, "Estimate", ssb_col],
  SSB_sd    = ran_vals[, "Est.Error", ssb_col],
  DSB_mu    = ran_vals[, "Estimate", dsb_col],
  DSB_sd    = ran_vals[, "Est.Error", dsb_col]
)

# ------------------------------------------------------------------------------
# 4. Prepare Paternal Lifespan & Right-Censoring
# ------------------------------------------------------------------------------
max_age        <- 25
min_age        <- 4
study_end_date <- lubridate::ymd("2025-03-01")

new <- data %>%
  dplyr::mutate(Animal_ID = as.character(Animal_ID)) %>%
  dplyr::distinct(Animal_ID, .keep_all = TRUE) %>%
  dplyr::filter(!is.na(SireGen)) %>%
  dplyr::mutate(
    CPRC_DOB        = as.Date(SireGen_DOB),
    CPRC_death_date = as.Date(SireGen_death_date),
    max_date = dplyr::if_else(
      !is.na(CPRC_death_date) & CPRC_death_date < study_end_date,
      CPRC_death_date,
      study_end_date
    ),
    father_age_last = as.numeric(difftime(max_date, CPRC_DOB, units = "days")) / 365.2425,
    breeding        = pmin(father_age_last, max_age) - min_age,
    breeding        = pmax(breeding, 0),
    Offspring       = as.integer(SireGen_num_offspring),
    censored        = dplyr::if_else(breeding < 21 & max_date == study_end_date, 1, 0),
    SireGen         = as.factor(SireGen)
  ) %>%
  dplyr::filter(!is.na(breeding), breeding > 0, !is.na(Offspring)) %>%
  dplyr::select(Animal_ID, SireGen, breeding, Offspring, censored)

# ------------------------------------------------------------------------------
# 5. Join Son Estimates with Paternal Demographics & Fit Fecundity Model
# ------------------------------------------------------------------------------
joined_data <- dplyr::inner_join(sons, new, by = "Animal_ID")

priors_dad <- c(
  prior(normal(0, 5),   class = Intercept),
  prior(normal(0, 0.5), class = b),
  prior(exponential(1), class = sd)
)

cat("Fitting Bayesian Paternal Fecundity Model (mod_dad1)...\n")
mod_dad1 <- brm(
  Offspring | cens(censored) ~ 
    me(SSB_mu, SSB_sd) + me(DSB_mu, DSB_sd) + offset(log(breeding)) + 
    (1 | SireGen),
  data      = joined_data, 
  family    = poisson(), #Doesn't run with negbin distribution
  prior     = priors_dad,
  chains    = 4, 
  cores     = 4,
  warmup    = 1000, 
  iter      = 4000, 
  control   = list(adapt_delta = 0.95),
  seed      = 42,
  save_pars = save_pars(all = TRUE)
)

summary(mod_dad1)
pp_check(mod_dad1, type = "bars") +  xlim(-1, 70)

# ------------------------------------------------------------------------------
# 6. Publication Visualizations (SSB & DSB vs Father Offspring)
# ------------------------------------------------------------------------------
# Panel A: Son SSB Propensity vs Father Offspring
ce_ssb <- conditional_effects(mod_dad1, "SSB_mu", prob = 0.95)
p_ssb  <- plot(ce_ssb, plot = FALSE, line_args = list(color = "black"))[[1]] +
  labs(
    x = "Son's SSB propensity (log rate)",
    y = "Father's lifetime offspring"
  ) +
  theme_cowplot() +
  theme(
    axis.title = element_text(size = 12, face = "plain", color = "black"),
    axis.text  = element_text(size = 11, face = "plain", color = "black")
  )
p_ssb$layers[[1]]$aes_params$fill  <- "goldenrod2"
p_ssb$layers[[1]]$aes_params$alpha <- 0.25

# Panel B: Son DSB Propensity vs Father Offspring
ce_dsb <- conditional_effects(mod_dad1, "DSB_mu", prob = 0.95)
p_dsb  <- plot(ce_dsb, plot = FALSE, line_args = list(color = "black"))[[1]] +
  labs(
    x = "Son's DSB propensity (log rate)",
    y = "Father's lifetime offspring"
  ) +
  theme_cowplot() +
  theme(
    axis.title = element_text(size = 12, face = "plain", color = "black"),
    axis.text  = element_text(size = 11, face = "plain", color = "black")
  )
p_dsb$layers[[1]]$aes_params$fill  <- "darkorchid2"
p_dsb$layers[[1]]$aes_params$alpha <- 0.25

# Combine panels side-by-side
fig_paternal_pleiotropy <- cowplot::plot_grid(
  p_ssb, 
  p_dsb, 
  nrow           = 1, 
  ncol           = 2, 
  labels         = c("A", "B"),
  label_size     = 14,
  label_fontface = "bold"
)

# Display plot in active session
print(fig_paternal_pleiotropy)

# ------------------------------------------------------------------------------
# 7. Save Figures and Model Outputs
# ------------------------------------------------------------------------------
output_dir <- "/SSB_delayed_fitness/"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# Save PDF & PNG
ggplot2::ggsave(
  filename = file.path(output_dir, "fig_paternal_pleiotropy.pdf"),
  plot     = fig_paternal_pleiotropy,
  width    = 8,
  height   = 4,
  device   = "pdf"
)

cat(sprintf("\nPipeline completed successfully! Output saved to:\n%s\n", output_dir))