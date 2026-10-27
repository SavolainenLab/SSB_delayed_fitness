####### Functional Benefits of SSB: Sexually Antagonistic Pleiotropy Analysis #######
## Author: Louis Bell-Roberts
## Date: October 2026
## Description: Tests whether male same-sex sexual behaviour (SSB) expression is 
## associated with sister lifetime reproductive success. Fits bivariate  
## models to estimate son behavioural propensities and right-censored mixed 
## models to evaluate sister fecundity.
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
data_path <- "/SSB_delayed_fitness/sister_pleiotropy_data.csv"

data <- readr::read_csv(data_path, show_col_types = FALSE) %>%
  dplyr::filter(Scan_Count >= 10, sex == "M", !is.na(Group)) %>%
  dplyr::mutate(Animal_ID = as.factor(Animal_ID))

# Partition age into within- and between-individual components (pure numeric scaling)
data_biv <- data %>%
  dplyr::distinct(Animal_ID, Year, .keep_all = TRUE) %>% # Ensure bivariate model fits on unique focal male annual observations
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
# 2. Fit Bivariate Model (Brother SSB & DSB Propensities)
# ------------------------------------------------------------------------------
bf_SSB <- bf(
  SSB_Raw ~ 1 + age_B_s + age_W_s + offset(log(Scan_Count)) + 
    (1 | p | Animal_ID) + (1 | Group_Year),
  family = negbinomial()
)

bf_DSB <- bf(
  DSB_Raw ~ 1 + age_B_s + age_W_s + offset(log(Scan_Count)) + 
    (1 | p | Animal_ID) + (1 | Group_Year),
  family = negbinomial()
)

cat("Fitting Bivariate Model for Brother SSB and DSB...\n")
biv_w_age_e2 <- brm(
  bf_SSB + bf_DSB + set_rescor(FALSE),
  data      = data_biv, 
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
# 3. Extract Brother Random Intercepts (BLUPs with Standard Errors)
# ------------------------------------------------------------------------------
ran_vals <- ranef(biv_w_age_e2)$Animal_ID

# Robust column matching (handles both 'SSBRaw_Intercept' and 'SSB_Raw_Intercept')
ssb_col <- grep("SSB.*Intercept", dimnames(ran_vals)[[3]], value = TRUE)[1]
dsb_col <- grep("DSB.*Intercept", dimnames(ran_vals)[[3]], value = TRUE)[1]

brothers <- tibble::tibble(
  Animal_ID = as.character(rownames(ran_vals)),
  SSB_mu    = ran_vals[, "Estimate", ssb_col],
  SSB_sd    = ran_vals[, "Est.Error", ssb_col],
  DSB_mu    = ran_vals[, "Estimate", dsb_col],
  DSB_sd    = ran_vals[, "Est.Error", dsb_col]
)

# ------------------------------------------------------------------------------
# 4. Prepare Sister Lifespan & Right-Censoring
# ------------------------------------------------------------------------------
max_age        <- 24
min_age        <- 3
study_end_date <- lubridate::ymd("2025-03-01")

new <- data %>%
  dplyr::mutate(Animal_ID = as.character(Animal_ID)) %>%
  dplyr::distinct(Animal_ID, Sister_ID, .keep_all = TRUE) %>%
  dplyr::filter(!is.na(Sister_ID)) %>%
  dplyr::mutate(
    CPRC_DOB        = as.Date(Sister_DOB),
    CPRC_death_date = as.Date(Sister_death_date),
    max_date = dplyr::if_else(
      !is.na(CPRC_death_date) & CPRC_death_date < study_end_date,
      CPRC_death_date,
      study_end_date
    ),
    sister_age_last = as.numeric(difftime(max_date, CPRC_DOB, units = "days")) / 365.2425,
    breeding        = pmin(sister_age_last, max_age) - min_age,
    breeding        = pmax(breeding, 0),
    Offspring       = as.integer(Sister_num_offspring),
    censored        = dplyr::if_else(breeding < 21 & max_date == study_end_date, 1, 0),
    DamGen          = as.factor(DamGen),
    Sister_ID       = as.factor(Sister_ID)
  ) %>%
  dplyr::filter(!is.na(breeding), breeding > 0, !is.na(Offspring)) %>%
  dplyr::select(Animal_ID, DamGen, Sister_ID, breeding, Offspring, censored)

# ------------------------------------------------------------------------------
# 5. Join Brother Estimates with Sister Demographics & Fit Fecundity Model
# ------------------------------------------------------------------------------
joined_data <- dplyr::inner_join(brothers, new, by = "Animal_ID")

priors_sis <- c(
  prior(normal(0, 5),   class = Intercept),
  prior(normal(0, 0.5), class = b),
  prior(exponential(1), class = sd)
)

cat("Fitting Bayesian Sister Fecundity Model (mod_sis1)...\n")
mod_sis1 <- brm(
  Offspring | cens(censored) ~ 
    me(SSB_mu, SSB_sd) + me(DSB_mu, DSB_sd) + offset(log(breeding)) + 
    (1 | DamGen),
  data      = joined_data, 
  family    = poisson(), #poisson better than negbin
  prior     = priors_sis,
  chains    = 4, 
  cores     = 4,
  warmup    = 1000, 
  iter      = 4000, 
  control   = list(adapt_delta = 0.95),
  seed      = 42,
  save_pars = save_pars(all = TRUE)
)

pp_check(mod_sis1, type = "bars")
summary(mod_sis1)

#Model including sister ID random effect
mod_sis2 <- brm(
  Offspring | cens(censored) ~ 
    me(SSB_mu, SSB_sd) + me(DSB_mu, DSB_sd) + offset(log(breeding)) + 
    (1 | DamGen) + (1 | Sister_ID),
  data      = joined_data, 
  family    = poisson(),
  prior     = priors_sis,
  chains    = 4, 
  cores     = 4,
  warmup    = 1000, 
  iter      = 4000, 
  control   = list(adapt_delta = 0.95),
  seed      = 42,
  save_pars = save_pars(all = TRUE)
)

pp_check(mod_sis2, type = "bars")
summary(mod_sis2) #Selecting this model structure for consistency across the 4 different familial relationships

# ------------------------------------------------------------------------------
# 6. Publication Visualizations (SSB & DSB vs Sister Offspring)
# ------------------------------------------------------------------------------
# Panel A: Brother SSB Propensity vs Sister Offspring
ce_ssb <- conditional_effects(mod_sis2, "SSB_mu", prob = 0.95)
p_ssb  <- plot(ce_ssb, plot = FALSE, line_args = list(color = "black"))[[1]] +
  labs(
    x = "Brother's SSB propensity (log rate)",
    y = "Sister's lifetime offspring"
  ) +
  theme_cowplot() +
  theme(
    axis.title = element_text(size = 12, face = "plain", color = "black"),
    axis.text  = element_text(size = 11, face = "plain", color = "black")
  )
p_ssb$layers[[1]]$aes_params$fill  <- "goldenrod2"
p_ssb$layers[[1]]$aes_params$alpha <- 0.25

# Panel B: Brother DSB Propensity vs Sister Offspring
ce_dsb <- conditional_effects(mod_sis2, "DSB_mu", prob = 0.95)
p_dsb  <- plot(ce_dsb, plot = FALSE, line_args = list(color = "black"))[[1]] +
  labs(
    x = "Brother's DSB propensity (log rate)",
    y = "Sister's lifetime offspring"
  ) +
  theme_cowplot() +
  theme(
    axis.title = element_text(size = 12, face = "plain", color = "black"),
    axis.text  = element_text(size = 11, face = "plain", color = "black")
  )
p_dsb$layers[[1]]$aes_params$fill  <- "darkorchid2"
p_dsb$layers[[1]]$aes_params$alpha <- 0.25

# Combine panels side-by-side
fig_sister_pleiotropy <- cowplot::plot_grid(
  p_ssb, 
  p_dsb, 
  nrow           = 1, 
  ncol           = 2, 
  labels         = c("A", "B"),
  label_size     = 14,
  label_fontface = "bold"
)

# Display plot in active session
print(fig_sister_pleiotropy)

# ------------------------------------------------------------------------------
# 7. Save Figures and Model Outputs
# ------------------------------------------------------------------------------
output_dir <- "/SSB_delayed_fitness/"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# Save PDF & PNG
ggplot2::ggsave(
  filename = file.path(output_dir, "fig_sister_pleiotropy.pdf"),
  plot     = fig_sister_pleiotropy,
  width    = 8,
  height   = 4,
  device   = "pdf"
)

cat(sprintf("\nPipeline completed successfully! Output saved to:\n%s\n", output_dir))