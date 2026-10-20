#### Analysis of Dyadic SSB (Age & Rank Differences)
##Louis Bell-Roberts
#9/10/26

library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

#Read in data
base_dir <- "/Volumes/Macintosh HD/Users/pn26026/Documents/Sandisk/Imperial_HPC/Dyadic_analyses/"

master_dyad_df <- readxl::read_excel(path = "/Users/pn26026/Library/CloudStorage/OneDrive-UniversityofBristol/SSB/Data/Supplementary_Tables_Collated_v1.3.xlsx", sheet = 17)

# ------------------------------------------------------------------------------
# 1. Sparsity & Dyad Overview ####
# ------------------------------------------------------------------------------
cat("\n=== 1. NETWORK SPARSITY & COUNT DISTRIBUTION ===\n")

total_dyads <- nrow(master_dyad_df)
ssb_dyads   <- master_dyad_df %>% filter(Total_Count_M > 0)
num_ssb     <- nrow(ssb_dyads)

cat(sprintf("Total available dyads: %d\n", total_dyads))
cat(sprintf("SSB dyads (Total_Count_M > 0): %d (%.2f%% of all dyads)\n", 
            num_ssb, (num_ssb / total_dyads) * 100))

# Distribution of mount counts among active dyads
mount_dist <- table(ssb_dyads$Total_Count_M)
cat("\nDistribution of mount counts among active dyads:\n")
print(mount_dist)

# Dyadic reciprocity (How many dyads mount in both directions?)
reciprocal_dyads <- ssb_dyads %>% 
  filter(ID1_to_ID2_M > 0 & ID2_to_ID1_M > 0) %>% 
  nrow()

cat(sprintf("\nDyads with bidirectional mounting: %d (%.2f%% of active dyads)\n",
            reciprocal_dyads, (reciprocal_dyads / num_ssb) * 100))


# ------------------------------------------------------------------------------
# 2. Dyadic Assortativity (Undirected Differences); SSB Dyads vs. Background Opportunity Pool ####
# ------------------------------------------------------------------------------
master_dyad_df <- master_dyad_df %>%
  mutate(
    SSB_Binary    = ifelse(Total_Count_M > 0, "SSB Dyad", "Non-SSB Dyad"),
    abs_Age_Diff  = abs(Age_Diff_ID1_minus_ID2),
    abs_Rank_Diff = abs(Rank_Diff_ID1_minus_ID2)
  )

# Summary statistics table
assort_summary <- master_dyad_df %>%
  group_by(SSB_Binary) %>%
  summarise(
    N_Dyads          = n(),
    Median_Age_Diff  = median(abs_Age_Diff, na.rm = TRUE),
    IQR_Age_Diff     = IQR(abs_Age_Diff, na.rm = TRUE),
    Mean_Age_Diff    = mean(abs_Age_Diff, na.rm = TRUE),
    SD_Age_Diff      = sd(abs_Age_Diff, na.rm = TRUE),
    Median_Rank_Diff = median(abs_Rank_Diff, na.rm = TRUE),
    IQR_Rank_Diff    = IQR(abs_Rank_Diff, na.rm = TRUE),
    Mean_Rank_Diff   = mean(abs_Rank_Diff, na.rm = TRUE),
    SD_Rank_Diff     = sd(abs_Rank_Diff, na.rm = TRUE),
    .groups = "drop"
  )

cat("\n=== 2. DYADIC ASSORTATIVITY SUMMARY STATISTICS ===\n")
print(as.data.frame(assort_summary))

# Non-parametric rank-sum tests (Observed vs. Background Pool)
test_age  <- wilcox.test(abs_Age_Diff ~ SSB_Binary, data = master_dyad_df)
test_rank <- wilcox.test(abs_Rank_Diff ~ SSB_Binary, data = master_dyad_df)

cat(sprintf("\nWilcoxon test - Absolute Age Diff (SSB vs Non-SSB): W = %.1f, p = %.4e\n", 
            test_age$statistic, test_age$p.value))
cat(sprintf("Wilcoxon test - Absolute Rank Diff (SSB vs Non-SSB): W = %.1f, p = %.4e\n", 
            test_rank$statistic, test_rank$p.value))


# ------------------------------------------------------------------------------
# 3. Directional Asymmetry (Who Mounts Whom in Observed Events) ####
# ------------------------------------------------------------------------------
# Reshape to directed mounting interaction level
directed_mounts <- bind_rows(
  # ID1 as actor, ID2 as recipient
  master_dyad_df %>%
    filter(ID1_to_ID2_M > 0) %>%
    transmute(
      Year, Group, dyad,
      actor_id       = ID1,
      recipient_id   = ID2,
      mount_count    = ID1_to_ID2_M,
      signed_age_diff  = Age_Diff_ID1_minus_ID2,   # Actor Age minus Recipient Age
      signed_rank_diff = Rank_Diff_ID1_minus_ID2   # Actor DS minus Recipient DS
    ),
  # ID2 as actor, ID1 as recipient
  master_dyad_df %>%
    filter(ID2_to_ID1_M > 0) %>%
    transmute(
      Year, Group, dyad,
      actor_id       = ID2,
      recipient_id   = ID1,
      mount_count    = ID2_to_ID1_M,
      signed_age_diff  = Age_Diff_ID2_minus_ID1,   # Actor Age minus Recipient Age
      signed_rank_diff = Rank_Diff_ID2_minus_ID1   # Actor DS minus Recipient DS
    )
)

total_directed_events <- nrow(directed_mounts)

# Percentages
pct_actor_older  <- mean(directed_mounts$signed_age_diff > 0, na.rm = TRUE) * 100
pct_actor_higher <- mean(directed_mounts$signed_rank_diff > 0, na.rm = TRUE) * 100

cat("\n=== 3. DIRECTIONAL ASYMMETRY IN OBSERVED MOUNTS ===\n")
cat(sprintf("Total directed mounting dyad-instances: %d\n", total_directed_events))
cat(sprintf("Actor is OLDER than recipient: %.1f%%\n", pct_actor_older))
cat(sprintf("Actor is HIGHER-RANKING than recipient: %.1f%%\n", pct_actor_higher))

# Binomial exact tests against 50% chance
binom_age  <- binom.test(sum(directed_mounts$signed_age_diff > 0, na.rm = TRUE), 
                         sum(!is.na(directed_mounts$signed_age_diff)), p = 0.5)
binom_rank <- binom.test(sum(directed_mounts$signed_rank_diff > 0, na.rm = TRUE), 
                         sum(!is.na(directed_mounts$signed_rank_diff)), p = 0.5)

cat(sprintf("Binomial test (Actor older = 50%%): p = %.4f\n", binom_age$p.value))
cat(sprintf("Binomial test (Actor higher-ranking = 50%%): p = %.4f\n", binom_rank$p.value))


# ------------------------------------------------------------------------------
# 4. Exploratory Visualizations ####
# ------------------------------------------------------------------------------
theme_pub <- theme_classic(base_size = 12) +
  theme(
    axis.text       = element_text(color = "black"),
    legend.position = "top",
    plot.title      = element_text(face = "bold", size = 11)
  )

palette_colors <- c("Non-SSB Dyad" = "#A0AAB2", "SSB Dyad" = "#D95F02")

# Panel A: Absolute Age Difference (Assortativity)
p_age_assort <- ggplot(master_dyad_df, aes(x = abs_Age_Diff, fill = SSB_Binary, color = SSB_Binary)) +
  geom_density(alpha = 0.45, linewidth = 0.7) +
  scale_fill_manual(values = palette_colors, name = "") +
  scale_color_manual(values = palette_colors, name = "") +
  labs(
    title = "A. Dyadic Age Assortativity",
    x     = expression(paste("Absolute Age Difference (|", Delta, "Age|, years)")),
    y     = "Density"
  ) +
  theme_pub

# Panel B: Absolute Rank Difference (Assortativity)
p_rank_assort <- ggplot(master_dyad_df, aes(x = abs_Rank_Diff, fill = SSB_Binary, color = SSB_Binary)) +
  geom_density(alpha = 0.45, linewidth = 0.7) +
  scale_fill_manual(values = palette_colors, name = "") +
  scale_color_manual(values = palette_colors, name = "") +
  labs(
    title = "B. Dyadic Rank Assortativity",
    x     = expression(paste("Absolute Rank Difference (|", Delta, "David's Score|)")),
    y     = "Density"
  ) +
  theme_pub

# Panel C: Directional Signed Age Difference (Who mounts whom?)
p_dir_age <- ggplot(directed_mounts, aes(x = signed_age_diff)) +
  geom_histogram(binwidth = 1, fill = "#386CB0", color = "white", alpha = 0.85) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.9) +
  annotate("text", x = -Inf, y = Inf, label = "Actor Younger", hjust = -0.1, vjust = 1.5, size = 3.5, fontface = "italic") +
  annotate("text", x = Inf,  y = Inf, label = "Actor Older",   hjust = 1.1,  vjust = 1.5, size = 3.5, fontface = "italic") +
  labs(
    title = "C. Directional Age Disparity",
    x     = expression(paste(Age[Actor] - Age[Recipient], " (years)")),
    y     = "Count of Mounts"
  ) +
  theme_pub

# Panel D: Directional Signed Rank Difference (Who mounts whom?)
p_dir_rank <- ggplot(directed_mounts, aes(x = signed_rank_diff)) +
  geom_histogram(binwidth = 5, fill = "#7FC97F", color = "white", alpha = 0.85) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.9) +
  annotate("text", x = -Inf, y = Inf, label = "Actor Lower Rank", hjust = -0.1, vjust = 1.5, size = 3.5, fontface = "italic") +
  annotate("text", x = Inf,  y = Inf, label = "Actor Higher Rank", hjust = 1.1,  vjust = 1.5, size = 3.5, fontface = "italic") +
  labs(
    title = "D. Directional Rank Disparity",
    x     = expression(paste(Rank[Actor] - Rank[Recipient], " (David's Score)")),
    y     = "Count of Mounts"
  ) +
  theme_pub

# Combine into a 2x2 multi-panel figure
eda_dyad_plot <- (p_age_assort + p_rank_assort) / (p_dir_age + p_dir_rank) +
  plot_annotation(
    title = "Exploratory Dyadic Structure of Male Same-Sex Mounting",
    theme = theme(plot.title = element_text(face = "bold", size = 14))
  )

# Display in active session
print(eda_dyad_plot)

# Save figure to PDF
plot_output_path <- file.path(base_dir, "Dyadic_SSB_Exploratory_Distributions.pdf")
# ggsave(plot_output_path, eda_dyad_plot, width = 9.5, height = 7.5)
cat(sprintf("\nPlot successfully saved to: %s\n", plot_output_path))

# ==============================================================================
# 1. EVALUATE PSEUDOREPLICATION & NODE REUSE ####
# ==============================================================================
library(tibble)

# A. Pseudoreplication among Active SSB Dyads (Total_Count_M > 0)
ssb_df <- master_dyad_df %>% filter(Total_Count_M > 0)

# Extract all male instances in SSB dyads (54 dyads * 2 = 108 node slots)
ssb_nodes_all <- c(ssb_df$ID1, ssb_df$ID2)
unique_ssb_males <- unique(ssb_nodes_all)

# Partner degree distribution (how many SSB dyads each male appears in)
ssb_degree_dist <- table(table(ssb_nodes_all))

# Across-year participation of SSB males
ssb_males_by_year <- bind_rows(
  ssb_df %>% select(Year, ID = ID1),
  ssb_df %>% select(Year, ID = ID2)
) %>% 
  distinct() %>% 
  count(ID, name = "Years_Active_in_SSB")

# B. Summary Table Output
cat("\n=== PSEUDOREPLICATION IN ACTIVE SSB DYADS ===\n")
cat(sprintf("Total active SSB dyads: %d (108 dyad slots)\n", nrow(ssb_df)))
cat(sprintf("Unique individual males involved: %d\n", length(unique_ssb_males)))
cat(sprintf("Node repetition ratio (Total slots / Unique males): %.2f\n\n", 
            length(ssb_nodes_all) / length(unique_ssb_males)))

cat("Distribution of SSB partner count per male:\n")
print(tibble(
  `Number of SSB Partners` = as.numeric(names(ssb_degree_dist)),
  `Number of Males`        = as.numeric(ssb_degree_dist)
))

cat("\nLongitudinal persistence in SSB:\n")
print(table(ssb_males_by_year$Years_Active_in_SSB))

# C. Background Pool Comparison
bg_nodes_all <- c(master_dyad_df$ID1, master_dyad_df$ID2)
cat(sprintf("\nTotal unique males across all 5,944 dyads: %d\n", length(unique(bg_nodes_all))))

# ==============================================================================
# VIF Diagnostic Check for Dyadic GLMM Predictors ####
# ==============================================================================
library(car)

# Fit a fixed-effects logistic regression on the identical dyadic dataframe
vif_glm <- glm(
  SSB_Binary ~ abs_Age_Diff + abs_Rank_Diff + mean_Age + mean_Rank + log_scans,
  data   = glmm_df,
  family = binomial(link = "logit")
)

# Calculate VIF scores
vif_scores <- car::vif(vif_glm)
print(round(vif_scores, 2))

# ==============================================================================
# Dyadic Logistic GLMM: Partner Choice (Age & Rank Disparity)
# Incorporating SRI, Observability Filtering, Offset, and Dyad Random Effect
# ==============================================================================
library(brms)
library(dplyr)
library(posterior)

# ------------------------------------------------------------------------------
# 1. Helper Function for Within-Cohort Scaling
# ------------------------------------------------------------------------------
safe_scale <- function(x) {
  if (all(is.na(x)) || sd(x, na.rm = TRUE) == 0) {
    return(rep(0, length(x)))
  } else {
    return(as.numeric(scale(x)))
  }
}

# ------------------------------------------------------------------------------
# 2. Filtering & Dyadic Feature Engineering
# ------------------------------------------------------------------------------
min_scan_threshold <- 10

glmm_refined_df <- master_dyad_df %>%
  # Apply minimum observability threshold: both males must meet >= 10 scans
  filter(
    ID1_yearly_scans >= min_scan_threshold,
    ID2_yearly_scans >= min_scan_threshold
  ) %>%
  mutate(
    # Binary response (0 = no mount, 1 = >= 1 mount)
    SSB_Binary = as.integer(Total_Count_M > 0),
    
    # Dyadic observation effort (total scans where either or both were observable)
    dyadic_effort     = ID1_yearly_scans + ID2_yearly_scans - Shared_Scans,
    log_dyadic_effort = log(ifelse(dyadic_effort <= 0, 1, dyadic_effort)),
    
    # Simple Ratio Index (SRI)
    SRI = Shared_Scans / ifelse(dyadic_effort <= 0, 1, dyadic_effort),
    
    # Absolute age and rank differences
    abs_age_diff  = abs(Age_Diff_ID1_minus_ID2),
    abs_rank_diff = abs(Rank_Diff_ID1_minus_ID2),
    
    # Dyad baseline means
    mean_age  = (Age_ID1 + Age_ID2) / 2,
    mean_rank = (DS_ID1 + DS_ID2) / 2,
    
    # Grouping variable for cohort intercepts
    Group_Year = paste(Group, Year, sep = "_"),
    
    # Factor conversions required for random effects
    ID1  = as.factor(ID1),
    ID2  = as.factor(ID2),
    dyad = as.factor(dyad)
  ) %>%
  filter(dyadic_effort > 0 & !is.nan(SRI)) %>%
  # Scale all continuous focal predictors (z-scored on 1-SD scale)
  mutate(
    abs_age_diff_s  = safe_scale(abs_age_diff),
    abs_rank_diff_s = safe_scale(abs_rank_diff),
    mean_age_s      = safe_scale(mean_age),
    mean_rank_s     = safe_scale(mean_rank),
    SRI_s           = safe_scale(SRI)
  )

cat(sprintf("Dyads remaining after scan filtering (>= %d scans): %d\n", 
            min_scan_threshold, nrow(glmm_refined_df)))
cat(sprintf("Active SSB dyads in filtered dataset: %d (%.2f%%)\n", 
            sum(glmm_refined_df$SSB_Binary), 
            mean(glmm_refined_df$SSB_Binary) * 100))

# ------------------------------------------------------------------------------
# 3. Model Formula Definition
# ------------------------------------------------------------------------------
bf_dyad_final <- bf(
  SSB_Binary ~ abs_age_diff_s + abs_rank_diff_s + 
    mean_age_s + mean_rank_s + 
    SRI_s + offset(log_dyadic_effort) + 
    (1 | mm(ID1, ID2)) + 
    #(1 | dyad) + dropped this random effect due to presence of divergent transitions
    (1 | Group_Year),
  family = bernoulli(link = "logit")
)

# ------------------------------------------------------------------------------
# 4. Priors Tailored to Standardized Scales & Unscaled Log Offset
# ------------------------------------------------------------------------------
priors_final <- c(
  prior(normal(-9, 2),     class = "Intercept"),
  prior(normal(0, 0.5),    class = "b"),
  prior(exponential(2),    class = "sd")
)

# ------------------------------------------------------------------------------
# 5. Fit Bayesian GLMM
# ------------------------------------------------------------------------------
cat("\nFitting Bayesian Dyadic Logistic GLMM...\n")

fit_dyad_final <- brm(
  formula   = bf_dyad_final,
  data      = glmm_refined_df,
  prior     = priors_final,
  chains    = 4,
  cores     = 4,
  warmup    = 1500,
  iter      = 3500,
  control   = list(adapt_delta = 0.99, max_treedepth = 15),
  seed      = 42,
  save_pars = save_pars(all = TRUE)
)

# ------------------------------------------------------------------------------
# 6. Model Output & Posterior Probabilities
# ------------------------------------------------------------------------------
summary(fit_dyad_final)

post_draws_final <- as_draws_df(fit_dyad_final)

# Posterior probabilities of directional hypotheses
prob_age_homophily  <- mean(post_draws_final$b_abs_age_diff_s < 0)
prob_rank_homophily <- mean(post_draws_final$b_abs_rank_diff_s < 0)

cat("\n=== HYPOTHESIS TESTING RESULTS ===\n")
cat(sprintf("Posterior P(b_abs_age_diff_s < 0)  [Age Homophily]:  %.4f\n", prob_age_homophily))
cat(sprintf("Posterior P(b_abs_rank_diff_s < 0) [Rank Homophily]: %.4f\n", prob_rank_homophily))

# saveRDS(fit_dyad_final, file = file.path(results_dir, "fit_dyad_logistic_brms.rds"))

# ==============================================================================
# Reduced Marginal Model for Age Homophily
# (Omits dyadic baseline means to estimate marginal partner assortativity)
# ==============================================================================
library(brms)

# 1. Define Formula
bf_dyad_reduced <- bf(
  SSB_Binary ~ abs_age_diff_s + abs_rank_diff_s + 
    SRI_s + offset(log_dyadic_effort) + 
    (1 | mm(ID1, ID2)) + 
    (1 | Group_Year),
  family = bernoulli(link = "logit")
)

# 2. Define Priors
priors_reduced <- c(
  prior(normal(-9, 2),     class = "Intercept"),
  prior(normal(0, 0.5),    class = "b"),
  prior(exponential(2),    class = "sd")
)

# 3. Fit Reduced Model
fit_dyad_reduced <- brm(
  formula   = bf_dyad_reduced,
  data      = glmm_refined_df,
  prior     = priors_reduced,
  chains    = 4, 
  cores     = 4, 
  warmup    = 1500, 
  iter      = 3500,
  control   = list(adapt_delta = 0.99, max_treedepth = 15),
  seed      = 42,
  save_pars = save_pars(all = TRUE)
)

# 4. View Results
summary(fit_dyad_reduced)

# 5. Extract Directional Probabilities
post_draws_reduced <- as_draws_df(fit_dyad_reduced)
cat(sprintf("\nPosterior P(beta_abs_age_diff_s < 0) [Age Homophily]: %.4f\n", 
            mean(post_draws_reduced$b_abs_age_diff_s < 0)))
cat(sprintf("Posterior P(beta_abs_rank_diff_s < 0) [Rank Homophily]: %.4f\n", 
            mean(post_draws_reduced$b_abs_rank_diff_s < 0)))

###
#Create output tables ####
###

library(brms)
library(dplyr)
library(tibble)
library(readr)
library(posterior)

extract_glmm_summary <- function(model, model_name = "Model") {
  
  # 1. Fixed Effects
  fe <- fixef(model) %>%
    as.data.frame() %>%
    rownames_to_column(var = "Parameter") %>%
    mutate(
      OR = ifelse(Parameter == "Intercept", NA_real_, exp(Estimate)),
      sig = ifelse((`Q2.5` > 0 & `Q97.5` > 0) | (`Q2.5` < 0 & `Q97.5` < 0), "*", ""),
      Est_CI = sprintf("%.2f [%.2f, %.2f]%s", Estimate, `Q2.5`, `Q97.5`, sig),
      OR_str = ifelse(is.na(OR), "—", sprintf("%.2f", OR)),
      Type = "Fixed Effect"
    ) %>%
    select(Type, Parameter, Est_CI, OR_str)
  
  # 2. Random Effects (Standard Deviations)
  vc <- VarCorr(model)
  re_list <- list()
  
  for (grp in names(vc)) {
    sd_est <- vc[[grp]]$sd[1, "Estimate"]
    sd_lci <- vc[[grp]]$sd[1, "Q2.5"]
    sd_uci <- vc[[grp]]$sd[1, "Q97.5"]
    
    re_list[[grp]] <- tibble(
      Type      = "Random Effect (SD)",
      Parameter = paste0("sd_", grp),
      Est_CI    = sprintf("%.2f [%.2f, %.2f]", sd_est, sd_lci, sd_uci),
      OR_str    = "—"
    )
  }
  re_df <- bind_rows(re_list)
  
  # 3. MCMC Diagnostics (Compatible with both rstan and cmdstanr)
  mod_sum     <- summary(model)
  nuts_d      <- nuts_params(model)
  n_divergent <- sum(subset(nuts_d, Parameter == "divergent__")$Value)
  rhat_vals   <- brms::rhat(model)
  max_rhat    <- max(rhat_vals, na.rm = TRUE)
  
  bulk_ess_fe  <- mod_sum$fixed$Bulk_ESS
  bulk_ess_re  <- unlist(lapply(mod_sum$random, function(x) x$Bulk_ESS))
  min_bulk_ess <- round(min(c(bulk_ess_fe, bulk_ess_re), na.rm = TRUE))
  
  iter_val   <- ifelse(!is.null(mod_sum$iter), mod_sum$iter, 3500)
  warmup_val <- ifelse(!is.null(mod_sum$warmup), mod_sum$warmup, 1500)
  
  diag_df <- tibble(
    Type      = "Diagnostic",
    Parameter = c("Iterations / Warmup", "Post-warmup Draws", 
                  "Divergent Transitions", "Maximum Rhat", "Minimum Bulk-ESS"),
    Est_CI    = c(
      sprintf("%s / %s", iter_val, warmup_val),
      as.character(ndraws(model)),
      as.character(n_divergent),
      sprintf("%.2f", max_rhat),
      as.character(min_bulk_ess)
    ),
    OR_str    = "—"
  )
  
  out <- bind_rows(fe, re_df, diag_df)
  colnames(out) <- c("Type", "Parameter", paste0(model_name, "_Est_CI"), paste0(model_name, "_OR"))
  return(out)
}

# Extract summaries dynamically
m1_summary <- extract_glmm_summary(fit_dyad_reduced, model_name = "Model1")
m2_summary <- extract_glmm_summary(fit_dyad_final,   model_name = "Model2")

# Master parameter list for display ordering
param_order <- tibble(
  Parameter = c(
    "Intercept",
    "abs_age_diff_s",
    "abs_rank_diff_s",
    "SRI_s",
    "mean_age_s",
    "mean_rank_s",
    "sd_mmID1ID2",
    "sd_Group_Year",
    "Iterations / Warmup",
    "Post-warmup Draws",
    "Divergent Transitions",
    "Maximum Rhat",
    "Minimum Bulk-ESS"
  ),
  Display_Name = c(
    "Intercept",
    "Age Disparity (|delta Age|, scaled)",
    "Rank Disparity (|delta Rank|, scaled)",
    "Spatial Association (SRI, scaled)",
    "Dyad Mean Age (scaled)",
    "Dyad Mean Rank (scaled)",
    "Male Identity (mmID1ID2, 142 levels)",
    "Group-Year (6 levels)",
    "Iterations / Warmup",
    "Post-warmup Draws",
    "Divergent Transitions",
    "Maximum Rhat",
    "Minimum Bulk-ESS"
  ),
  Component = c(
    rep("Fixed Effect", 6),
    rep("Random Effect (SD)", 2),
    rep("Diagnostic", 5)
  )
)

# Merge extracted estimates into publication table
supp_table_dyadic <- param_order %>%
  left_join(m1_summary %>% select(-Type), by = "Parameter") %>%
  left_join(m2_summary %>% select(-Type), by = "Parameter") %>%
  mutate(
    Model1_Est_CI = ifelse(is.na(Model1_Est_CI), "—", Model1_Est_CI),
    Model1_OR     = ifelse(is.na(Model1_OR), "—", Model1_OR),
    Model2_Est_CI = ifelse(is.na(Model2_Est_CI), "—", Model2_Est_CI),
    Model2_OR     = ifelse(is.na(Model2_OR), "—", Model2_OR)
  ) %>%
  select(
    Component,
    Predictor = Display_Name,
    Model1_Estimate_95CI = Model1_Est_CI,
    Model1_Odds_Ratio    = Model1_OR,
    Model2_Estimate_95CI = Model2_Est_CI,
    Model2_Odds_Ratio    = Model2_OR
  )

# Display in console
print(as.data.frame(supp_table_dyadic))

# Export to CSV
output_path <- file.path("/Users/pn26026/OneDrive - University of Bristol/SSB/Code/Supp_Table_Dyadic_age_rank_GLMM.csv")
write_csv(supp_table_dyadic, output_path)
