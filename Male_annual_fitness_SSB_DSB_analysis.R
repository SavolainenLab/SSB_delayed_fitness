#### Testing the annual association between SSB and number of offspring sired ####
##Louis Bell-Roberts
#12/05/2026

library(tidyverse)
library(brms)
library(ggstance)
library(ggplot2)
library(coda)
library(bayesplot)
library(tidybayes)
library(patchwork)
library(lubridate)
library(broom.mixed)
library(posterior)
library(cowplot)

#Read in data

# File path (Local working path; convert to relative path prior to final archive)
supp_data_path <- "/SSB_delayed_fitness/Supplementary_Tables_Collated_v1.4.xlsx"

# Sheet 18: Combined cleaned collapsed data with observation effort
# Sheet 19: Cleaned pedigree data
data <- readxl::read_excel(path = supp_data_path, sheet = 19)
pedigree_pre <- readxl::read_excel(path = supp_data_path, sheet = 20)

# --- Create Analysis DF ---
analysis_df <- data %>% 
  filter(
    sex == "M",
    Scan_Count >= 10,
    !is.na(assigned_group)
  ) %>% 
  group_by(Animal_ID) %>% 
  mutate( #Calculate age variables
    age_B = mean(Age),
    age_W = Age - age_B
  ) %>%
  ungroup() %>%
  mutate( #Scale age variables
    Animal_ID = as.factor(Animal_ID),
    Age_s = scale(Age),
    age_B_s = scale(age_B),
    age_W_s = scale(age_W),
    Group_Year = as.factor(paste(assigned_group, Year, sep = "_")),
    adjusted_SSB = (SSB_Raw)/log(Scan_Count),
    adjusted_DSB = (DSB_Raw)/log(Scan_Count)
  )

########################
#Multivariate model####
########################

priors_SSB_DSB_offsp <- c(
  # --- FIXED EFFECTS ---
  prior(normal(0, 5),   class = Intercept, resp = "SSBRaw"),
  prior(normal(0, 5),   class = Intercept, resp = "DSBRaw"),
  prior(normal(0, 5),   class = Intercept, resp = "numoffspringallbirths"),
  prior(normal(0, 0.5), class = b,         resp = "SSBRaw"),
  prior(normal(0, 0.5), class = b,         resp = "DSBRaw"),
  prior(normal(0, 0.5), class = b,         resp = "numoffspringallbirths"),
  
  # --- RANDOM EFFECTS ---
  # 1. SD for the INTERCEPTS
  prior(exponential(1), class = sd, coef = "Intercept", group = "Animal_ID", resp = "SSBRaw"),
  prior(exponential(1), class = sd, coef = "Intercept", group = "Animal_ID", resp = "DSBRaw"),
  prior(exponential(1), class = sd, coef = "Intercept", group = "Animal_ID", resp = "numoffspringallbirths"),
  prior(lkj(2), class = cor, group = "Animal_ID"),
  
  # assigned_group:Year (Intercepts only)
  prior(exponential(1), class = sd, group = "assigned_group:Year", resp = "SSBRaw"),
  prior(exponential(1), class = sd, group = "assigned_group:Year", resp = "DSBRaw"),
  prior(exponential(1), class = sd, group = "assigned_group:Year", resp = "numoffspringallbirths"),
  
  # --- DISTRIBUTION PARAMETERS ---
  prior(exponential(1), class = shape, resp = "SSBRaw"),
  prior(exponential(1), class = shape, resp = "DSBRaw"),
  prior(exponential(1), class = shape, resp = "numoffspringallbirths")
)

bf_SSB <- bf(SSB_Raw ~ 1 + age_W_s + age_B_s + offset(log(Scan_Count)) + (1|p|Animal_ID) + (1|assigned_group:Year)) + negbinomial()

bf_DSB <- bf(DSB_Raw ~ 1 + age_W_s + age_B_s + offset(log(Scan_Count)) + (1|p|Animal_ID) + (1|assigned_group:Year)) + negbinomial()

bf_offspring <- bf(num_offspring_all_births ~ 1 + age_W_s + age_B_s + (1|p|Animal_ID) + (1|assigned_group:Year)) + negbinomial()

biv_SSB_DSB_off <- brm(
  bf_SSB + bf_DSB + bf_offspring + set_rescor(F),
  prior = priors_SSB_DSB_offsp,
  data = analysis_df, seed = 1,
  warmup = 1000, iter = 4000, thin = 1, chains = 4, cores = 4,
  control = list(adapt_delta = 0.90)
)

summary(biv_SSB_DSB_off)

#Calculate pd for the covariance estimates
# Extract all post-warmup draws from the brms model into a data frame
posterior_draws <- as_draws_df(biv_SSB_DSB_off)

# Calculate the Probability of Direction (pd) for the three Animal_ID correlations.
# The formula is simple: what proportion of the 12,000 draws fall below zero (if the mean is negative) or above zero (if the mean is positive)?
pd_correlations <- data.frame(
  Correlation = c("SSB <-> DSB", "SSB <-> Offspring", "DSB <-> Offspring"),
  
  # 1. SSB <-> DSB (Mean is -ve, so we calculate the proportion of draws < 0)
  Pd = c(
    mean(posterior_draws$`cor_Animal_ID__SSBRaw_Intercept__DSBRaw_Intercept` < 0),
    
    # 2. SSB <-> Offspring (Mean is -ve, so we calculate the proportion of draws < 0)
    mean(posterior_draws$`cor_Animal_ID__SSBRaw_Intercept__numoffspringallbirths_Intercept` < 0),
    
    # 3. DSB <-> Offspring (Mean is +ve, so we calculate the proportion of draws > 0)
    mean(posterior_draws$`cor_Animal_ID__DSBRaw_Intercept__numoffspringallbirths_Intercept` > 0)
  )
)

# Convert to percentage and round
pd_correlations$Pd_Percentage <- paste0(round(pd_correlations$Pd * 100, 2), "%")

# Print the results
print(pd_correlations)

pp_check(biv_SSB_DSB_off, resp='SSBRaw', type = "bars") + xlim(c(-1, 60))
pp_check(biv_SSB_DSB_off, resp='DSBRaw', type = "bars") + xlim(c(-1, 50))
pp_check(biv_SSB_DSB_off, resp='numoffspringallbirths', type = "bars") + xlim(c(-1, 40))

#Generate model output tables
# --- 1. Extract brms Summary ---
# The summary contains ESS and Rhat values
mod_sum <- summary(biv_SSB_DSB_off)

# --- 2. Extract and Format Fixed Effects ---
# Convert the fixed effects matrix to a data frame and add metadata columns
fixed_df <- as.data.frame(mod_sum$fixed) %>%
  tibble::rownames_to_column("Raw_Parameter") %>%
  dplyr::mutate(
    Effect_Type = "Fixed",
    Group = NA_character_ # Fixed effects don't have a grouping variable
  )

# --- 3. Extract and Format Random Effects ---
# Loop through all grouping variables (Animal_ID, assigned_group:Year) 
# and extract their respective parameter matrices (SDs and Correlations)
ran_df_list <- lapply(names(mod_sum$random), function(grp) {
  as.data.frame(mod_sum$random[[grp]]) %>%
    tibble::rownames_to_column("Raw_Parameter") %>%
    dplyr::mutate(
      Effect_Type = "Random",
      Group = grp
    )
})
# Bind all random effect groups together into one data frame
ran_df <- dplyr::bind_rows(ran_df_list)

# --- 4. Combine and Select Columns (Including ESS and Rhat) ---
# Merge fixed and random effects, pulling the convergence diagnostics
raw_combined <- dplyr::bind_rows(fixed_df, ran_df) %>%
  dplyr::select(
    Effect_Type,
    Group,
    Raw_Parameter,
    Estimate,
    Est_Error = `Est.Error`,       # Native brms column name
    Lower_95_CI = `l-95% CI`,      # Native brms column name
    Upper_95_CI = `u-95% CI`,      # Native brms column name
    ESS_Bulk = Bulk_ESS,           # Extract Bulk Effective Sample Size
    ESS_Tail = Tail_ESS,           # Extract Tail Effective Sample Size
    Rhat                           # Extract Rhat convergence diagnostic
  )

# --- 5. Clean, Label, and Format the Single Combined Table ---
combined_table <- raw_combined %>%
  dplyr::mutate(
    # Identify which trait the parameter belongs to, flagging correlations specifically
    Response_Variable = dplyr::case_when(
      str_detect(Raw_Parameter, "^cor\\(") ~ "Multivariate (Correlations)",
      str_detect(Raw_Parameter, "SSBRaw") ~ "SSB",
      str_detect(Raw_Parameter, "DSBRaw") ~ "DSB",
      str_detect(Raw_Parameter, "numoffspring") ~ "Offspring",
      TRUE ~ "Other"
    ),
    
    # Translate syntax into readable manuscript labels
    Parameter = Raw_Parameter,
    
    # 1. Clean random effect correlation wrappers first
    Parameter = str_replace_all(Parameter, "^cor\\((.*)\\)$", "Cor: \\1"),
    # Replace the comma separator in correlation terms with a bidirectional arrow
    Parameter = ifelse(str_detect(Parameter, "^Cor: "), str_replace_all(Parameter, ",", " <-> "), Parameter),
    
    # 2. Clean random effect Standard Deviation wrappers
    Parameter = str_replace_all(Parameter, "^sd\\((.*)\\)$", "SD: \\1"), 
    
    # 3. Carefully strip trait prefixes to avoid blanking out correlation names
    Parameter = ifelse(
      Response_Variable == "Multivariate (Correlations)",
      str_replace_all(Parameter, "SSBRaw_", "SSB ") %>%
        str_replace_all("DSBRaw_", "DSB ") %>%
        str_replace_all("numoffspringallbirths_", "Offspring "),
      str_replace_all(Parameter, "(SSBRaw_|DSBRaw_|numoffspringallbirths_)", "")
    ),
    
    # 4. Clean specific predictor labels (including the squared terms for the offspring model)
    Parameter = str_replace_all(Parameter, "Iage_W_sE2", "Within-individual age^2 (scaled)"),
    Parameter = str_replace_all(Parameter, "Iage_B_sE2", "Between-individual age^2 (scaled)"),
    Parameter = str_replace_all(Parameter, "age_W_s", "Within-individual age (scaled)"),
    Parameter = str_replace_all(Parameter, "age_B_s", "Between-individual age (scaled)")
  ) %>%
  # Reorder columns for the final table
  dplyr::select(
    Effect_Type, Response_Variable, Group, Parameter, 
    Estimate, Est_Error, Lower_95_CI, Upper_95_CI, 
    ESS_Bulk, ESS_Tail, Rhat
  ) %>%
  # Round numeric values for consistency
  dplyr::mutate(
    # Round estimates, errors, CIs, and Rhat to 3 decimal places
    dplyr::across(c(Estimate, Est_Error, Lower_95_CI, Upper_95_CI, Rhat), ~ round(.x, 3)),
    # Round ESS to whole numbers since they represent sample sizes
    dplyr::across(c(ESS_Bulk, ESS_Tail), ~ round(.x, 0)) 
  )

# --- 6. Export the Final Unified Table ---
# Save the combined dataframe as a single CSV file
# write.csv(combined_table, "/SSB_delayed_fitness/Output_table_SSB_DSB_offspring.csv", row.names = FALSE)

# Plotting ####

# ==============================================================================
# TWO-PANEL AGE TRAJECTORY FIGURE (SSB & DSB)
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. Prediction Grid (Age: 3.5 to 25)
# ------------------------------------------------------------------------------
# Directly use the mean and SD from your data
mean_age <- mean(analysis_df$Age, na.rm = TRUE)
sd_age_w <- sd(analysis_df$age_W, na.rm = TRUE)

age_grid_base <- data.frame(
  Age = seq(3.5, 25, length.out = 100)
) %>%
  mutate(
    # Project ages onto the within-individual standardised scale
    age_W_s    = (Age - mean_age) / sd_age_w,
    age_B_s    = 0,   # Population average cohort
    Scan_Count = 100  # Baseline exposure for rate standardisation
  )

# ------------------------------------------------------------------------------
# 2. Spaghetti Plotting Function
# ------------------------------------------------------------------------------
create_age_spaghetti_panel <- function(model, 
                                       grid_data, 
                                       y_resp, 
                                       y_label, 
                                       line_color, 
                                       y_limits = c(0, 15)) {
  
  # A. Generate 300 posterior draws
  post_draws <- posterior_epred(
    model,
    newdata    = grid_data,
    resp       = y_resp,
    ndraws     = 300,
    re_formula = NA
  )
  
  # B. Population mean trajectory
  grid_data$estimate__ <- colMeans(post_draws)
  
  # C. Reshape draws for ggplot
  draws_df <- as.data.frame(t(post_draws))
  colnames(draws_df) <- paste0("draw_", seq_len(ncol(draws_df)))
  
  draws_long <- cbind(grid_data, draws_df) %>%
    pivot_longer(
      cols      = starts_with("draw_"),
      names_to  = "draw",
      values_to = "prediction"
    )
  
  # D. Plot
  ggplot() +
    geom_line(
      data      = draws_long,
      aes(x = Age, y = prediction, group = draw),
      color     = line_color,
      alpha     = 0.10,
      linewidth = 0.30
    ) +
    geom_line(
      data      = grid_data,
      aes(x = Age, y = estimate__),
      color     = "black",
      linewidth = 0.8
    ) +
    labs(x = "Age (years)", y = y_label) +
    theme_cowplot() +
    theme(
      panel.background = element_rect(fill = "white", color = NA),
      plot.background  = element_rect(fill = "white", color = NA),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line        = element_line(color = "black", linewidth = 0.6),
      axis.title       = element_text(size = 13, face = "bold"),
      axis.text        = element_text(size = 11, color = "black")
    ) +
    # Ticks start cleanly at 5, but coordinate limits start at 3.5
    scale_x_continuous(
      breaks = seq(5, 25, by = 5),
      expand = c(0, 0)
    ) +
    coord_cartesian(
      xlim = c(3.5, 25), 
      ylim = y_limits
    )
}

# ------------------------------------------------------------------------------
# 3. Render Panels
# ------------------------------------------------------------------------------
# Panel A: SSB (Yellow)
panel_a_ssb <- create_age_spaghetti_panel(
  model      = biv_SSB_DSB_off,
  grid_data  = age_grid_base,
  y_resp     = "SSBRaw",
  y_label    = "Annual SSB mounts",
  line_color = "goldenrod2",
  y_limits   = c(0, 25)
)

# Panel B: DSB (Purple)
panel_b_dsb <- create_age_spaghetti_panel(
  model      = biv_SSB_DSB_off,
  grid_data  = age_grid_base,
  y_resp     = "DSBRaw",
  y_label    = "Annual DSB mounts",
  line_color = "darkorchid2",
  y_limits   = c(0, 60)
)

# ------------------------------------------------------------------------------
# 4. Multi-Panel Composite
# ------------------------------------------------------------------------------
two_panel_age_plot <- plot_grid(
  panel_a_ssb,
  panel_b_dsb,
  nrow           = 1,
  ncol           = 2,
  labels         = c("A", "B"),
  label_size     = 14,
  label_fontface = "bold",
  align          = "h",
  axis           = "tb"
)

print(two_panel_age_plot)

# ==============================================================================
# SECTION: RANDOM INTERCEPTS PLOTS (PANELS C, D, E)
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. Output Directory Configuration & Aesthetic Settings
# ------------------------------------------------------------------------------
output_plots_dir <- "/SSB_delayed_fitness/"
if (!dir.exists(output_plots_dir)) dir.create(output_plots_dir, recursive = TRUE)

# Global panel label font size
panel_label_size <- 28

# ------------------------------------------------------------------------------
# 2. Extract Raw Animal_ID Random Intercepts
# ------------------------------------------------------------------------------
ran_vals <- ranef(biv_SSB_DSB_off)$Animal_ID

intercepts_df <- tibble(
  Animal_ID     = rownames(ran_vals),
  SSB_std       = as.numeric(ran_vals[, "Estimate", "SSBRaw_Intercept"]),
  DSB_std       = as.numeric(ran_vals[, "Estimate", "DSBRaw_Intercept"]),
  Offspring_std = as.numeric(ran_vals[, "Estimate", "numoffspringallbirths_Intercept"])
)

# Extract posterior draws for the correlation parameters
posterior_draws <- as_draws_df(biv_SSB_DSB_off)

# ------------------------------------------------------------------------------
# 3. Helper Functions: Inset Density & Scatterplot
# ------------------------------------------------------------------------------
create_density_inset <- function(draws_vector, fill_color) {
  df <- data.frame(val = draws_vector)
  
  # Calculate 95% Credible Interval and Posterior Median
  ci_95  <- quantile(draws_vector, probs = c(0.025, 0.975))
  pt_est <- median(draws_vector)
  
  # Peak density for coordinate scaling
  dens     <- density(draws_vector)
  max_dens <- max(dens$y)
  
  ggplot(df, aes(x = val)) +
    # Solid density fill with no black outline
    geom_density(
      fill  = fill_color, 
      alpha = 1, 
      color = NA
    ) +
    # Vertical reference line at zero
    geom_vline(
      xintercept = 0, 
      linetype   = "dashed", 
      color      = "black", 
      linewidth  = 0.55
    ) +
    # Thicker dark 95% CI bar touching base of curve with square edges
    annotate(
      "segment",
      x = ci_95[1], xend = ci_95[2],
      y = 0, yend = 0,
      color = "grey20", linewidth = 1.35, lineend = "square"
    ) +
    # Circle dot representing median point estimate
    annotate(
      "point",
      x = pt_est,
      y = 0,
      color = "grey20", size = 1.8, shape = 16
    ) +
    coord_cartesian(ylim = c(-0.02 * max_dens, max_dens * 1.05), clip = "off") +
    theme_void() +
    theme(
      panel.background = element_blank(),
      plot.background  = element_blank(),
      plot.margin      = margin(0, 0, 0, 0)
    )
}

# Main scatterplot with slightly thicker dashed regression guide
create_scatter_panel <- function(data, x_var, y_var, x_label, y_label, point_color) {
  ggplot(data, aes(x = .data[[x_var]], y = .data[[y_var]])) +
    geom_point(color = point_color, alpha = 0.55, size = 1.2) +
    # Slightly thicker dashed regression line
    geom_smooth(
      method    = "lm", 
      se        = FALSE, 
      linetype  = "dashed", 
      color     = "black", 
      linewidth = 0.88
    ) +
    labs(x = x_label, y = y_label) +
    theme_cowplot() +
    theme(
      panel.background = element_rect(fill = "white", color = NA),
      plot.background  = element_rect(fill = "white", color = NA),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line        = element_line(color = "black", linewidth = 0.6),
      axis.title       = element_text(size = 13.5, face = "plain", color = "black"),
      axis.text        = element_text(size = 11.5, face = "plain", color = "black"),
      axis.title.x     = element_text(size = 13.5, face = "plain", color = "black"),
      axis.title.y     = element_text(size = 13.5, face = "plain", color = "black"),
      axis.text.x      = element_text(size = 11.5, face = "plain", color = "black"),
      axis.text.y      = element_text(size = 11.5, face = "plain", color = "black")
    )
}

# ------------------------------------------------------------------------------
# 4. Panel C: Offspring vs. SSB (Yellow / goldenrod2) - Inset Top Right
# ------------------------------------------------------------------------------
scatter_c <- create_scatter_panel(
  data        = intercepts_df,
  x_var       = "SSB_std",
  y_var       = "Offspring_std",
  x_label     = "Annual SSB Propensity (std)",
  y_label     = "Annual Offspring Sired (std)",
  point_color = "goldenrod2"
)

inset_c <- create_density_inset(
  draws_vector = posterior_draws$cor_Animal_ID__SSBRaw_Intercept__numoffspringallbirths_Intercept,
  fill_color   = "goldenrod2"
)

panel_c <- ggdraw(scatter_c) +
  draw_plot(inset_c, x = 0.64, y = 0.70, width = 0.31, height = 0.26) +
  draw_plot_label(label = "C", x = 0.02, y = 0.98, size = panel_label_size, fontface = "bold")

# ------------------------------------------------------------------------------
# 5. Panel D: Offspring vs. DSB (Purple / darkorchid2) - Inset Top Left
# ------------------------------------------------------------------------------
scatter_d <- create_scatter_panel(
  data        = intercepts_df,
  x_var       = "DSB_std",
  y_var       = "Offspring_std",
  x_label     = "Annual DSB Propensity (std)",
  y_label     = "Annual Offspring Sired (std)",
  point_color = "darkorchid2"
)

inset_d <- create_density_inset(
  draws_vector = posterior_draws$cor_Animal_ID__DSBRaw_Intercept__numoffspringallbirths_Intercept,
  fill_color   = "darkorchid2"
)

panel_d <- ggdraw(scatter_d) +
  draw_plot(inset_d, x = 0.18, y = 0.70, width = 0.31, height = 0.26) +
  draw_plot_label(label = "D", x = 0.02, y = 0.98, size = panel_label_size, fontface = "bold")

# ------------------------------------------------------------------------------
# 6. Panel E: SSB vs. DSB (darkgrey) - Inset Top Right
# ------------------------------------------------------------------------------
scatter_e <- create_scatter_panel(
  data        = intercepts_df,
  x_var       = "DSB_std",
  y_var       = "SSB_std",
  x_label     = "Annual DSB Propensity (std)",
  y_label     = "Annual SSB Propensity (std)",
  point_color = "darkgrey"
)

inset_e <- create_density_inset(
  draws_vector = posterior_draws$cor_Animal_ID__SSBRaw_Intercept__DSBRaw_Intercept,
  fill_color   = "darkgrey"
)

panel_e <- ggdraw(scatter_e) +
  draw_plot(inset_e, x = 0.64, y = 0.70, width = 0.31, height = 0.26) +
  draw_plot_label(label = "E", x = 0.02, y = 0.98, size = panel_label_size, fontface = "bold")

# ------------------------------------------------------------------------------
# 7. Save Individual Panels as RDS Files
# ------------------------------------------------------------------------------
saveRDS(panel_c, file = file.path(output_plots_dir, "panel_c_offspring_ssb.rds"))
saveRDS(panel_d, file = file.path(output_plots_dir, "panel_d_offspring_dsb.rds"))
saveRDS(panel_e, file = file.path(output_plots_dir, "panel_e_ssb_dsb.rds"))

cat(sprintf("Panels C, D, and E saved successfully as .rds to:\n%s\n", output_plots_dir))

# ------------------------------------------------------------------------------
# 8. Interactive Session Preview (Combine C, D, E)
# ------------------------------------------------------------------------------
composite_c_d_e <- plot_grid(
  panel_c,
  panel_d,
  panel_e,
  nrow  = 1,
  ncol  = 3,
  align = "h",
  axis  = "tb"
)

print(composite_c_d_e)

ggsave(
  filename = "/SSB_delayed_fitness/Fig1_cde.pdf", 
  plot     = composite_c_d_e, 
  height   = 4, 
  width    = 12
)