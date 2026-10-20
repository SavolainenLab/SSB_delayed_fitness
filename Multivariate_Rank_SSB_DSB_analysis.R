#### Testing the annual association between SSB, DSB and rank ####
##Louis Bell-Roberts
#23/01/2026


library(tidyverse)
library(brms)
library(ggstance)
library(ggplot2)
library(coda)
library(bayesplot)
library(tidybayes)
library(patchwork)
library(lubridate)
library(posterior)
library(cowplot)
library(magick)
library(pdftools)

#Read in data
# File path (Local working path; convert to relative path prior to final archive)
supp_data_path <- "/Users/pn26026/Library/CloudStorage/OneDrive-UniversityofBristol/SSB/Data/Supplementary_Tables_Collated_v1.3.xlsx"

# Sheet 18: Combined cleaned collapsed data with observation effort
# Sheet 19: Cleaned pedigree data
data <- readxl::read_excel(path = supp_data_path, sheet = 18)
pedigree_pre <- readxl::read_excel(path = supp_data_path, sheet = 19)

# --- 2. Create Analysis DF and JOIN the lookup ---
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
  mutate( #Scale age variables and rank (DS)
    Animal_ID = as.factor(Animal_ID),
    Age_s = scale(Age),
    age_B_s = scale(age_B),
    age_W_s = scale(age_W),
    Group_Year = as.factor(paste(assigned_group, Year, sep = "_")),
    adjusted_SSB = (SSB_Raw)/log(Scan_Count),
    adjusted_DSB = (DSB_Raw)/log(Scan_Count)
  ) %>%
  mutate(
    assigned_group = as.factor(assigned_group),
    Year = as.factor(Year)
  ) %>% 
  group_by(Year, Group_for_DS) %>%
  mutate(
    DS_s = scale(DS)
  ) %>%
  ungroup()

########################
#Multivariate model####
########################
priors_SSB_DSB_rank <- c(
  # --- FIXED EFFECTS ---
  prior(normal(0, 5),   class = Intercept, resp = "SSBRaw"),
  prior(normal(0, 5),   class = Intercept, resp = "DSBRaw"),
  prior(normal(0, 5),   class = Intercept, resp = "DSs"),
  prior(normal(0, 0.5), class = b,         resp = "SSBRaw"),
  prior(normal(0, 0.5), class = b,         resp = "DSBRaw"),
  prior(normal(0, 0.5), class = b,         resp = "DSs"),
  
  # --- RANDOM EFFECTS ---
  # 1. SD for the INTERCEPTS
  prior(exponential(1), class = sd, coef = "Intercept", group = "Animal_ID", resp = "SSBRaw"),
  prior(exponential(1), class = sd, coef = "Intercept", group = "Animal_ID", resp = "DSBRaw"),
  prior(exponential(1), class = sd, coef = "Intercept", group = "Animal_ID", resp = "DSs"),
  prior(lkj(2), class = cor, group = "Animal_ID"),
  
  # assigned_group:Year (Intercepts only)
  prior(exponential(1), class = sd, group = "assigned_group:Year", resp = "SSBRaw"),
  prior(exponential(1), class = sd, group = "assigned_group:Year", resp = "DSBRaw"),
  prior(exponential(1), class = sd, group = "assigned_group:Year", resp = "DSs"),

  # --- DISTRIBUTION PARAMETERS ---
  prior(exponential(1), class = shape, resp = "SSBRaw"),
  prior(exponential(1), class = shape, resp = "DSBRaw"),
  prior(exponential(1), class = sigma, resp = "DSs"))


bf_SSB <- bf(SSB_Raw ~ 1 + age_W_s + age_B_s + offset(log(Scan_Count)) + (1|p|Animal_ID) + (1|assigned_group:Year)) + negbinomial()

bf_DSB <- bf(DSB_Raw ~ 1 + age_W_s + age_B_s + offset(log(Scan_Count)) + (1|p|Animal_ID) + (1|assigned_group:Year)) + negbinomial()

bf_rank <- bf(DS_s ~ 1 + age_W_s + age_B_s + (1|p|Animal_ID) + (1|assigned_group:Year)) + gaussian()

biv_SSB_DSB_rank <- brm(
  bf_SSB + bf_DSB + bf_rank + set_rescor(F),
  prior = priors_SSB_DSB_rank,
  data = analysis_df, seed = 1,
  warmup = 1000, iter = 4000, thin = 1, chains = 4, cores = 4,
  control = list(adapt_delta = 0.9)
)

summary(biv_SSB_DSB_rank)

biv_SSB_DSB_rank <- readRDS(file = "/Users/pn26026/OneDrive - University of Bristol/SSB/Rscripts_local/Male_fitness/All_births/biv_SSB_DSB_rank_mod.rds")

pp_check(biv_SSB_DSB_rank, resp='DSs') + xlim(c(-3, 2.5))
pp_check(biv_SSB_DSB_rank, resp='SSBRaw', type = "bars") + xlim(c(-1, 60))
pp_check(biv_SSB_DSB_rank, resp='DSBRaw', type = "bars") + xlim(c(-1, 50))

# --- 1. Extract brms Summary ---
# Extract summary from the new rank model. 
# The summary contains ESS and Rhat values automatically.
mod_sum <- summary(biv_SSB_DSB_rank)

# --- 2. Extract and Format Fixed Effects ---
# Convert the fixed effects matrix to a data frame and add metadata columns
fixed_df <- as.data.frame(mod_sum$fixed) %>%
  tibble::rownames_to_column("Raw_Parameter") %>%
  dplyr::mutate(
    Effect_Type = "Fixed",
    Group = NA_character_ # Fixed effects do not have a grouping level
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
    # Identify which trait the parameter belongs to, flagging correlations specifically.
    Response_Variable = dplyr::case_when(
      str_detect(Raw_Parameter, "^cor\\(") ~ "Multivariate (Correlations)",
      str_detect(Raw_Parameter, "SSBRaw") ~ "SSB",
      str_detect(Raw_Parameter, "DSBRaw") ~ "DSB",
      str_detect(Raw_Parameter, "DSs") ~ "Rank", 
      TRUE ~ "Other"
    ),
    
    # Translate syntax into readable labels
    Parameter = Raw_Parameter,
    
    # 1. Clean random effect correlation wrappers first
    Parameter = str_replace_all(Parameter, "^cor\\((.*)\\)$", "Cor: \\1"),
    # Replace the comma separator in correlation terms with a bidirectional arrow
    Parameter = ifelse(str_detect(Parameter, "^Cor: "), str_replace_all(Parameter, ",", " <-> "), Parameter),
    
    # 2. Clean random effect Standard Deviation wrappers
    Parameter = str_replace_all(Parameter, "^sd\\((.*)\\)$", "SD: \\1"), 
    
    # 3. Strip trait prefixes
    Parameter = ifelse(
      Response_Variable == "Multivariate (Correlations)",
      str_replace_all(Parameter, "SSBRaw_", "SSB ") %>%
        str_replace_all("DSBRaw_", "DSB ") %>%
        str_replace_all("DSs_", "Rank "),
      str_replace_all(Parameter, "(SSBRaw_|DSBRaw_|DSs_)", "")
    ),
    
    # 4. Clean predictor labels
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
    # Round ESS to whole numbers
    dplyr::across(c(ESS_Bulk, ESS_Tail), ~ round(.x, 0)) 
  )

# --- 6. Export the Final Unified Table ---
# Save the combined dataframe as a single CSV file
# write.csv(combined_table, "/Users/louis.bell-roberts/Documents/SSB_Macaques/Rscripts_local/Male_fitness/All_births/Output_table_SSB_DSB_rank.csv", row.names = FALSE)

# Plotting ####
# ==============================================================================
# Plotting: 3-Panel Longitudinal Age Trajectories (SSB, DSB, Rank)
# All panels derived from single model: biv_SSB_DSB_rank
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. Prediction Grid: Age 3.5 to 25 years
# ------------------------------------------------------------------------------
mean_age <- mean(analysis_df$Age, na.rm = TRUE)
sd_age_w <- sd(analysis_df$age_W, na.rm = TRUE)

age_grid_base <- data.frame(
  Age = seq(3.5, 25, length.out = 100)
) %>%
  mutate(
    # Standardised within-individual deviation
    age_W_s    = (Age - mean_age) / sd_age_w,
    # Population average cohort effect
    age_B_s    = 0,
    # 100 scans for count exposure (ignored by Gaussian rank sub-model)
    Scan_Count = 100
  )

# ------------------------------------------------------------------------------
# 2. Reusable Spaghetti Plotting Function (Plain text axes, 300 draws)
# ------------------------------------------------------------------------------
create_age_spaghetti_panel <- function(model, 
                                       grid_data, 
                                       y_resp, 
                                       y_label, 
                                       line_color, 
                                       y_limits = NULL,
                                       mean_linewidth = 0.8) {
  
  # A. Generate 300 posterior draws (population-level prediction)
  post_draws <- posterior_epred(
    model,
    newdata    = grid_data,
    resp       = y_resp,
    ndraws     = 300,
    re_formula = NA
  )
  
  # B. Calculate population mean trajectory
  grid_data$estimate__ <- colMeans(post_draws)
  
  # C. Reshape draws to long format
  draws_df <- as.data.frame(t(post_draws))
  colnames(draws_df) <- paste0("draw_", seq_len(ncol(draws_df)))
  
  draws_long <- cbind(grid_data, draws_df) %>%
    pivot_longer(
      cols      = starts_with("draw_"),
      names_to  = "draw",
      values_to = "prediction"
    )
  
  # D. Build ggplot with plain (non-bold) text styling
  p <- ggplot() +
    # 300 posterior draws
    geom_line(
      data      = draws_long,
      aes(x = Age, y = prediction, group = draw),
      color     = line_color,
      alpha     = 0.10,
      linewidth = 0.30
    ) +
    # Population mean trajectory
    geom_line(
      data      = grid_data,
      aes(x = Age, y = estimate__),
      color     = "black",
      linewidth = mean_linewidth,
      linetype  = "solid"
    ) +
    labs(x = "Age (years)", y = y_label) +
    theme_cowplot() +
    theme(
      panel.background = element_rect(fill = "white", color = NA),
      plot.background  = element_rect(fill = "white", color = NA),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line        = element_line(color = "black", linewidth = 0.6),
      # Strict plain text for titles and tick numbers (no bold)
      axis.title       = element_text(size = 12, face = "plain", color = "black"),
      axis.text        = element_text(size = 11, face = "plain", color = "black"),
      axis.title.x     = element_text(size = 12, face = "plain", color = "black"),
      axis.title.y     = element_text(size = 12, face = "plain", color = "black"),
      axis.text.x      = element_text(size = 11, face = "plain", color = "black"),
      axis.text.y      = element_text(size = 11, face = "plain", color = "black")
    ) +
    # Ticks start at 5, but the axis coordinate framing starts at 3.5
    scale_x_continuous(
      breaks = seq(5, 25, by = 5),
      expand = c(0, 0)
    )
  
  if (!is.null(y_limits)) {
    p <- p + coord_cartesian(xlim = c(3.5, 25), ylim = y_limits)
  } else {
    p <- p + coord_cartesian(xlim = c(3.5, 25))
  }
  
  return(p)
}

# ------------------------------------------------------------------------------
# 3. Generate All Three Panels from biv_SSB_DSB_rank
# ------------------------------------------------------------------------------
# Panel A: SSB (Yellow / goldenrod2)
panel_a_ssb <- create_age_spaghetti_panel(
  model      = biv_SSB_DSB_rank,
  grid_data  = age_grid_base,
  y_resp     = "SSBRaw",
  y_label    = "Annual SSB mounts",
  line_color = "goldenrod2",
  y_limits   = c(0, 25)
)

# Panel B: DSB (Purple / darkorchid2)
panel_b_dsb <- create_age_spaghetti_panel(
  model      = biv_SSB_DSB_rank,
  grid_data  = age_grid_base,
  y_resp     = "DSBRaw",
  y_label    = "Annual DSB mounts",
  line_color = "darkorchid2",
  y_limits   = c(0, 25)
)

# Panel C: Dominance Rank (Blue)
panel_c_rank <- create_age_spaghetti_panel(
  model      = biv_SSB_DSB_rank,
  grid_data  = age_grid_base,
  y_resp     = "DSs",
  y_label    = "Dominance rank (std)",
  line_color = "blue",
  y_limits   = c(-1.5, 2.5)
)

# ------------------------------------------------------------------------------
# 4. Assemble Multi-Panel Figure (A, B, C)
# ------------------------------------------------------------------------------
three_panel_age_plot <- plot_grid(
  panel_a_ssb,
  panel_b_dsb,
  panel_c_rank,
  nrow           = 1,
  ncol           = 3,
  labels         = c("A", "B", "C"),
  label_size     = 14,
  label_fontface = "bold",   # Retains bold for panel identifier tags
  align          = "h",
  axis           = "tb"
)

print(three_panel_age_plot)

# ------------------------------------------------------------------------------
# 5. Export PDF
# ------------------------------------------------------------------------------
# ggsave(
#   filename = "/Users/pn26026/OneDrive - University of Bristol/SSB/Code/Fig2_panels_a_b_c_vector.pdf",
#   plot     = three_panel_age_plot,
#   width    = 12,
#   height   = 4,
#   dpi      = 300
# )

# ==============================================================================
# SECTION: RANDOM INTERCEPTS PLOTS (PANELS F & G) & COMPOSITE C–G
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. Output Directory Configuration & Aesthetic Settings
# ------------------------------------------------------------------------------
output_plots_dir <- "/Users/pn26026/OneDrive - University of Bristol/SSB/Code/"
if (!dir.exists(output_plots_dir)) dir.create(output_plots_dir, recursive = TRUE)

# Panel tag size
panel_label_size <- 28

# ------------------------------------------------------------------------------
# 2. Extract Raw (Unscaled) Animal_ID Random Intercepts
# ------------------------------------------------------------------------------
ran_vals_rank <- ranef(biv_SSB_DSB_rank)$Animal_ID

intercepts_rank_df <- tibble(
  Animal_ID = rownames(ran_vals_rank),
  SSB_std   = as.numeric(ran_vals_rank[, "Estimate", "SSBRaw_Intercept"]),
  DSB_std   = as.numeric(ran_vals_rank[, "Estimate", "DSBRaw_Intercept"]),
  Rank_std  = as.numeric(ran_vals_rank[, "Estimate", "DSs_Intercept"])
)

posterior_draws_rank <- as_draws_df(biv_SSB_DSB_rank)

cor_ssb_rank_col <- names(posterior_draws_rank)[
  str_detect(names(posterior_draws_rank), "cor_Animal_ID") & 
    str_detect(names(posterior_draws_rank), "SSBRaw_Intercept") & 
    str_detect(names(posterior_draws_rank), "DSs_Intercept")
][1]

cor_dsb_rank_col <- names(posterior_draws_rank)[
  str_detect(names(posterior_draws_rank), "cor_Animal_ID") & 
    str_detect(names(posterior_draws_rank), "DSBRaw_Intercept") & 
    str_detect(names(posterior_draws_rank), "DSs_Intercept")
][1]

# ------------------------------------------------------------------------------
# 3. Helper Functions: Inset Density & Scatterplot
# ------------------------------------------------------------------------------
create_density_inset <- function(draws_vector, fill_color) {
  df <- data.frame(val = draws_vector)
  
  ci_95  <- quantile(draws_vector, probs = c(0.025, 0.975))
  pt_est <- median(draws_vector)
  
  dens     <- density(draws_vector)
  max_dens <- max(dens$y)
  
  ggplot(df, aes(x = val)) +
    geom_density(
      fill  = fill_color, 
      alpha = 1, 
      color = NA
    ) +
    geom_vline(
      xintercept = 0, 
      linetype   = "dashed", 
      color      = "black", 
      linewidth  = 0.55
    ) +
    annotate(
      "segment",
      x = ci_95[1], xend = ci_95[2],
      y = 0, yend = 0,
      color = "grey20", linewidth = 1.35, lineend = "square"
    ) +
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

create_scatter_panel <- function(data, x_var, y_var, x_label, y_label, point_color) {
  ggplot(data, aes(x = .data[[x_var]], y = .data[[y_var]])) +
    geom_point(color = point_color, alpha = 0.55, size = 1.2) +
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
# 4. Panel F: Rank vs. SSB (Yellow / goldenrod2) - Inset Top Right
# ------------------------------------------------------------------------------
scatter_f <- create_scatter_panel(
  data        = intercepts_rank_df,
  x_var       = "SSB_std",
  y_var       = "Rank_std",
  x_label     = "Annual SSB Propensity (std)",
  y_label     = "Annual Rank (std)",
  point_color = "goldenrod2"
)

inset_f <- create_density_inset(
  draws_vector = posterior_draws_rank[[cor_ssb_rank_col]],
  fill_color   = "goldenrod2"
)

panel_f <- ggdraw(scatter_f) +
  draw_plot(inset_f, x = 0.64, y = 0.70, width = 0.31, height = 0.26) +
  draw_plot_label(label = "F", x = 0.02, y = 0.98, size = panel_label_size, fontface = "bold")

# ------------------------------------------------------------------------------
# 5. Panel G: Rank vs. DSB (Purple / darkorchid2) - Inset Top Left
# ------------------------------------------------------------------------------
scatter_g <- create_scatter_panel(
  data        = intercepts_rank_df,
  x_var       = "DSB_std",
  y_var       = "Rank_std",
  x_label     = "Annual DSB Propensity (std)",
  y_label     = "Annual Rank (std)",
  point_color = "darkorchid2"
)

inset_g <- create_density_inset(
  draws_vector = posterior_draws_rank[[cor_dsb_rank_col]],
  fill_color   = "darkorchid2"
)

panel_g <- ggdraw(scatter_g) +
  draw_plot(inset_g, x = 0.18, y = 0.70, width = 0.31, height = 0.26) +
  draw_plot_label(label = "G", x = 0.02, y = 0.98, size = panel_label_size, fontface = "bold")

# Save Panels F & G
saveRDS(panel_f, file = file.path(output_plots_dir, "panel_f_rank_ssb.rds"))
saveRDS(panel_g, file = file.path(output_plots_dir, "panel_g_rank_dsb.rds"))

# ------------------------------------------------------------------------------
# 6. Load Panels C, D, E & Assemble Multi-Panel Figure (C–G)
# ------------------------------------------------------------------------------
panel_c <- readRDS(file.path(output_plots_dir, "panel_c_offspring_ssb.rds"))
panel_d <- readRDS(file.path(output_plots_dir, "panel_d_offspring_dsb.rds"))
panel_e <- readRDS(file.path(output_plots_dir, "panel_e_ssb_dsb.rds"))

top_row <- plot_grid(
  panel_c,
  panel_d,
  panel_e,
  nrow  = 1,
  ncol  = 3,
  align = "h",
  axis  = "tb"
)

bottom_row <- plot_grid(
  panel_f,
  panel_g,
  NULL,
  nrow  = 1,
  ncol  = 3,
  align = "h",
  axis  = "tb"
)

composite_c_to_g <- plot_grid(
  top_row,
  bottom_row,
  nrow        = 2,
  ncol        = 1,
  rel_heights = c(1, 1)
)

print(composite_c_to_g)

ggsave(
  filename = file.path(output_plots_dir, "Fig1_panels_c_d_e_f_g_vector.pdf"),
  plot     = composite_c_to_g,
  width    = 12,
  height   = 8,
  dpi      = 300
)