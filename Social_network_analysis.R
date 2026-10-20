####### Functional Benefits of SSB: Network Nodal Regressions #######
## Author: Louis Bell-Roberts
## Date: 26/05/2026
## Description: Complete pipeline. Generates networks, runs 3 
## permutations of nodal regressions dynamically, and auto-generates 
## publication-ready supplementary tables.

# ==============================================================================
# 1. SETUP & DATA
# ==============================================================================
library(posterior)
library(bisonR)
library(igraph) 
library(brms)
library(magrittr)
library(tidyr)
library(dplyr)
library(stringi)
library(future)
library(purrr)
library(ggplot2)
library(cowplot)

plan(multisession, workers = 10)
base_dir <- "/Volumes/Macintosh HD/Users/pn26026/Documents/Sandisk/Imperial_HPC/Dyadic_analyses/"
results_dir <- file.path(base_dir, "Social_network_results_v6")

if(!dir.exists(results_dir)) dir.create(results_dir)

master_dyad_df <- readxl::read_excel(path = "/Users/pn26026/Library/CloudStorage/OneDrive-UniversityofBristol/SSB/Data/Supplementary_Tables_Collated_v1.3.xlsx", sheet = 17)

# ==============================================================================
# 2. FILTERING & PRE-PROCESSING
# ==============================================================================
# min_scan_threshold >= 10

indiv_observability <- readxl::read_excel(path = "/Users/pn26026/Library/CloudStorage/OneDrive-UniversityofBristol/SSB/Data/Supplementary_Tables_Collated_v1.3.xlsx", sheet = 20)

valid_monkeys <- readxl::read_excel(path = "/Users/pn26026/Library/CloudStorage/OneDrive-UniversityofBristol/SSB/Data/Supplementary_Tables_Collated_v1.3.xlsx", sheet = 21)

group_membership <- readxl::read_excel(path = "/Users/pn26026/Library/CloudStorage/OneDrive-UniversityofBristol/SSB/Data/Supplementary_Tables_Collated_v1.3.xlsx", sheet = 22)

directed_interactions_long <- bind_rows(
  master_dyad_df %>%
    dplyr::select(Year, Group, dyad, Shared_Scans,
                  node_1_id = ID1, node_2_id = ID2,
                  actor_rank = DS_ID1, actor_age = Age_ID1,
                  Grooming = ID1_to_ID2_Gr, Aggression = ID1_to_ID2_Ag,
                  Mounting = ID1_to_ID2_M, Coalition = Coalition_count),
  
  master_dyad_df %>%
    dplyr::select(Year, Group, dyad, Shared_Scans,
                  node_1_id = ID2, node_2_id = ID1,
                  actor_rank = DS_ID2, actor_age = Age_ID2,
                  Grooming = ID2_to_ID1_Gr, Aggression = ID2_to_ID1_Ag,
                  Mounting = ID2_to_ID1_M, Coalition = Coalition_count)
)

# ==============================================================================
# 3. NETWORK GENERATION LOOP
# ==============================================================================
bison_modfit_list <- list()
bison_strength_list <- list()
bison_eigen_list <- list()
bison_degree_list <- list()
degree_thresholds <- c(0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05)

years <- c(2023, 2024, 2025)
groups <- c("R", "V")
behaviors <- c("Grooming", "Aggression", "Mounting", "Coalition")

for (yr in years) {
  for (grp in groups) {
    for (beh in behaviors) {
      model_key <- paste(yr, grp, beh, sep = "_")
      valid_males <- group_membership %>% filter(Year == yr, Group == grp) %>% pull(ID) %>% sort()
      if(length(valid_males) < 2) next
      
      is_directed <- ifelse(beh == "Coalition", FALSE, TRUE)
      
      full_grid <- expand.grid(node_1_id = valid_males, node_2_id = valid_males, stringsAsFactors = FALSE) %>%
        mutate(Year = yr)
      if(is_directed) { full_grid <- full_grid %>% filter(node_1_id != node_2_id) } else { full_grid <- full_grid %>% filter(node_1_id < node_2_id) }
      
      bison_ready_df <- full_grid %>%
        left_join(directed_interactions_long, by = c("Year", "node_1_id", "node_2_id")) %>%
        mutate(event = replace_na(.data[[beh]], 0), Shared_Scans = replace_na(Shared_Scans, 0)) %>%
        left_join(indiv_observability, by = c("Year", "node_1_id" = "Individual_ID")) %>% rename(Total_Scans_Sender = Total_Indiv_Scans) %>%
        left_join(indiv_observability, by = c("Year", "node_2_id" = "Individual_ID")) %>% rename(Total_Scans_Receiver = Total_Indiv_Scans) %>%
        mutate(Total_Scans_Sender = replace_na(Total_Scans_Sender, 0), Total_Scans_Receiver = replace_na(Total_Scans_Receiver, 0),
               duration = Total_Scans_Sender + Total_Scans_Receiver - Shared_Scans) %>%
        dplyr::select(node_1_id, node_2_id, event, duration)
      
      priors <- get_default_priors("count_conjugate")
      priors$edge <- "gamma(0.1, 0.1)"
      
      bison_model_fit <- bison_model(event | duration ~ dyad(node_1_id, node_2_id), data = bison_ready_df, model_type = "count_conjugate", priors = priors, directed = is_directed)
      
      bison_modfit_list[[model_key]] <- bison_model_fit
      bison_strength_list[[model_key]] <- extract_metric(bison_model_fit, "node_strength")[1:100, ]
      bison_eigen_list[[model_key]] <- extract_metric(bison_model_fit, "node_eigen")[1:100, ]
      
      bison_degree_list[[model_key]] <- list()
      for (th in degree_thresholds) {
        bison_degree_list[[model_key]][[as.character(th)]] <- extract_metric(bison_model_fit, paste0("node_degree[", th, "]"))[1:100, ]
      }
    }
  }
}

# ==============================================================================
# 4. PREPARE UNIFIED NODAL DATAFRAME
# ==============================================================================
nodal_summary <- directed_interactions_long %>% group_by(Year, Group, node_1_id) %>%
  summarise(Mounts_Given = sum(Mounting, na.rm = TRUE), Grooming_Given = sum(Grooming, na.rm = TRUE), Aggression_Given = sum(Aggression, na.rm = TRUE), Total_Coalition = sum(Coalition, na.rm = TRUE), Rank = first(actor_rank), Age = first(actor_age), .groups = "drop") %>% rename(ID = node_1_id)

received_summary <- directed_interactions_long %>% group_by(Year, Group, node_2_id) %>%
  summarise(Mounts_Received = sum(Mounting, na.rm = TRUE), Grooming_Received = sum(Grooming, na.rm = TRUE), Aggression_Received = sum(Aggression, na.rm = TRUE), .groups = "drop") %>% rename(ID = node_2_id)

full_nodal_data <- full_join(nodal_summary, received_summary, by = c("Year", "Group", "ID")) %>%
  mutate(Mounts_Given = replace_na(Mounts_Given, 0), Mounts_Received = replace_na(Mounts_Received, 0), Total_Mounts = Mounts_Given + Mounts_Received,
         Grooming_Given = replace_na(Grooming_Given, 0), Grooming_Received = replace_na(Grooming_Received, 0), Total_Grooming = Grooming_Given + Grooming_Received,
         Aggression_Given = replace_na(Aggression_Given, 0), Aggression_Received = replace_na(Aggression_Received, 0), Total_Aggression = Aggression_Given + Aggression_Received,
         Total_Coalition = replace_na(Total_Coalition, 0)) %>%
  left_join(indiv_observability, by = c("Year", "ID" = "Individual_ID")) %>% rename(Total_Scans_Actor = Total_Indiv_Scans) %>% arrange(ID)

# ==============================================================================
# 5. BUILD IMPUTED DATASETS PER BEHAVIOR
# ==============================================================================
num_draws <- 100
target_threshold <- "0.01"
safe_scale <- function(x) { if (sd(x, na.rm = TRUE) == 0) return(rep(0, length(x))) else return(as.numeric(scale(x))) }

net_vs_net_datasets <- list()

for (i in 1:num_draws) {
  draw_list <- list()
  for (yr in years) {
    for (grp in groups) {
      base_data <- full_nodal_data %>% filter(Year == yr, Group == grp) %>% dplyr::select(ID, Year, Group, Rank, Age)
      gr_key <- paste(yr, grp, "Grooming", sep = "_")
      ag_key <- paste(yr, grp, "Aggression", sep = "_")
      co_key <- paste(yr, grp, "Coalition", sep = "_")
      mo_key <- paste(yr, grp, "Mounting", sep = "_")
      
      if (all(c(gr_key, ag_key, co_key, mo_key) %in% names(bison_modfit_list))) {
        metrics_df <- data.frame(
          ID = names(bison_modfit_list[[mo_key]]$node_to_idx),
          Groom_Strength = as.numeric(bison_strength_list[[gr_key]][i, ]), Groom_Degree = as.numeric(bison_degree_list[[gr_key]][[target_threshold]][i, ]), Groom_Eigen = as.numeric(bison_eigen_list[[gr_key]][i, ]),
          Agg_Strength   = as.numeric(bison_strength_list[[ag_key]][i, ]), Agg_Degree   = as.numeric(bison_degree_list[[ag_key]][[target_threshold]][i, ]), Agg_Eigen   = as.numeric(bison_eigen_list[[ag_key]][i, ]),
          Coal_Strength  = as.numeric(bison_strength_list[[co_key]][i, ]), Coal_Degree  = as.numeric(bison_degree_list[[co_key]][[target_threshold]][i, ]), Coal_Eigen  = as.numeric(bison_eigen_list[[co_key]][i, ]),
          Mount_Strength = as.numeric(bison_strength_list[[mo_key]][i, ]), Mount_Degree = as.numeric(bison_degree_list[[mo_key]][[target_threshold]][i, ]), Mount_Eigen = as.numeric(bison_eigen_list[[mo_key]][i, ])
        )
        
        draw_list[[paste(yr, grp)]] <- left_join(base_data, metrics_df, by = "ID") %>%
          mutate(
            Rank_z = safe_scale(Rank),
            Group_Year_Fct = as.factor(paste(Group, Year, sep = "_")),
            
            Groom_Strength_z = safe_scale(Groom_Strength), Groom_Degree_z = safe_scale(Groom_Degree), Groom_Eigen_z = safe_scale(Groom_Eigen),
            Agg_Strength_z   = safe_scale(Agg_Strength), Agg_Degree_z   = safe_scale(Agg_Degree), Agg_Eigen_z   = safe_scale(Agg_Eigen),
            Coal_Strength_z  = safe_scale(Coal_Strength), Coal_Degree_z  = safe_scale(Coal_Degree), Coal_Eigen_z  = safe_scale(Coal_Eigen),
            Mount_Strength_z = safe_scale(Mount_Strength), Mount_Eigen_z  = safe_scale(Mount_Eigen)
          )
      }
    }
  }
  net_vs_net_datasets[[i]] <- bind_rows(draw_list)
}

filter_imputations <- function(dataset_list, exclusions) { lapply(dataset_list, function(df) { df %>% filter(!(paste(Year, Group, sep = "_") %in% exclusions)) }) }

# ==============================================================================
# 6. DYNAMIC MODEL PERMUTATION (GLOBAL MIXED ONLY)
# ==============================================================================

# Define Mixed-Effect Priors (Fixed-effect priors removed as they are no longer needed)
priors_list <- list(
  NegBin_Mixed = c(
    set_prior("normal(0,1)",      class = "Intercept"),
    set_prior("normal(0,0.5)",    class = "b"),
    set_prior("exponential(2)",   class = "sd"),
    set_prior("exponential(0.1)", class = "shape")
  ),
  Gauss_Mixed = c(
    set_prior("normal(0,1)",      class = "Intercept"),
    set_prior("normal(0,0.5)",    class = "b"),
    set_prior("exponential(2)",   class = "sd"),
    set_prior("exponential(1)",   class = "sigma")
  )
)

# Base Predictor Formula
global_rhs <- paste(
  "1 + Rank_z + scale(Age) +",
  "Groom_Degree_z + Agg_Degree_z + Coal_Degree_z +",
  "Groom_Strength_z + Agg_Strength_z + Coal_Strength_z +",
  "Groom_Eigen_z + Agg_Eigen_z + Coal_Eigen_z +",
  "Group_Year_Fct + (1 | ID)"
)

# Define the 3 Responses with their families and mixed priors
response_vars <- list(
  Mount_Degree     = list(family = negbinomial(), prior = priors_list$NegBin_Mixed),
  Mount_Strength_z = list(family = gaussian(),    prior = priors_list$Gauss_Mixed),
  Mount_Eigen_z    = list(family = gaussian(),    prior = priors_list$Gauss_Mixed)
)

# Build list containing the 3 Global Mixed models
model_manifest <- list()

for (resp in names(response_vars)) {
  mod_key <- paste0(resp, "_Global_Mixed")
  
  model_manifest[[mod_key]] <- list(
    formula  = as.formula(paste(resp, "~", global_rhs)),
    family   = response_vars[[resp]]$family,
    prior    = response_vars[[resp]]$prior,
    filename = paste0("mdl_", mod_key, ".rds")
  )
}


# ==============================================================================
# 7. EXECUTE BRMS MODELS
# ==============================================================================
con <- file(file.path(base_dir, "Social_network_results_v6/model_output_summaries.txt"), open = "wt")
sink(con)
sink(con, type = "message")

print("Starting execution of 3 permutations...")

for (mod_name in names(model_manifest)) {
  cfg <- model_manifest[[mod_name]]
  filepath <- file.path(results_dir, cfg$filename)
  
  if(!file.exists(filepath)) {
    cat(paste("\n--- FITTING MODEL:", mod_name, "---\n"))
    mdl <- brm_multiple(
      formula = cfg$formula,
      data = filter_imputations(net_vs_net_datasets, c("2024_V", "2025_R")),
      family = cfg$family,
      prior = cfg$prior,
      init = 0,                                 
      chains = 1, cores = 1, iter = 4000, warmup = 2000,
      control = list(adapt_delta = 0.999, max_treedepth = 15)
    )
    saveRDS(mdl, filepath)
    rm(mdl)
    gc()
  } else {
    cat(paste("\n--- SKIPPING:", mod_name, "(Already exists on disk) ---\n"))
  }
}

# Shut down the multisession workers before compiling brms models to free up RAM
future::plan(sequential)
gc()

# ==============================================================================
# 9. EXTRACT FINAL PARAMETERS & DIAGNOSTIC TABLES
# ==============================================================================
print("Extracting parameters and sampler diagnostics...")

extract_model_all_data <- function(mod_key, cfg) {
  filepath <- file.path(results_dir, cfg$filename)
  if (!file.exists(filepath)) return(NULL)
  
  mdl <- readRDS(filepath)
  model_formula <- paste(trimws(deparse(formula(mdl)$formula)), collapse = " ")
  mdl_summ <- summary(mdl)
  
  # A. Pool Parameters
  p_df <- bind_rows(
    if(!is.null(mdl_summ$fixed)) { as.data.frame(mdl_summ$fixed) %>% tibble::rownames_to_column("Parameter") %>% mutate(Effect_Type = "Fixed") } else { data.frame() },
    if(length(mdl_summ$random) > 0) { map_dfr(names(mdl_summ$random), ~ as.data.frame(mdl_summ$random[[.x]]) %>% tibble::rownames_to_column("Parameter") %>% mutate(Parameter = paste0(Parameter, " [", .x, "]"), Effect_Type = "Random")) } else { data.frame() },
    if(!is.null(mdl_summ$spec_pars) && nrow(mdl_summ$spec_pars) > 0) { as.data.frame(mdl_summ$spec_pars) %>% tibble::rownames_to_column("Parameter") %>% mutate(Effect_Type = "Distributional") } else { data.frame() }
  ) %>%
    mutate(Model = mod_key, Model_Structure = model_formula, Estimate = round(Estimate, 3), Est.Error = round(Est.Error, 3), Lower_95_CI = round(`l-95% CI`, 3), Upper_95_CI = round(`u-95% CI`, 3)) %>%
    select(Model, Model_Structure, Effect_Type, Parameter, Estimate, Est.Error, Lower_95_CI, Upper_95_CI)
  
  # B. Convergence Diagnostics (Rhat/ESS)
  current_draws <- as_draws_array(mdl)
  n_imps <- dim(current_draws)[2]
  diags_df <- map_dfr(1:n_imps, ~ summarise_draws(subset_draws(current_draws, chain = .x), default_convergence_measures()) %>% mutate(Imputation_Dataset = .x)) %>%
    mutate(Model = mod_key) %>% select(Model, Imputation_Dataset, variable, rhat, ess_bulk, ess_tail)
  
  # C. Sampler Divergences
  div_df <- brms::nuts_params(mdl) %>% filter(Parameter == "divergent__") %>% group_by(Chain) %>% summarise(Divergences = sum(Value), .groups = "drop")
  diag_sum <- data.frame(Model = mod_key, Total_Divergent_Transitions = sum(div_df$Divergences), Mean_Divergent_Transitions_per_Imputation = round(mean(div_df$Divergences), 2), Max_Divergent_Transitions_in_Single_Imputation = max(div_df$Divergences))
  
  return(list(params = p_df, archive = diags_df, summary = diag_sum))
}

all_data <- compact(imap(model_manifest, ~extract_model_all_data(.y, .x)))

final_params <- map_dfr(all_data, ~ .x$params)
final_archive <- map_dfr(all_data, ~ .x$archive)

mi_performance <- final_archive %>% group_by(Model) %>%
  summarise(Total_Imputations_Run = n_distinct(Imputation_Dataset), Mean_Rhat = round(mean(rhat, na.rm = TRUE), 4), Max_Rhat = round(max(rhat, na.rm = TRUE), 4), Mean_Bulk_ESS = round(mean(ess_bulk, na.rm = TRUE), 0), Min_Bulk_ESS = round(min(ess_bulk, na.rm = TRUE), 0), Mean_Tail_ESS = round(mean(ess_tail, na.rm = TRUE), 0), Min_Tail_ESS = round(min(ess_tail, na.rm = TRUE), 0), .groups = "drop") %>%
  left_join(map_dfr(all_data, ~ .x$summary), by = "Model")

write.csv(final_params, file.path(results_dir, "Supp_Table_1_Parameters_Final.csv"), row.names = FALSE)
write.csv(mi_performance, file.path(results_dir, "Supp_Table_2_MI_Diagnostics_Final.csv"), row.names = FALSE)

print("PIPELINE COMPLETE.")


# ==============================================================================
# 8. VARIANCE INFLATION FACTOR (VIF) EXPORT
# ==============================================================================
library(car)
vif_data <- filter_imputations(net_vs_net_datasets, c("2024_V", "2025_R"))[[1]]
dummy_lm <- lm(Mount_Degree ~ Rank_z + scale(Age) + Group_Year_Fct +
                 Groom_Degree_z + Agg_Degree_z + Coal_Degree_z + 
                 Groom_Strength_z + Agg_Strength_z + Coal_Strength_z + 
                 Groom_Eigen_z + Agg_Eigen_z + Coal_Eigen_z, data = vif_data)

vif_results <- car::vif(dummy_lm)

# Standardize VIF extraction
vif_df <- data.frame(Predictor = rownames(as.data.frame(vif_results)), VIF_Score = NA)
if(is.matrix(vif_results)) { vif_df$VIF_Score <- round(vif_results[, "GVIF"], 2) } else { vif_df$VIF_Score <- round(vif_results, 2) }
vif_df$Threshold_Flag <- ifelse(vif_df$VIF_Score > 5, "High Collinearity (>5)", "Acceptable")

write.csv(vif_df, file.path(results_dir, "Supp_Table_3_VIF_Diagnostics.csv"), row.names = FALSE)

# Revert outputs back to the console
sink(type = "message")
sink()

# ==============================================================================
#  SPAGHETTI PLOT (SSB STRENGTH VS. GROOMING STRENGTH) - UNLABELLED PANEL
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. Configuration & Directories
# ------------------------------------------------------------------------------
# Output directory for saved figure RDS files
output_plots_dir <- "/Users/pn26026/OneDrive - University of Bristol/SSB/Code/"
if (!dir.exists(output_plots_dir)) dir.create(output_plots_dir, recursive = TRUE)

# Color for SSB trajectory draws (consistent with earlier SSB panels)
line_color <- "goldenrod2"

# ------------------------------------------------------------------------------
# 2. Load Model & Construct Prediction Grid
# ------------------------------------------------------------------------------
# Load the Mount_Strength_z model saved in Section 7
model_path <- file.path(results_dir, "mdl_Mount_Strength_z_Global_Mixed.rds")
mdl_mount_str <- readRDS(model_path)

# Extract empirical range and reference levels from vif_data (Section 8)
groom_range <- seq(
  min(vif_data$Groom_Strength_z, na.rm = TRUE),
  max(vif_data$Groom_Strength_z, na.rm = TRUE),
  length.out = 100
)
mean_age <- mean(vif_data$Age, na.rm = TRUE)

# Extract the exact factor levels the model was trained on
allowed_levels <- levels(mdl_mount_str$data$Group_Year_Fct)
ref_group_year <- allowed_levels[1]

# Prediction grid holding all other covariates at population mean / baseline
grid_data <- tibble(
  Groom_Strength_z = groom_range,
  Rank_z           = 0,
  Age              = mean_age,
  Groom_Degree_z   = 0,
  Agg_Degree_z     = 0,
  Coal_Degree_z    = 0,
  Agg_Strength_z   = 0,
  Coal_Strength_z  = 0,
  Groom_Eigen_z    = 0,
  Agg_Eigen_z      = 0,
  Coal_Eigen_z     = 0,
  Group_Year_Fct   = factor(ref_group_year, levels = allowed_levels),
  ID               = NA
)

# ------------------------------------------------------------------------------
# 3. Generate 300 Posterior Draws & Population Mean
# ------------------------------------------------------------------------------
post_draws <- posterior_epred(
  mdl_mount_str,
  newdata    = grid_data,
  ndraws     = 300,
  re_formula = NA
)

# Population mean trajectory
grid_data$estimate__ <- colMeans(post_draws)

# Reshape draws to long format for spaghetti plotting
draws_df <- as.data.frame(t(post_draws))
colnames(draws_df) <- paste0("draw_", seq_len(ncol(draws_df)))

draws_long <- cbind(grid_data, draws_df) %>%
  pivot_longer(
    cols      = starts_with("draw_"),
    names_to  = "draw",
    values_to = "prediction"
  )

# ------------------------------------------------------------------------------
# 4. Build Spaghetti Plot
# ------------------------------------------------------------------------------
p_spaghetti <- ggplot() +
  # 300 colored posterior draws
  geom_line(
    data      = draws_long,
    aes(x = Groom_Strength_z, y = prediction, group = draw),
    color     = line_color,
    alpha     = 0.10,
    linewidth = 0.30
  ) +
  # Solid black population mean line
  geom_line(
    data      = grid_data,
    aes(x = Groom_Strength_z, y = estimate__),
    color     = "black",
    linewidth = 0.88,
    linetype  = "solid"
  ) +
  labs(
    x = "Grooming Strength (std)",
    y = "SSB Strength (std)"
  ) +
  theme_cowplot() +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.background  = element_rect(fill = "white", color = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line        = element_line(color = "black", linewidth = 0.6),
    # Plain text typography (no bold)
    axis.title       = element_text(size = 13.5, face = "plain", color = "black"),
    axis.text        = element_text(size = 11.5, face = "plain", color = "black"),
    axis.title.x     = element_text(size = 13.5, face = "plain", color = "black"),
    axis.title.y     = element_text(size = 13.5, face = "plain", color = "black"),
    axis.text.x      = element_text(size = 11.5, face = "plain", color = "black"),
    axis.text.y      = element_text(size = 11.5, face = "plain", color = "black")
  )

# ------------------------------------------------------------------------------
# 5. Save Unlabelled Plot as RDS
# ------------------------------------------------------------------------------
# Save pure ggplot object directly so labels/sizes can be added dynamically later
output_rds_path <- file.path(output_plots_dir, "panel_a_ssb_grooming_strength.rds")
saveRDS(p_spaghetti, file = d)

cat(sprintf("\nUnlabelled Panel A spaghetti plot saved successfully to:\n%s\n", output_rds_path))

# Interactive session preview
print(p_spaghetti)
ggsave(filename = "/Users/pn26026/OneDrive - University of Bristol/SSB/Code//panel_a_ssb_grooming_strength.pdf", plot = p_spaghetti, height = 4, width = 4)
