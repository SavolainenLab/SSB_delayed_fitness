####### Investigate Mount Time Difference Distributions #######
## Author: Louis Bell-Roberts
## Description: Calculates inter-mount intervals and bout distributions for
## directional dyads and generates comparison histograms for DSB and SSB.
# ==============================================================================

# Load necessary libraries
library(tidyverse)
library(data.table)
library(lubridate)
library(readxl)
library(patchwork)

# ------------------------------------------------------------------------------
# 0. Repository Paths
# ------------------------------------------------------------------------------
repo_root    <- "/Users/pn26026/OneDrive - University of Bristol/SSB/Code/Github_upload_v1.1/"
supp_xlsx    <- file.path(repo_root, "Data_gen/Data/Raw_inputs/Supplementary_Tables_Collated_v1.4.xlsx")
plot_out_pdf <- file.path(repo_root, "Mount_time_diff_panel_plot.pdf")

# ------------------------------------------------------------------------------
# 1. Ingest Pedigree & Define Anonymisation Placeholders
# ------------------------------------------------------------------------------
pedigree <- as.data.frame(readxl::read_excel(supp_xlsx, sheet = 20, col_types = "text")) %>%
  mutate(
    AnimalID = trimws(as.character(AnimalID)),
    Sex      = toupper(trimws(as.character(Sex)))
  )

all_males <- pedigree %>%
  filter(Sex == "M") %>%
  pull(AnimalID) %>%
  unique()

all_females <- pedigree %>%
  filter(Sex == "F") %>%
  pull(AnimalID) %>%
  unique()

# Placeholder restoration vector to match field ad-lib notation
placeholder_restore <- c(
  "MNQE" = "FUN",
  "LTFO" = "MUN",
  "VPBE" = "SAMUN",
  "GLXD" = "JVUN",
  "HYCZ" = "UN"
)

# ------------------------------------------------------------------------------
# 2. Ingest and Clean Ad-Lib Mounting Data (Sheet 13)
# ------------------------------------------------------------------------------
master_data_pre <- as.data.frame(readxl::read_excel(supp_xlsx, sheet = 13, col_types = "text")) %>%
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

# Filter for valid mounts and clean timestamps
mounts_raw <- master_data_pre %>%
  filter(
    Action == "M",
    !is.na(Date), !is.na(Time),
    !is.na(ACTOR1), !is.na(RECIPIENT),
    ACTOR1 != "", RECIPIENT != "",
    ACTOR1 != RECIPIENT
  ) %>%
  mutate(
    Date_clean = substr(as.character(Date), 1, 10),
    Time_clean = trimws(as.character(Time)),
    Time_clean = ifelse(grepl(" ", Time_clean), sub(".* ", "", Time_clean), Time_clean),
    DateTime   = parse_date_time(paste(Date_clean, Time_clean), orders = c("ymd HMS", "ymd HM"), tz = "UTC")
  ) %>%
  filter(!is.na(DateTime))

# ------------------------------------------------------------------------------
# 3. Inter-Mount Intervals & Directional Bout Identification
# ------------------------------------------------------------------------------
mounts_df <- setDT(mounts_raw)

# Directional dyad: Actor (Giver) to Recipient (Receiver)
mounts_df[, dyad := paste(ACTOR1, RECIPIENT, sep = "_")]

# Chronological sort within each directional dyad
setorder(mounts_df, dyad, DateTime)

# Time difference in minutes between consecutive mounts by the same actor on the same recipient
mounts_df[, time_diff_mins := as.numeric(difftime(DateTime, shift(DateTime, type = "lag"), units = "mins")), by = dyad]

# Temporal bout threshold (10 minutes)
t_thresh <- 10
mounts_df[, new_bout := (is.na(time_diff_mins) | time_diff_mins > t_thresh)]
mounts_df[, bout_id  := cumsum(new_bout)]

# Summary of mounts per directional bout
bout_summary      <- mounts_df[, .(n_mounts = .N), by = .(dyad, bout_id)]
multi_mount_bouts <- bout_summary[n_mounts > 1]

# Extract records belonging to serial mounting bouts (>= 2 mounts in the bout)
multi_mounts_df <- mounts_df[multi_mount_bouts, on = .(dyad, bout_id), nomatch = 0]

# Exclude initial baseline mount (interval is NA or > 10m from earlier interaction)
multi_mounts_df <- multi_mounts_df[new_bout == FALSE]

# Within-bout average intervals
avg_time_diff_per_bout <- multi_mounts_df[!is.na(time_diff_mins), 
                                          .(avg_time_diff = mean(time_diff_mins), n_intervals = .N), 
                                          by = .(dyad, bout_id)]

cat("==================================================================\n")
cat(" DIRECTIONAL MOUNTING INTERVAL METRICS (10-min Bout Cutoff)\n")
cat("==================================================================\n")
cat("Mean of within-bout averages:    ", round(mean(avg_time_diff_per_bout$avg_time_diff), 3), "mins\n")
cat("Median of within-bout averages:  ", round(median(avg_time_diff_per_bout$avg_time_diff), 3), "mins\n")
cat("Pooled raw interval mean:        ", round(mean(multi_mounts_df$time_diff_mins), 3), "mins\n")
cat("Total serial bouts (>=2 mounts): ", nrow(avg_time_diff_per_bout), "\n\n")

# ------------------------------------------------------------------------------
# 4. Classify DSB vs. SSB Mounts (Directional)
# ------------------------------------------------------------------------------
female_UNs <- c("FUN", "JVF", "JVFUN")
male_UNs   <- c("SAMUN", "MUN", "JVM", "JVMUN")

# Different-Sex Behaviour (DSB): Male actor to Female recipient
mounts_DSB <- mounts_df[ACTOR1 %in% c(all_males, male_UNs) & 
                          RECIPIENT %in% c(all_females, female_UNs)]

# Same-Sex Behaviour (SSB): Male actor to Male recipient
mounts_SSB <- mounts_df[ACTOR1 %in% c(all_males, male_UNs) & 
                          RECIPIENT %in% c(all_males, male_UNs)]

cat("==================================================================\n")
cat(" COHORT MOUNT COUNTS (DIRECTIONAL)\n")
cat("==================================================================\n")
cat("DSB Mounts (Male -> Female): ", nrow(mounts_DSB), "\n")
cat("SSB Mounts (Male -> Male):   ", nrow(mounts_SSB), "\n")
cat("Other / Unclassified:        ", nrow(mounts_df) - (nrow(mounts_DSB) + nrow(mounts_SSB)), "\n\n")

cat("==================================================================\n")
cat(" TIME DIFFERENCE QUANTILES (Consecutive Mounts <= 60 mins)\n")
cat("==================================================================\n")
cat("DSB Quantiles:\n")
print(round(quantile(mounts_DSB[!is.na(time_diff_mins) & time_diff_mins <= 60]$time_diff_mins, probs = seq(0, 1, by = 0.1)), 2))
cat("\nSSB Quantiles:\n")
print(round(quantile(mounts_SSB[!is.na(time_diff_mins) & time_diff_mins <= 60]$time_diff_mins, probs = seq(0, 1, by = 0.1)), 2))

# ------------------------------------------------------------------------------
# 5. Visualisation Panels
# ------------------------------------------------------------------------------
theme_custom <- theme_minimal(base_size = 11) +
  theme(
    plot.title       = element_text(size = 12),
    plot.subtitle    = element_text(color = "gray30", size = 9),
    panel.grid.minor = element_blank()
  )

# Plot 3: DSB specific intervals (<= 60 mins)
p3 <- ggplot(mounts_DSB[!is.na(time_diff_mins) & time_diff_mins <= 60], aes(x = time_diff_mins)) +
  geom_histogram(binwidth = 0.25, color = "black", fill = "skyblue", alpha = 0.85, linewidth = 0.1) +
  coord_cartesian(xlim = c(0, 60)) +
  scale_x_continuous(breaks = seq(0, 60, by = 10)) +
  labs(
    title = "Different-Sex Behaviour (DSB)",
    x     = "Time Difference (minutes)",
    y     = "Count"
  ) +
  theme_custom

# Plot 4: SSB specific intervals (<= 60 mins)
p4 <- ggplot(mounts_SSB[!is.na(time_diff_mins) & time_diff_mins <= 60], aes(x = time_diff_mins)) +
  geom_histogram(binwidth = 0.25, color = "black", fill = "skyblue", alpha = 0.85, linewidth = 0.1) +
  coord_cartesian(xlim = c(0, 60)) +
  scale_x_continuous(breaks = seq(0, 60, by = 10)) +
  labs(
    title = "Same-Sex Behaviour (SSB)",
    x     = "Time Difference (minutes)",
    y     = "Count"
  ) +
  theme_custom

# Combine into a 1x2 panel layout
panel_plot <- (p3 | p4)
print(panel_plot)

# Save to PDF
ggsave(plot_out_pdf, panel_plot, width = 11, height = 5)

cat(sprintf("\nPlot successfully saved to:\n - %s\n", plot_out_pdf))