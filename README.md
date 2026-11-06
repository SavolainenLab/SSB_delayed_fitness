# Male Same-Sex Sexual Behaviour (SSB) and Delayed Reproductive Success in Rhesus Macaques

This repository contains the curated datasets, Bayesian statistical modelling scripts, data wrangling pipelines, and visualisation workflows supporting the manuscript:

> **"Delayed fitness benefits of same-sex sexual behaviour in a nonhuman primate species"**

---

## 1. System Requirements

### Operating Systems
- **Tested platform:** macOS (Tahoe, Apple Silicon arm64) running R version 4.6.1.
- **Compatibility:** Fully compatible with Linux and Microsoft Windows 10/11.

### Software Dependencies
The analytical pipeline runs in **R (version >= 4.4.2)**. An active C++ compiler toolchain is required by Stan and brms to compile Bayesian models.

### Core R Packages and Verified Versions
- **Data manipulation and utilities:**
  * `tidyverse` (v2.0.0; including `dplyr` v1.2.1, `ggplot2` v4.0.3, `readr` v2.2.0, `tidyr` v1.3.2, `purrr` v1.2.2, `stringr` v1.6.0, `lubridate` v1.9.5)
  * `readxl` (v1.4.3), `openxlsx` (v4.2.5), `writexl` (v1.5.0), `janitor` (v2.2.1), `data.table` (v1.18.6)
  * `future` (v1.33.0), `qpdf` (v1.3.2), `magick` (v2.8.0), `pdftools` (v3.4.0)
- **Dominance rank and network analysis:**
  * `EloRating` (v0.46.11), `bisonR` (v0.1.0; GitHub: `robayed/bisonR`), `igraph` (v1.5.1)
- **Bayesian modelling and diagnostics:**
  * `brms` (v2.23.0), `rstan` (v2.32.7)
  * `posterior` (v1.7.0), `tidybayes` (v3.0.6), `bayesplot` (v1.16.0), `bayestestR` (v0.13.2), `coda` (v0.19.4)
  * `cowplot` (v1.1.1), `patchwork` (v1.2.0), `ggstance` (v0.3.5)

### Non-Standard Hardware
- **Hardware requirements:** No specialised hardware is required.
- **Recommended specifications:** A standard desktop computer or laptop with a multi-core processor (>= 4 physical CPU cores) and at least 16 GB RAM is recommended to support multi-chain MCMC sampling in parallel.

---

## 2. Installation Guide

### Installation from GitHub
Clone the repository to your local computer:

```bash
git clone [https://github.com/SavolainenLab/SSB_delayed_fitness.git](https://github.com/SavolainenLab/SSB_delayed_fitness.git)
cd SSB_delayed_fitness
```

### Installing Package Dependencies
Open an R session from the project directory and install the required packages:

```R
# 1. Install CRAN dependencies
cran_packages <- c(
  "tidyverse", "readxl", "writexl", "openxlsx", "brms", "rstan",
  "tidybayes", "bayesplot", "posterior", "bayestestR", "coda",
  "cowplot", "patchwork", "ggstance", "igraph", "future", "qpdf",
  "broom.mixed", "lubridate", "EloRating", "janitor", "data.table"
)
install.packages(cran_packages)

# 2. Install GitHub dependencies (bisonR for network analysis)
if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
remotes::install_github("robayed/bisonR")
```

### Typical Install Time
- Typical install time is less than 15 minutes on a standard desktop computer with an active internet connection.

---

## 3. Instructions for Use

### Data Organisation
The primary data source is the collated supplementary workbook:
- `Supplementary_Tables.xlsx`: Contains all supplementary tables referenced in the main text. Table descriptions are provided in the sheet `'Table Index'`.
- All animal identifiers have been replaced with unique four-letter anonymised codes. This preserves relational integrity across behavioural observations, dyadic networks, and pedigree records while maintaining privacy.

### Data Generation Pipeline
All processed individual and dyadic analysis datasets can be fully generated from scratch using only the anonymised raw tables in `Supplementary_Tables.xlsx`:

- **Execution Script:** `Data_cleaning_pipeline/Combine_cleaned_data_v1.4.R`

Executing this singular script runs the modular wrangling scripts in sequence within isolated environments:
1. `Scans/Generate_scans_per_year_v1.2.R` (Individual observation scan effort)
2. `Mounts_and_Age/Generate_behav_data_from_master_repeated_measures.R` (Mount counts and chronological age)
3. `Rank/Generate_rank_David_v1.2_group_by_observer.R` (David's score dominance rank calculations)
4. `Fecundity/Calculate_annual_fecundity_v1.3.R` (Annual siring success)
5. `Death_date/Death_dates.R` (Demographic exit dates)
6. `Dyadic/Generate_dyadic_data.R` (Dyadic behavioral and association matrices)

The script collates all modular outputs, resolves social group assignments, computes future reproductive output (t+1 to t+6) and future standardized rank (t+1 to t+8), and exports the final analysis files to `Data/Cleaned_outputs/`:
- `Combined_cleaned_collapsed_w_OBS.csv` (Used in downstream analyses; Focal male cohort; 5-minute bout-collapsed mounts)
- `Combined_cleaned_uncollapsed_w_OBS.csv` (Focal male cohort; uncollapsed mounts)
- `Dyadic_counts.csv` (Within-group dyadic behavioral and attribute dataset)

### Running the Analysis
1. Open the repository directory in R or RStudio, or set your working directory to the project root:

```R
setwd("path/to/SSB_delayed_fitness")
```

2. Execute the self-contained R scripts corresponding to each section of the manuscript:

| Manuscript Section | Primary Script | Description and Outputs |
| :--- | :--- | :--- |
| **Immediate Reproductive Output** | `Male_annual_fitness_SSB_v1.2_publication_ready.R` | Multivariate GLMM assessing annual associations between SSB, DSB, and offspring sired; outputs parameter summaries and Figure 1 panels. |
| **Dominance Rank Associations** | `Multivariate_Rank_SSB_DSB_v1.2_publication_ready.R` | Multivariate model evaluating age trajectories for SSB, DSB and dominance rank; generates Figure 2 panels. |
| **Affiliative Male-Male Bonds** | `Social_network_analysis_publication_ready_v1.2.R` | BISoN edge-weight models and nodal regressions for grooming, coalition, and aggression networks; outputs network metrics and Figure 3a panel. |
| **Delayed Fitness and Rank Ascent** | `Future_repro_rank_v1.4.8_noDS_no_mi_death_accounted_pub_ready.R` | Longitudinal models testing delayed siring success ($t+1$ to $t+4$) and rank ascent, including model comparisons; generates Figure 3 & 4 panels. |
| **Observability Sensitivity Analysis** | `Future_repro_rank_v1.4.5.2_brm_multiple_mi_death_accounted_pub_ready.R` | Multiple-imputation models evaluating robustness of delayed fitness associations accounting for variation in male observability. |

### Computational Execution Notes
- **Fitting models from scratch:** Full Bayesian mixed-effects models require substantial computational time (several hours depending on processor speed and allocated CPU cores).

---

## 4. Citation and Licence

### Citation
If you use the code or data from this repository, please cite:

Bell-Roberts, L. et al. Delayed fitness benefits of same-sex sexual behaviour in a nonhuman primate species (under review).

### Licence
- **Code Licence:** MIT Licence
- **Data Licence:** Creative Commons Attribution 4.0 International Licence (CC-BY 4.0)
