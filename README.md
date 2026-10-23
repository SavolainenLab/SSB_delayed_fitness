# Male Same-Sex Sexual Behaviour (SSB) and Delayed Reproductive Success in Rhesus Macaques

This repository contains the curated datasets, Bayesian statistical modelling scripts, data wrangling, and visualisation workflows supporting the manuscript:

> **"Delayed fitness benefits of same-sex sexual behaviour in a nonhuman primate species"**

---

## 1. System Requirements

### Operating Systems
* **Tested platform:** macOS (Tahoe, Apple Silicon arm64) running R version 4.6.1.
* **Compatibility:** Fully compatible with Linux and Microsoft Windows 10/11.

### Software Dependencies
The analytical pipeline runs in **R (version >= 4.4.2)**. An active C++ compiler toolchain is required by Stan/brms to compile Bayesian models:

### Core R Packages and Verified Versions
* **Data manipulation and utilities:**
  * `tidyverse` (v2.0.0; including `dplyr` v1.2.1, `ggplot2` v4.0.3, `readr` v2.2.0, `tidyr` v1.3.2, `purrr` v1.2.2, `stringr` v1.6.0, `lubridate` v1.9.5)
  * `readxl` (v1.4.3), `openxlsx` (v4.2.5), `writexl` (v1.5.0)
  * `future` (v1.33.0), `qpdf` (v1.3.2), `magick` (v2.8.0), `pdftools` (v3.4.0)
* **Bayesian modelling, network inference, and diagnostics:**
  * `brms` (v2.23.0)
  * `posterior` (v1.7.0), `tidybayes` (v3.0.6), `bayesplot` (v1.16.0), `bayestestR` (v0.13.2), `coda` (v0.19.4)
  * `bisonR` (v0.1.0; GitHub: `robayed/bisonR`), `igraph` (v1.5.1)
  * `cowplot` (v1.1.1), `patchwork` (v1.2.0), `ggstance` (v0.3.5)

### Non-Standard Hardware
* **Hardware requirements:** No specialised hardware is required.
* **Recommended specifications:** A standard desktop computer or laptop with a multi-core processor (>= 4 physical CPU cores) and at least 16 GB RAM is recommended to support multi-chain MCMC sampling in parallel.

---

## 2. Installation Guide

### Installation from GitHub
Clone the repository to your local computer:

```bash
git clone https://github.com/SavolainenLab/SSB_delayed_fitness.git
cd SSB_delayed_fitness
