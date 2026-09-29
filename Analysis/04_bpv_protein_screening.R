# ==============================================================================
# METADATA
# ==============================================================================
#
# Author: Chaoyang Lin, MD
# Created: 2026-09-19
# Last modified: 2026-09-20
# Contact: chaoyanglint@163.com / Chaoyang Lin
#
# Dependencies:
#   R >= 4.1.0
#   Required packages: data.table, dplyr, tidyr, tibble, fst, limma
#
# Goal:
#   Estimate BPV-protein associations and apply the prespecified protein screen.
#
# Required input data:
#   - ./BPV_ASCVD/data/cohort_olink.fst
#   - ./BPV_ASCVD/data/olink_proteins.rds
#
# Main outputs:
#   - ./BPV_ASCVD/data/bpv_protein_results.csv
#   - ./BPV_ASCVD/tables/supp_bpv_proteins.csv
#   - ./BPV_ASCVD/data/bpv_protein_plot_data.rds
#
# Notes:
# Selected protein: 5-year VIM FDR <0.05; consistent direction across 5-year VIM/SD/CV/ARV and 10-year VIM;
# 10-year VIM FDR <0.05; and at least one of 5-year SD/CV/ARV FDR <0.05.
# Each model adjusts for base_covars and the corresponding mean BP
# ==============================================================================

library(data.table)
library(tidyr)
library(limma)
library(dplyr)

cohort_olink <- fst::read_fst(
  "./BPV_ASCVD/data/cohort_olink.fst",
  as.data.table = FALSE
)
olink_cols <- readRDS("./BPV_ASCVD/data/olink_proteins.rds")

base_covars <- c(
  "Age_i0", "Sex", "Ethnic_White_i0", "Education_College_i0",
  "TownsendDeprivationIndex", "SmokingStatus_i0", "Alcohol_Heavy_i0",
  "PAshort_MV_i0", "BMI_i0", "FamilyHistory_HeartDiseaseORStroke_i0",
  "Baseline_AntiHypertensive_i0"
)

bpv_specs <- tibble::tribble(
  ~spec_id, ~pressure, ~window, ~metric, ~exposure_var, ~mean_bp_var,
  "SBPV_5y_VIM",  "SBPV", "5y",  "VIM", "SBPV_VIM_z",      "mean_SBP_0to5y_i0",
  "SBPV_5y_SD",   "SBPV", "5y",  "SD",  "SBPV_SD_z",       "mean_SBP_0to5y_i0",
  "SBPV_5y_CV",   "SBPV", "5y",  "CV",  "SBPV_CV_z",       "mean_SBP_0to5y_i0",
  "SBPV_5y_ARV",  "SBPV", "5y",  "ARV", "SBPV_ARV_z",      "mean_SBP_0to5y_i0",
  "SBPV_10y_VIM", "SBPV", "10y", "VIM", "SBPV_VIM_10y_z",  "mean_SBP_0to10y_i0",
  "DBPV_5y_VIM",  "DBPV", "5y",  "VIM", "DBPV_VIM_z",      "mean_DBP_0to5y_i0",
  "DBPV_5y_SD",   "DBPV", "5y",  "SD",  "DBPV_SD_z",       "mean_DBP_0to5y_i0",
  "DBPV_5y_CV",   "DBPV", "5y",  "CV",  "DBPV_CV_z",       "mean_DBP_0to5y_i0",
  "DBPV_5y_ARV",  "DBPV", "5y",  "ARV", "DBPV_ARV_z",      "mean_DBP_0to5y_i0",
  "DBPV_10y_VIM", "DBPV", "10y", "VIM", "DBPV_VIM_10y_z",  "mean_DBP_0to10y_i0"
)

# Fit one prespecified BPV metric against all available protein assays with limma.
run_limma_pwas <- function(data, proteins, spec) {
  exposure_var <- spec$exposure_var[[1]]
  mean_bp_var <- spec$mean_bp_var[[1]]
  model_vars <- unique(c(exposure_var, mean_bp_var, base_covars))

  keep <- complete.cases(data[, model_vars, drop = FALSE])
  model_data <- droplevels(data[keep, , drop = FALSE])

  protein_raw <- as.matrix(model_data[, proteins, drop = FALSE])
  protein_sd <- apply(protein_raw, 2, sd, na.rm = TRUE)
  protein_raw <- protein_raw[, is.finite(protein_sd) & protein_sd > 0, drop = FALSE]
  # Standardize each assay across the complete-case model sample; limma
  # estimates the protein change per 1-SD increment in the BPV metric.
  protein_z <- t(scale(protein_raw))

  design <- model.matrix(
    reformulate(c(exposure_var, mean_bp_var, base_covars)),
    data = model_data
  )
  fit <- eBayes(lmFit(protein_z, design))
  tab <- topTable(
    fit,
    coef = exposure_var,
    number = Inf,
    adjust.method = "BH",
    sort.by = "P",
    confint = TRUE
  )
  protein_n <- rowSums(!is.na(protein_z))

  data.frame(
    feature = rownames(tab),
    beta = tab$logFC,
    ci_low = tab$CI.L,
    ci_high = tab$CI.R,
    p_value = tab$P.Value,
    fdr = tab$adj.P.Val,
    n_nonmissing = unname(protein_n[rownames(tab)]),
    n_analysis = nrow(model_data),
    spec_id = spec$spec_id[[1]],
    pressure = spec$pressure[[1]],
    window = spec$window[[1]],
    metric = spec$metric[[1]],
    stringsAsFactors = FALSE
  )
}

# Retain proteins meeting the five-year, ten-year, alternative-metric, and direction rules.
make_selection_table <- function(pwas_long, pressure) {
  wide <- pwas_long |>
    filter(.data$pressure == .env$pressure) |>
    select(feature, spec_id, beta, fdr) |>
    pivot_wider(
      names_from = spec_id,
      values_from = c(beta, fdr),
      names_glue = "{.value}_{spec_id}"
    ) |>
    as.data.frame()

  model_ids <- paste0(
    pressure,
    c("_5y_VIM", "_5y_SD", "_5y_CV", "_5y_ARV", "_10y_VIM")
  )
  beta_cols <- paste0("beta_", model_ids)
  alternative_fdr_cols <- paste0(
    "fdr_", pressure, c("_5y_SD", "_5y_CV", "_5y_ARV")
  )

  # 5-year VIM and 10-year VIM must pass FDR, all five signs agree, and
  # at least one of 5-year SD/CV/ARV result pass FDR.
  complete_direction <- complete.cases(wide[, beta_cols, drop = FALSE])
  same_direction <- (
    rowSums(wide[, beta_cols, drop = FALSE] > 0, na.rm = TRUE) == 5L |
      rowSums(wide[, beta_cols, drop = FALSE] < 0, na.rm = TRUE) == 5L
  )

  data.frame(
    pressure = pressure,
    feature = wide$feature,
    passes_protein_screen =
      !is.na(wide[[paste0("fdr_", pressure, "_5y_VIM")]]) &
      wide[[paste0("fdr_", pressure, "_5y_VIM")]] < 0.05 &
      complete_direction & same_direction &
      !is.na(wide[[paste0("fdr_", pressure, "_10y_VIM")]]) &
      wide[[paste0("fdr_", pressure, "_10y_VIM")]] < 0.05 &
      rowSums(wide[, alternative_fdr_cols, drop = FALSE] < 0.05,
              na.rm = TRUE) >= 1L
  )
}

pwas_all <- bind_rows(lapply(seq_len(nrow(bpv_specs)), function(i) {
  run_limma_pwas(cohort_olink, olink_cols, bpv_specs[i, , drop = FALSE])
}))

protein_selection <- bind_rows(
  make_selection_table(pwas_all, "SBPV"),
  make_selection_table(pwas_all, "DBPV")
)
pwas_all <- pwas_all |>
  left_join(
    protein_selection, by = c("pressure", "feature"),
    relationship = "many-to-one"
  )

selected_proteins <- pwas_all |>
  filter(
    spec_id %in% c("SBPV_5y_VIM", "DBPV_5y_VIM"),
    passes_protein_screen
  ) |>
  select(-passes_protein_screen)

supplementary_proteins <- selected_proteins |>
  transmute(
    `Protein name` = feature,
    Beta = beta,
    `95% CI lower` = ci_low,
    `95% CI upper` = ci_high,
    `P-value` = p_value,
    `FDR-adjusted P-value` = fdr,
    `N available` = n_nonmissing,
    `N analyzed` = n_analysis,
    Exposure = paste0(pressure, "-VIM Z-Score"),
    `Analysis method` = "limma_PWAS"
  )

fwrite(
  selected_proteins,
  "./BPV_ASCVD/data/bpv_protein_results.csv"
)
fwrite(
  supplementary_proteins,
  "./BPV_ASCVD/tables/supp_bpv_proteins.csv"
)
saveRDS(
  pwas_all |>
    filter(spec_id %in% c("SBPV_5y_VIM", "DBPV_5y_VIM")) |>
    select(feature, pressure, beta, p_value, fdr, passes_protein_screen),
  "./BPV_ASCVD/data/bpv_protein_plot_data.rds"
)
