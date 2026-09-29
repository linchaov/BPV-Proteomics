# ==============================================================================
# METADATA
# ==============================================================================
# Author: Chaoyang Lin, MD
#
# Required packages: data.table, dplyr, fst, survival, broom, writexl
#
# Inputs: data/cohort_olink.fst; tables/supp_bpv_proteins.csv
# Output: tables/supp_protein_outcomes.xlsx (one worksheet per outcome)
#
# Goal: Test the union of BPV-selected proteins against incident outcomes.
# Notes: Each Cox model adjusts for mean SBP and the clinical covariates. 
#        BH-FDR is calculated separately within each outcome across all tested proteins.
# ==============================================================================

library(data.table)
library(dplyr)
library(fst)
library(survival)
library(broom)

output_dir <- "./BPV_ASCVD/tables"
proteins <- fread(file.path(output_dir, "supp_bpv_proteins.csv")) |>
  pull(`Protein name`) |>
  unique()

clinical_covars <- c(
  "Age_i0", "Sex", "Ethnic_White_i0", "Education_College_i0",
  "TownsendDeprivationIndex", "SmokingStatus_i0", "Alcohol_Heavy_i0",
  "PAshort_MV_i0", "BMI_i0", "FamilyHistory_HeartDiseaseORStroke_i0",
  "Baseline_AntiHypertensive_i0"
)
adjusters <- c("mean_SBP_0to5y_i0", clinical_covars)

outcomes <- data.frame(
  outcome = c(
    "Myocardial infarction", "Coronary revascularization",
    "Ischaemic stroke", "Peripheral arterial disease",
    "Major adverse limb event"
  ),
  stem = c("MyocardialInfarction", "CoronaryRevascularization",
           "IschaemicStroke", "PeripheralArterialDisease", "MajorAdverseLimbEvent")
)

cohort <- read_fst(
  "./BPV_ASCVD/data/cohort_olink.fst",
  columns = unique(c(
    proteins, adjusters,
    paste0("Baseline_", outcomes$stem, "_i0"),
    paste0("Incident_", outcomes$stem, "_i0"),
    paste0("Time_", outcomes$stem, "_i0")
  )),
  as.data.table = FALSE
)

# Convert the incident-event field to the 0/1 coding used by Surv().
event01 <- function(x) {
  case_when(
    as.character(x) %in% c("1", "Yes", "TRUE") ~ 1L,
    as.character(x) %in% c("0", "No", "FALSE") ~ 0L,
    TRUE ~ NA_integer_
  )
}

# Fit one adjusted Cox model per protein in the outcome-specific risk set.
fit_protein <- function(protein, risk_set) {
  model_vars <- c(protein, adjusters)
  d <- droplevels(risk_set[complete.cases(risk_set[, model_vars, drop = FALSE]),
                              c(model_vars, ".time", ".event"), drop = FALSE])
  d$.protein_z <- as.numeric(scale(d[[protein]]))
  fit <- coxph(
    reformulate(c(".protein_z", adjusters), response = "Surv(.time, .event)"),
    data = d
  )
  tidy(fit, exponentiate = FALSE, conf.int = TRUE) |>
    filter(term == ".protein_z") |>
    transmute(
      protein = protein, n = fit$n, events = fit$nevent,
      log_hr = estimate, hr = exp(estimate),
      ci_low = exp(conf.low), ci_high = exp(conf.high), p_value = p.value
    )
}

# Apply BH to all screened proteins.
fit_outcome <- function(i) {
  stem <- outcomes$stem[i]
  baseline <- paste0("Baseline_", stem, "_i0")
  event <- paste0("Incident_", stem, "_i0")
  time <- paste0("Time_", stem, "_i0")
  risk_set <- cohort[!is.na(cohort[[baseline]]) &
                       as.character(cohort[[baseline]]) == "No", , drop = FALSE]
  risk_set$.event <- event01(risk_set[[event]])
  risk_set$.time <- suppressWarnings(as.numeric(as.character(risk_set[[time]])))
  risk_set <- risk_set[!is.na(risk_set$.event) &
                         is.finite(risk_set$.time) & risk_set$.time > 0, , drop = FALSE]

  bind_rows(lapply(proteins, fit_protein, risk_set = risk_set)) |>
    mutate(fdr = p.adjust(p_value, method = "BH"), outcome = outcomes$outcome[i])
}

all_results <- bind_rows(lapply(seq_len(nrow(outcomes)), fit_outcome))
fwrite(all_results, file.path(output_dir, "protein_outcome_cox_all.csv"))

tables <- lapply(outcomes$outcome, function(name) {
  all_results |>
    filter(outcome == name, !is.na(fdr), fdr < 0.05) |>
    arrange(fdr, p_value, protein) |>
    transmute(
      `Protein name` = protein,
      `Sample size, n` = n,
      `Events, n` = events,
      `ln(HR)` = sprintf("%.5f", log_hr),
      `HR (95% CI)` = sprintf("%.2f (%.2f-%.2f)", hr, ci_low, ci_high),
      `P-value` = p_value,
      `FDR-adjusted P-value` = fdr
    )
})
names(tables) <- outcomes$outcome
writexl::write_xlsx(tables, file.path(output_dir, "supp_protein_outcomes.xlsx"))
