# ==============================================================================
# METADATA
# ==============================================================================
#
# Author: Chaoyang Lin, MD
# Contact: chaoyanglint@163.com / Chaoyang Lin
#
# Dependencies:
#   R >= 4.1.0
#   Required packages: dplyr, tibble, fst, survival, broom, data.table
#
# Goal:
#   Estimate total-effect Cox associations of SBPV and DBPV with cardiovascular outcomes.
#
# Required input data:
#   - ./BPV_ASCVD/data/cohort_olink.fst
#
# Main outputs:
#   - ./BPV_ASCVD/data/bpv_outcomes.csv
#   - ./BPV_ASCVD/tables/supp_bpv_outcomes.csv
# ==============================================================================

library(survival)
library(data.table)
library(dplyr)

cohort_olink <- fst::read_fst(
  "./BPV_ASCVD/data/cohort_olink.fst",
  as.data.table = FALSE
)

covars <- c(
  "Age_i0", "Sex", "Ethnic_White_i0", "Education_College_i0",
  "TownsendDeprivationIndex", "SmokingStatus_i0", "Alcohol_Heavy_i0",
  "PAshort_MV_i0", "BMI_i0", "FamilyHistory_HeartDiseaseORStroke_i0",
  "Baseline_AntiHypertensive_i0"
)

outcomes <- tibble::tribble(
  ~outcome, ~baseline_var, ~event_var, ~time_var,
  "Myocardial infarction", "Baseline_MyocardialInfarction_i0",
    "Incident_MyocardialInfarction_i0", "Time_MyocardialInfarction_i0",
  "Coronary revascularization", "Baseline_CoronaryRevascularization_i0",
    "Incident_CoronaryRevascularization_i0", "Time_CoronaryRevascularization_i0",
  "Ischaemic stroke", "Baseline_IschaemicStroke_i0",
    "Incident_IschaemicStroke_i0", "Time_IschaemicStroke_i0",
  "Peripheral arterial disease", "Baseline_PeripheralArterialDisease_i0",
    "Incident_PeripheralArterialDisease_i0", "Time_PeripheralArterialDisease_i0",
  "Major adverse limb event", "Baseline_MajorAdverseLimbEvent_i0",
    "Incident_MajorAdverseLimbEvent_i0", "Time_MajorAdverseLimbEvent_i0"
)

# Convert incident-event flags to 0/1, leaving unrecognized values missing.
to_total_event01 <- function(x) {
  x <- as.character(x)
  case_when(
    x %in% c("1", "Yes", "TRUE") ~ 1L,
    x %in% c("0", "No", "FALSE") ~ 0L,
    TRUE ~ NA_integer_
  )
}

# Display P values to three decimals without printing small values as zero.
format_table_p <- function(p) {
  ifelse(is.na(p), NA_character_, ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))
}

# Fit the adjusted BPV-outcome Cox models on each outcome's event-free cohort.
run_bpv_total_effect <- function(data, trait, exposure_var, mean_bp_var) {
  adjustment_vars <- c(exposure_var, mean_bp_var, covars)

  bind_rows(lapply(seq_len(nrow(outcomes)), function(i) {
    o <- outcomes[i, ]
    vars <- unique(c(adjustment_vars, o$baseline_var, o$event_var, o$time_var))
    d <- data[, vars, drop = FALSE]
    d <- d[!is.na(d[[o$baseline_var]]) & as.character(d[[o$baseline_var]]) == "No", , drop = FALSE]
    d$.event <- to_total_event01(d[[o$event_var]])
    d$.time <- suppressWarnings(as.numeric(as.character(d[[o$time_var]])))
    keep <- complete.cases(d[, adjustment_vars, drop = FALSE]) &
      !is.na(d$.event) & is.finite(d$.time) & d$.time > 0
    d <- droplevels(d[keep, , drop = FALSE])

    # The exposure is modeled separately for each BPV trait. cox.zph tests
    # its Schoenfeld-residual trend and the full model's PH assumption.
    fit <- coxph(
      reformulate(adjustment_vars, response = "Surv(.time, .event)"),
      data = d
    )
    ph <- cox.zph(fit, transform = "km")$table
    person_years <- sum(d$.time) / 365.25
    incidence_ci <- poisson.test(fit$nevent, T = person_years)$conf.int * 1000
    broom::tidy(fit, exponentiate = TRUE, conf.int = TRUE) |>
      filter(term == exposure_var) |>
      transmute(
        trait = .env$trait, outcome = o$outcome, n = fit$n, events = fit$nevent,
        hr = estimate, ci_low = conf.low, ci_high = conf.high, p_value = p.value,
        ph_p_bpv = ph[exposure_var, "p"], ph_p_global = ph["GLOBAL", "p"],
        person_years = person_years,
        incidence_rate = fit$nevent / person_years * 1000,
        incidence_ci_low = incidence_ci[1], incidence_ci_high = incidence_ci[2]
      )
  }))
}

# One-at-a-time BPV models; BH-FDR is applied separately within each trait.
bpv_total_effect_cox <- bind_rows(
  run_bpv_total_effect(cohort_olink, "SBPV", "SBPV_VIM_z", "mean_SBP_0to5y_i0"),
  run_bpv_total_effect(cohort_olink, "DBPV", "DBPV_VIM_z", "mean_DBP_0to5y_i0")
)

fwrite(
  bpv_total_effect_cox |>
    arrange(match(outcome, outcomes$outcome), match(trait, c("SBPV", "DBPV"))) |>
    transmute(
      Outcome = outcome,
      Exposure = paste0(trait, "-VIM (per 1 SD)"),
      `N analyzed` = n,
      `Events / person-years` = sprintf("%d / %.0f", events, person_years),
      `Incidence rate (95% CI) per 1,000 person-years` = sprintf(
        "%.2f (%.2f-%.2f)", incidence_rate, incidence_ci_low, incidence_ci_high
      ),
      `HR (95% CI)` = sprintf("%.2f (%.2f-%.2f)", hr, ci_low, ci_high),
      `P-value` = format_table_p(p_value),
      `BPV PH P-value` = format_table_p(ph_p_bpv),
      `Global PH P-value` = format_table_p(ph_p_global)
    ),
  "./BPV_ASCVD/tables/supp_bpv_outcomes.csv"
)

