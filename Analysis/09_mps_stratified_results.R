# ==============================================================================
# METADATA
# ==============================================================================
# Author: Chaoyang Lin, MD
# Dependencies: R >= 4.2.0; fst, survival, compareC, nricens, EValue,
#               future.apply, future, writexl
# Goal: Evaluate alpha-0.9 elastic-net MPS associations and apparent
#       incremental value for both BPV traits across five outcomes.
# Strata: all; hypertension yes/no; antihypertensive or statin yes/no.
# Outputs: prediction/mps_stratified/mps_results.csv
#          and mps_tables.xlsx.
# Notes: Fine-Gray uses all-cause death before the target outcome as the
#        competing event. C is Harrell's C; NRI is category-free at 10 years.
#        Interaction P values come from fitted MPS-by-group terms in the cohort. 
#       Runtime: 500 NRI bootstrap samples.
# ==============================================================================

library(fst)
library(survival)
library(compareC)
library(nricens)
library(future)
library(future.apply)
library(writexl)

data_root <- "./BPV_ASCVD"
output_dir <- file.path(data_root, "prediction/mps_stratified")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

horizon <- 10 * 365.25

clinical_covars <- c(
  "Age_i0", "Sex", "Ethnic_White_i0", "Education_College_i0",
  "TownsendDeprivationIndex", "SmokingStatus_i0", "Alcohol_Heavy_i0",
  "PAshort_MV_i0", "BMI_i0", "FamilyHistory_HeartDiseaseORStroke_i0",
  "Baseline_AntiHypertensive_i0"
)
outcomes <- data.frame(
  Outcome = c("MI", "Coronary revascularization", "Ischemic stroke", "PAD", "MALE"),
  stem = c("MyocardialInfarction", "CoronaryRevascularization",
           "IschaemicStroke", "PeripheralArterialDisease", "MajorAdverseLimbEvent")
)
strata <- c("All", "Hypertension", "No_hypertension",
            "Antihypertensive_or_statin", "Neither_medication")

outcome_columns <- unlist(lapply(c("Baseline_", "Incident_", "Time_"),
  function(prefix) paste0(prefix, outcomes$stem, "_i0")))
cohort <- read_fst(file.path(data_root, "data/cohort_olink.fst"),
  columns = unique(c("ID", clinical_covars, "Baseline_Hypertension_i0",
    "Baseline_LipidLowering_Statin_i0", "Death_Any_AllCause_i0",
    "Time_Death_Summary_i0", "SBPV_VIM_z", "DBPV_VIM_z",
    "mean_SBP_0to5y_i0", "mean_DBP_0to5y_i0", outcome_columns)),
  as.data.table = FALSE)
scores <- read_fst(file.path(data_root, "data/overall_mps/mps_scores.fst"),
                   as.data.table = FALSE)
data <- cbind(cohort, scores[match(cohort$ID, scores$ID), -1L, drop = FALSE])

# Convert Yes/No or 1/0 fields to numeric event indicators; keep missing values.
event01 <- function(x) {
  x <- as.character(x)
  ifelse(x %in% c("Yes", "1", "TRUE"), 1L,
         ifelse(x %in% c("No", "0", "FALSE"), 0L, NA_integer_))
}

# Mark use of either antihypertensive medication or statins.
medication_group <- function(d) {
  anti <- event01(d$Baseline_AntiHypertensive_i0)
  statin <- event01(d$Baseline_LipidLowering_Statin_i0)
  ifelse(anti == 1L | statin == 1L, 1L,
         ifelse(anti == 0L & statin == 0L, 0L, NA_integer_))
}

# Select the cohort or one of the hypertension and medication strata.
make_stratum <- function(d, group) {
  htn <- event01(d$Baseline_Hypertension_i0)
  medication <- medication_group(d)
  switch(group,
    All = rep(TRUE, nrow(d)),
    Hypertension = htn == 1L,
    No_hypertension = htn == 0L,
    Antihypertensive_or_statin = medication == 1L,
    Neither_medication = medication == 0L
  )
}

# Convert a fitted Cox model to individual 10-year risks for the NRI analysis.
risk_10y <- function(fit) {
  baseline <- survival::basehaz(fit, centered = TRUE)
  idx <- findInterval(horizon, baseline$time)
  h0 <- if (idx == 0L) 0 else baseline$hazard[idx]
  1 - exp(-h0 * exp(fit$linear.predictors))
}

# Build an outcome-specific complete-case sample with BPV, MPS and death times.
prepare_model_data <- function(trait, outcome, group) {
  exposure <- paste0(trait, "_VIM_z")
  score <- paste0(trait, "_MPS")
  mean_bp <- if (trait == "SBPV") "mean_SBP_0to5y_i0" else "mean_DBP_0to5y_i0"
  baseline <- paste0("Baseline_", outcome$stem, "_i0")
  event <- paste0("Incident_", outcome$stem, "_i0")
  time <- paste0("Time_", outcome$stem, "_i0")
  adjusters <- c(mean_bp, exposure, clinical_covars)
  # Include antihypertensive medication only in the full-cohort model.
  if (group != "All") {
    adjusters <- setdiff(adjusters, "Baseline_AntiHypertensive_i0")
  }

  included <- make_stratum(data, group) &
    as.character(data[[baseline]]) %in% c("No", "0", "FALSE")
  d <- data[which(included), unique(c("ID", score, adjusters, event, time,
    "Death_Any_AllCause_i0", "Time_Death_Summary_i0")), drop = FALSE]
  d$.event <- event01(d[[event]])
  d$.time <- as.numeric(as.character(d[[time]]))
  d$.death <- event01(d$Death_Any_AllCause_i0)
  d$.death_time <- as.numeric(as.character(d$Time_Death_Summary_i0))
  d <- droplevels(d[complete.cases(d[, c(score, adjusters, ".event", ".time",
    ".death", ".death_time")]) & d$.time > 0 & d$.death_time > 0, , drop = FALSE])
  # Remove covariates constant within a stratum.
  adjusters <- adjusters[vapply(d[, adjusters, drop = FALSE],
    function(x) length(unique(x)) > 1L, logical(1))]
  list(data = d, score = score, adjusters = adjusters)
}

# Expand follow-up into Fine-Gray risk intervals with all-cause death competing.
make_finegray_data <- function(d, score, adjusters) {
  # An outcome recorded on the death date takes priority over death.
  death_first <- d$.death == 1L & d$.death_time <= d$.time &
    (d$.event == 0L | d$.death_time < d$.time)
  d$.fg_time <- ifelse(death_first, d$.death_time, d$.time)
  d$.fg_status <- factor(ifelse(death_first, "death",
    ifelse(d$.event == 1L, "target", "censor")),
    levels = c("censor", "target", "death"))
  fg <- survival::finegray(
    reformulate(c(score, adjusters, "ID"),
                response = "Surv(.fg_time, .fg_status)"),
    data = d, etype = "target"
  )
  list(data = fg, competing_deaths = sum(death_first))
}

# Fit the nested Cox models and the weighted Fine-Gray model in one stratum.
# Return association estimates and apparent incremental-value statistics.
run_model <- function(trait, outcome, group, seed) {
  prepared <- prepare_model_data(trait, outcome, group)
  d <- prepared$data
  score <- prepared$score
  adjusters <- prepared$adjusters

  base <- coxph(reformulate(adjusters, response = "Surv(.time, .event)"),
                data = d, ties = "efron", x = TRUE)
  full <- coxph(reformulate(c(adjusters, score), response = "Surv(.time, .event)"),
                data = d, ties = "efron", x = TRUE)
  fg <- make_finegray_data(d, score, adjusters)
  fg_fit <- coxph(reformulate(c(adjusters, score),
                              response = "Surv(fgstart, fgstop, fgstatus)"),
                  data = fg$data, weights = fgwt, cluster = ID, ties = "efron")
  fg_b <- coef(fg_fit)[[score]]
  fg_se <- sqrt(vcov(fg_fit)[score, score])
  fg_ci <- exp(fg_b + c(-1, 1) * qnorm(0.975) * fg_se)

  # compareC scores increase with survival; Cox linear predictors increase with risk.
  c_test <- compareC(d$.time, d$.event,
                     -full$linear.predictors, -base$linear.predictors)
  ll0 <- logLik(base)
  ll1 <- logLik(full)
  lrt_chi2 <- 2 * (as.numeric(ll1) - as.numeric(ll0))
  lrt_df <- attr(ll1, "df") - attr(ll0, "df")

  set.seed(seed)
  invisible(capture.output({
    nri_fit <- suppressMessages(nricens::nricens(
      time = d$.time, event = d$.event,
      p.std = risk_10y(base), p.new = risk_10y(full),
      t0 = horizon, updown = "diff", cut = 0,
      point.method = "km", niter = 500L, msg = FALSE
    ))
  }))
  nri <- nri_fit$nri
  boot_se <- sd(nri_fit$bootstrapsample[, "NRI"], na.rm = TRUE)
  nri_p <- if (is.finite(boot_se) && boot_se > 0) {
    2 * pnorm(-abs(nri["NRI", "Estimate"] / boot_se))
  } else NA_real_

  data.frame(
    Trait = trait, Outcome = outcome$Outcome, Stratum = group,
    Events = sum(d$.event), Competing_deaths = fg$competing_deaths,
    MPS_SHR = exp(fg_b), MPS_SHR_CI_low = fg_ci[1],
    MPS_SHR_CI_high = fg_ci[2],
    MPS_SHR_P = summary(fg_fit)$coefficients[score, "Pr(>|z|)"],
    C_reference = c_test$est.c[2], C_with_MPS = c_test$est.c[1],
    Delta_C = c_test$est.diff_c, Delta_C_P = c_test$pval,
    NRI_10y = nri["NRI", "Estimate"],
    NRI_10y_P_bootstrap_Wald = nri_p,
    NRI_events = nri["NRI+", "Estimate"],
    NRI_nonevents = nri["NRI-", "Estimate"],
    LRT_chi_square = lrt_chi2,
    LRT_P = pchisq(lrt_chi2, df = lrt_df, lower.tail = FALSE)
  )
}

jobs <- list()
for (i in seq_len(nrow(outcomes))) {
  for (group in strata) for (trait in c("SBPV", "DBPV")) {
    jobs[[length(jobs) + 1L]] <- list(
      trait = trait, outcome = outcomes[i, ], group = group
    )
  }
}
future::plan(future::multisession, workers = 6L)
results <- future.apply::future_lapply(seq_along(jobs), function(i) {
  job <- jobs[[i]]
  run_model(job$trait, job$outcome, job$group, 20261001L + i)
}, future.seed = TRUE, future.packages = c("survival", "compareC", "nricens"),
   future.chunk.size = 1)
future::plan(future::sequential)

result_table <- do.call(rbind, results)

# Arrange each outcome as All (SBPV, DBPV), then each subgroup (SBPV, DBPV).
result_table <- result_table[order(
  match(result_table$Outcome, outcomes$Outcome),
  match(result_table$Stratum, strata),
  match(result_table$Trait, c("SBPV", "DBPV"))
), , drop = FALSE]
row.names(result_table) <- NULL
write.csv(result_table, file.path(output_dir, "mps_results.csv"),
          row.names = FALSE, na = "")

# Mark P values in the compact main table; keep numeric P values elsewhere.
format_significance <- function(p) {
  ifelse(is.na(p), "", ifelse(p < 0.001, "***",
    ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "ns"))))
}

# Format Cox HR and Fine-Gray SHR with two-decimal confidence limits.
format_ratio <- function(est, low, high) {
  sprintf("%.2f (%.2f\u2013%.2f)", est, low, high)
}

# Show the C-index change with a significance symbol in the compact table.
format_c_symbols <- function(d) {
  arrow <- ifelse(d$Delta_C >= 0, "\u2191", "\u2193")
  sprintf("%.3f (%s%.3f %s)", d$C_with_MPS, arrow,
          abs(d$Delta_C), format_significance(d$Delta_C_P))
}

# Replace internal stratum codes with table-ready labels.
display_stratum <- function(x) {
  unname(c(
    All = "All", Hypertension = "Hypertension",
    No_hypertension = "No hypertension",
    Antihypertensive_or_statin = "Antihypertensive or statin",
    Neither_medication = "Neither medication"
  )[x])
}

# Combine effect estimates and significance symbols in a compact main table.
format_main_symbols <- function(d) {
  data.frame(
    Outcome = d$Outcome, Stratum = display_stratum(d$Stratum),
    MPS = paste0("MPS-", d$Trait), Events = d$Events,
    `Competing deaths` = d$Competing_deaths,
    `SHR (95% CI), significance` = paste0(
      format_ratio(d$MPS_SHR, d$MPS_SHR_CI_low, d$MPS_SHR_CI_high),
      " ", format_significance(d$MPS_SHR_P)),
    `Cox C-index (change), significance` = format_c_symbols(d),
    `Cox LRT chi-square` = paste0(
      sprintf("%.1f", d$LRT_chi_square), " ", format_significance(d$LRT_P)),
    `NRI events` = sprintf("%.1f%%", 100 * d$NRI_events),
    `NRI nonevents` = sprintf("%.1f%%", 100 * d$NRI_nonevents),
    `Total NRI, significance` = paste0(
      sprintf("%.1f%%", 100 * d$NRI_10y), " ",
      format_significance(d$NRI_10y_P_bootstrap_Wald)),
    check.names = FALSE
  )
}

writexl::write_xlsx(
  list(Table3 = format_main_symbols(result_table)),
  path = file.path(output_dir, "mps_tables.xlsx")
)
