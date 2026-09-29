# ==============================================================================
# METADATA
# ==============================================================================
# Author: Chaoyang Lin, MD
# Dependencies: R >= 4.1.0; dplyr, tibble, fst, data.table, glmnet,
#               missMDA, survival, broom
# Goal: Build alpha-0.9 elastic-net overall MPSs and test outcome associations.
# Inputs: data/cohort_olink.fst; data/bpv_protein_results.csv.
# Outputs: data/overall_mps/mps_scores.fst;
#          mps_coefficients.csv; mps_model_summary.csv;
#          mps_formulas.csv and .txt; mps_outcomes.csv in the same directory.
# ==============================================================================

library(glmnet)
library(missMDA)
library(survival)
library(data.table)
library(dplyr)

root <- "./BPV_ASCVD"
output_dir <- file.path(root, "data", "overall_mps")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
cohort <- fst::read_fst(file.path(root, "data", "cohort_olink.fst"),
                        as.data.table = FALSE)
selected <- fread(file.path(root, "data", "bpv_protein_results.csv")) |>
  filter(spec_id %in% c("SBPV_5y_VIM", "DBPV_5y_VIM"))

clinical_covars <- c(
  "Age_i0", "Sex", "Ethnic_White_i0", "Education_College_i0",
  "TownsendDeprivationIndex", "SmokingStatus_i0", "Alcohol_Heavy_i0",
  "PAshort_MV_i0", "BMI_i0", "FamilyHistory_HeartDiseaseORStroke_i0",
  "Baseline_AntiHypertensive_i0"
)
outcomes <- tibble::tribble(
  ~outcome, ~baseline, ~event, ~time,
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

# Keep cross-validation folds reproducible for each BPV trait.
make_foldid <- function(n, seed) {
  set.seed(seed)
  sample(rep(seq_len(10L), length.out = n))
}
# Exclude the intercept and coefficients shrunk exactly to zero.
nonzero_coefficients <- function(fit, penalty, features) {
  beta <- as.numeric(coef(fit, s = penalty))[-1L]
  names(beta) <- features
  beta[is.finite(beta) & beta != 0]
}

# Fit the score to BPV residualized on the other BPV trait.
fit_mps <- function(trait, seed) {
  exposure <- paste0(trait, "_VIM_z")
  other <- paste0(if (trait == "SBPV") "DBPV" else "SBPV", "_VIM_z")
  features <- intersect(unique(selected$feature[selected$pressure == trait]),
                        names(cohort))
  d <- cohort[complete.cases(cohort[, c("ID", exposure, other)]), , drop = FALSE]
  x <- as.matrix(d[, features, drop = FALSE])
  storage.mode(x) <- "double"
  usable <- vapply(seq_len(ncol(x)), function(j) {
    assay_sd <- sd(x[, j], na.rm = TRUE)
    sum(!is.na(x[, j])) >= 3L && is.finite(assay_sd) && assay_sd > 0
  }, logical(1))
  x <- x[, usable, drop = FALSE]
  if (anyNA(x)) {
    set.seed(seed)
    x <- imputePCA(x, ncp = min(2L, ncol(x) - 1L), scale = TRUE)$completeObs
  }
  x <- scale(x)
  target <- residuals(lm(reformulate(other, response = exposure), data = d))
  fit <- cv.glmnet(
    x, target, alpha = 0.9, family = "gaussian", nfolds = 10L,
    foldid = make_foldid(nrow(x), seed), intercept = TRUE,
    standardize = FALSE, parallel = FALSE
  )
  beta_1se <- nonzero_coefficients(fit, "lambda.1se", colnames(x))
  beta_min <- nonzero_coefficients(fit, "lambda.min", colnames(x))
  use_1se <- length(beta_1se) > 0L
  beta <- if (use_1se) beta_1se else beta_min
  raw <- as.numeric(x[, names(beta), drop = FALSE] %*% beta)
  if (cor(raw, d[[exposure]], use = "complete.obs") < 0) {
    raw <- -raw
    beta <- -beta
  }
  score <- as.numeric(scale(raw))
  scores <- data.frame(ID = d$ID, value = score)
  names(scores)[2L] <- paste0(trait, "_MPS")
  coefficients <- data.table(
    trait = trait, feature = names(beta), coefficient = unname(beta)
  )[order(-abs(coefficient), feature)]
  list(scores = scores, coefficients = coefficients)
}

sbpv <- fit_mps("SBPV", 20261001L)
dbpv <- fit_mps("DBPV", 20261002L)
scores <- full_join(sbpv$scores, dbpv$scores, by = "ID")
fst::write_fst(scores, file.path(output_dir, "mps_scores.fst"), compress = 50)
# Write the complete weighted-sum formulas.
format_formula <- function(d) {
  terms <- sprintf("%.16g × z(%s)", abs(d$coefficient), d$feature)
  signs <- ifelse(d$coefficient < 0, " - ", " + ")
  signs[1L] <- if (d$coefficient[1L] < 0) "-" else ""
  paste0("z(", paste0(signs, terms, collapse = ""), ")")
}
formulas <- data.frame(
  Trait = c("SBPV", "DBPV"),
  `N of proteins` = c(nrow(sbpv$coefficients), nrow(dbpv$coefficients)),
  Formula = c(format_formula(sbpv$coefficients),
              format_formula(dbpv$coefficients)),
  check.names = FALSE
)
fwrite(formulas, file.path(output_dir, "mps_formulas.csv"))
# Fit outcome Cox models with and without the corresponding measured BPV.
# Convert the cohort's event coding to the numeric indicator used by Surv().
event01 <- function(x) {
  x <- as.character(x)
  ifelse(x %in% c("1", "Yes", "TRUE"), 1L,
         ifelse(x %in% c("0", "No", "FALSE"), 0L, NA_integer_))
}
analysis <- left_join(cohort, scores, by = "ID")

# Estimate the score HR in the outcome-specific, baseline-free risk set.
fit_outcome <- function(trait, o, add_bpv) {
  score_name <- paste0(trait, "_MPS")
  mean_bp <- if (trait == "SBPV") "mean_SBP_0to5y_i0" else "mean_DBP_0to5y_i0"
  adjusters <- c(mean_bp, clinical_covars,
                 if (add_bpv) paste0(trait, "_VIM_z"))
  d <- analysis[as.character(analysis[[o$baseline]]) == "No",
                unique(c(score_name, adjusters, o$event, o$time)), drop = FALSE]
  d$.event <- event01(d[[o$event]])
  d$.time <- suppressWarnings(as.numeric(as.character(d[[o$time]])))
  d <- droplevels(d[complete.cases(d[, c(score_name, adjusters), drop = FALSE]) &
                       !is.na(d$.event) & is.finite(d$.time) & d$.time > 0,
                     , drop = FALSE])
  fit <- coxph(reformulate(c(score_name, adjusters),
                          response = "Surv(.time, .event)"), data = d)
  term <- broom::tidy(fit, exponentiate = TRUE, conf.int = TRUE) |>
    filter(term == score_name)
  data.table(trait = trait, outcome = o$outcome,
             model = if (add_bpv) "Clinical + measured BPV" else "Clinical",
             n = fit$n, events = fit$nevent, hr = term$estimate,
             ci_low = term$conf.low, ci_high = term$conf.high,
             p_value = term$p.value)
}
results <- list()
for (trait in c("SBPV", "DBPV")) {
  outcome_ids <- if (trait == "SBPV") seq_len(5L) else c(3L, 4L)
  for (i in outcome_ids) {
    for (add_bpv in c(FALSE, TRUE)) {
      results[[length(results) + 1L]] <- fit_outcome(trait, outcomes[i, ], add_bpv)
    }
  }
}
fwrite(rbindlist(results), file.path(output_dir, "mps_outcomes.csv"))
