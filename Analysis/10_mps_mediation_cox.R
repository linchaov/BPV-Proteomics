# ==============================================================================
# METADATA
# ==============================================================================
# Author: Chaoyang Lin, MD
#
# Dependencies: R >= 4.2.0; fst, regmedint
#
# Goal: Explore single-MPS mediation with a Cox outcome model.
# Exposure: 1-SD BPV-VIM difference (a1 = 1 versus a0 = 0).
# Output: mediation/mps_cox/mediation_summary.csv
# ==============================================================================

library(fst)
library(regmedint)

root <- "./BPV_ASCVD"
out_dir <- file.path(root, "mediation/mps_cox")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

covars <- c(
  "Age_i0", "Sex", "Ethnic_White_i0", "Education_College_i0",
  "TownsendDeprivationIndex", "SmokingStatus_i0", "Alcohol_Heavy_i0",
  "PAshort_MV_i0", "BMI_i0", "FamilyHistory_HeartDiseaseORStroke_i0",
  "Baseline_AntiHypertensive_i0"
)
outcomes <- c(MI = "MyocardialInfarction", IS = "IschaemicStroke",
              PAD = "PeripheralArterialDisease")
jobs <- data.frame(
  Trait = c("SBPV", "SBPV", "SBPV", "DBPV", "DBPV"),
  Outcome = c("MI", "IS", "PAD", "IS", "PAD")
)

# Join the cohort to the overall MPSs.
outcome_vars <- unlist(lapply(c("Baseline_", "Incident_", "Time_"),
                              function(prefix) paste0(prefix, outcomes, "_i0")))
vars <- unique(c("ID", "SBPV_VIM_z", "DBPV_VIM_z",
                 "mean_SBP_0to5y_i0", "mean_DBP_0to5y_i0",
                 covars, outcome_vars))
cohort <- fst::read_fst(file.path(root, "data/cohort_olink.fst"),
                        columns = vars, as.data.table = FALSE)
scores <- fst::read_fst(file.path(root, "data/overall_mps/mps_scores.fst"),
                        as.data.table = FALSE)
scores <- scores[match(cohort$ID, scores$ID), , drop = FALSE]
data <- cbind(cohort, scores[, -1L, drop = FALSE])

# Restrict to participants free of the outcome at baseline and use complete cases.
prepare_model <- function(data, trait, outcome) {
  exposure <- paste0(trait, "_VIM_z")
  mediator <- paste0(trait, "_MPS")
  adjusters <- c(if (trait == "SBPV") "mean_SBP_0to5y_i0"
                 else "mean_DBP_0to5y_i0", covars)
  suffix <- paste0(outcomes[[outcome]], "_i0")
  baseline <- paste0("Baseline_", suffix)
  event <- paste0("Incident_", suffix)
  time <- paste0("Time_", suffix)
  d <- data[, unique(c(exposure, mediator, adjusters, baseline, event, time)),
            drop = FALSE]
  d <- d[as.character(d[[baseline]]) %in% c("No", "0", "FALSE"), , drop = FALSE]
  e <- as.character(d[[event]])
  d$event <- ifelse(e %in% c("Yes", "1", "TRUE"), 1L,
                    ifelse(e %in% c("No", "0", "FALSE"), 0L, NA_integer_))
  d$years <- as.numeric(as.character(d[[time]])) / 365.25
  d <- droplevels(d[complete.cases(d[, c(exposure, mediator, adjusters,
                                         "event", "years")]) & d$years > 0,
                    , drop = FALSE])
  x <- model.matrix(reformulate(adjusters), data = d)[, -1L, drop = FALSE]
  colnames(x) <- make.names(colnames(x), unique = TRUE)
  model_data <- cbind(d[, c(exposure, mediator, "years", "event"), drop = FALSE],
                      as.data.frame(x))
  list(data = model_data, exposure = exposure, mediator = mediator,
       cvar = colnames(x))
}

# Extract natural direct/indirect and total effects from the Cox model.
# Estimates are hazard-ratio-scale approximations.
fit_model <- function(data, trait, outcome) {
  input <- prepare_model(data, trait, outcome)
  d <- input$data
  fit <- regmedint::regmedint(
    data = d, yvar = "years", eventvar = "event",
    avar = input$exposure, mvar = input$mediator, cvar = input$cvar,
    a0 = 0, a1 = 1, m_cde = median(d[[input$mediator]]),
    c_cond = colMeans(d[, input$cvar, drop = FALSE]),
    mreg = "linear", yreg = "survCox", interaction = FALSE
  )
  s <- summary(fit, exponentiate = TRUE)$summary_myreg
  extract <- function(effect, label) {
    z <- s[effect, ]
    data.frame(Effect = label, Estimate = z["exp(est)"],
               CI_low = z["exp(lower)"], CI_high = z["exp(upper)"])
  }
  effects <- rbind(extract("pnde", "Direct_HR"),
                   extract("tnie", "Indirect_HR"),
                   extract("te", "Total_HR"))
  pm <- s["pm", ]
  effects <- rbind(effects,
                   data.frame(Effect = "PM_percent", Estimate = 100 * pm["est"],
                              CI_low = 100 * pm["lower"],
                              CI_high = 100 * pm["upper"]))
  data.frame(Trait = trait, Outcome = outcome, N = nrow(d),
             Events = sum(d$event), effects, row.names = NULL)
}

results <- lapply(seq_len(nrow(jobs)), function(i) {
  fit_model(data, jobs$Trait[i], jobs$Outcome[i])
})
summary_table <- do.call(rbind, results)
write.csv(summary_table, file.path(out_dir, "mediation_summary.csv"),
          row.names = FALSE, na = "")
