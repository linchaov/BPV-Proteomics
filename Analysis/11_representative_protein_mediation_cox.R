# ==============================================================================
# METADATA
# ==============================================================================
# Author: Chaoyang Lin, MD
# Dependencies: R >= 4.2.0; fst, regmedint
# Goal: Test the Figure 4 functional-group representative proteins individually
#       as potential mediators of BPV associations with MI, ischemic stroke, PAD.
# Input: tables/supp_functional_proteins.csv; data/cohort_olink.fst.
# Output: mediation/functional_representative_proteins/results.csv.
# ==============================================================================

library(fst)
library(regmedint)

root <- "./BPV_ASCVD"
out_dir <- file.path(root, "mediation/functional_representative_proteins")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

covars <- c(
  "Age_i0", "Sex", "Ethnic_White_i0", "Education_College_i0",
  "TownsendDeprivationIndex", "SmokingStatus_i0", "Alcohol_Heavy_i0",
  "PAshort_MV_i0", "BMI_i0", "FamilyHistory_HeartDiseaseORStroke_i0",
  "Baseline_AntiHypertensive_i0"
)
outcomes <- c(MI = "MyocardialInfarction", IS = "IschaemicStroke",
              PAD = "PeripheralArterialDisease")

# Use the protein list, including its F-group mapping.
representatives <- unique(read.csv(
  file.path(root, "data/representative_proteins.csv"),
  check.names = FALSE
)[, c("feature", "belongs_to", "membership")])

outcome_vars <- unlist(lapply(c("Baseline_", "Incident_", "Time_"),
                              function(prefix) paste0(prefix, outcomes, "_i0")))
vars <- unique(c(representatives$feature, "SBPV_VIM_z", "DBPV_VIM_z",
                 "mean_SBP_0to5y_i0", "mean_DBP_0to5y_i0",
                 covars, outcome_vars))
cohort <- fst::read_fst(file.path(root, "data/cohort_olink.fst"),
                        columns = vars, as.data.table = FALSE)

# Use outcome-specific baseline-free cohorts and the covariates.
# Standardize each protein within its model's complete-case sample.
prepare_model <- function(trait, outcome, protein) {
  exposure <- paste0(trait, "_VIM_z")
  adjusters <- c(if (trait == "SBPV") "mean_SBP_0to5y_i0"
                 else "mean_DBP_0to5y_i0", covars)
  suffix <- paste0(outcomes[[outcome]], "_i0")
  baseline <- paste0("Baseline_", suffix)
  event <- paste0("Incident_", suffix)
  time <- paste0("Time_", suffix)
  d <- cohort[, unique(c(exposure, protein, adjusters, baseline, event, time)),
              drop = FALSE]
  d <- d[as.character(d[[baseline]]) %in% c("No", "0", "FALSE"), , drop = FALSE]
  e <- as.character(d[[event]])
  d$event <- ifelse(e %in% c("Yes", "1", "TRUE"), 1L,
                    ifelse(e %in% c("No", "0", "FALSE"), 0L, NA_integer_))
  d$years <- as.numeric(as.character(d[[time]])) / 365.25
  d <- droplevels(d[complete.cases(d[, c(exposure, protein, adjusters,
                                         "event", "years")]) & d$years > 0,
                    , drop = FALSE])
  d$protein_z <- as.numeric(scale(d[[protein]]))
  x <- model.matrix(reformulate(adjusters), data = d)[, -1L, drop = FALSE]
  colnames(x) <- make.names(colnames(x), unique = TRUE)
  model_data <- cbind(d[, c(exposure, "protein_z", "years", "event"),
                        drop = FALSE], as.data.frame(x))
  list(data = model_data, exposure = exposure, cvar = colnames(x))
}

# The PM test and confidence interval come from regmedint's Cox mediation model.
# Retain direct, indirect and total HR-scale estimates for interpretation.
fit_model <- function(trait, outcome, protein) {
  input <- prepare_model(trait, outcome, protein)
  d <- input$data
  fit <- regmedint::regmedint(
    data = d, yvar = "years", eventvar = "event", avar = input$exposure,
    mvar = "protein_z", cvar = input$cvar,
    a0 = 0, a1 = 1, m_cde = 0,
    c_cond = colMeans(d[, input$cvar, drop = FALSE]),
    mreg = "linear", yreg = "survCox", interaction = FALSE
  )
  s <- summary(fit, exponentiate = TRUE)$summary_myreg
  # Read the named regmedint summary row and column for each effect.
  hr <- function(effect, column) unname(s[effect, column])
  data.frame(
    Trait = trait, Outcome = outcome, Protein = protein,
    N = nrow(d), Events = sum(d$event),
    Indirect_HR = hr("tnie", "exp(est)"),
    Indirect_CI_low = hr("tnie", "exp(lower)"),
    Indirect_CI_high = hr("tnie", "exp(upper)"),
    Indirect_p = hr("tnie", "p"),
    PM_percent = 100 * hr("pm", "est"),
    PM_CI_low = 100 * hr("pm", "lower"),
    PM_CI_high = 100 * hr("pm", "upper"),
    PM_p = hr("pm", "p")
  )
}

jobs <- expand.grid(
  Trait = c("SBPV", "DBPV"), Outcome = names(outcomes),
  Protein = representatives$feature, stringsAsFactors = FALSE
)
jobs <- merge(jobs, representatives, by.x = "Protein", by.y = "feature",
              sort = FALSE)
jobs <- jobs[(jobs$Trait == "SBPV" | jobs$Outcome != "MI") &
             (jobs$Trait == "SBPV" | jobs$membership != "SBPV only"),
             , drop = FALSE]
jobs <- jobs[order(match(jobs$Trait, c("SBPV", "DBPV")),
                   match(jobs$Outcome, names(outcomes)),
                   match(jobs$Protein, representatives$feature)), ]

results <- lapply(seq_len(nrow(jobs)), function(i) {
  job <- jobs[i, ]
  cbind(belongs_to = job$belongs_to,
        fit_model(job$Trait, job$Outcome, job$Protein))
})
results <- do.call(rbind, results)
rownames(results) <- NULL
write.csv(results, file.path(out_dir, "results.csv"), row.names = FALSE)
