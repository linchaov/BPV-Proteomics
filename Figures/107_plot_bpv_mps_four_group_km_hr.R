# ==============================================================================
# METADATA
# ==============================================================================
# Author: Chaoyang Lin, MD
# Goal: Draw four-group KM curves and adjusted HR forest plots for measured
#       SBPV/DBPV and their overall multi-protein scores (MPSs).
# Input: Analysis cohort and MPS scores from scripts 01 and 08.
# Output: Four PDFs in ./BPV_ASCVD/figures/.
# ==============================================================================
suppressPackageStartupMessages({
  library(fst)
  library(survival)
  library(ggplot2)
  library(patchwork)
})

data_dir <- "./BPV_ASCVD/data"
figure_dir <- "./BPV_ASCVD/figures"

dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

clinical_covars <- c(
  "Age_i0", "Sex", "Ethnic_White_i0", "Education_College_i0",
  "TownsendDeprivationIndex", "SmokingStatus_i0", "Alcohol_Heavy_i0",
  "PAshort_MV_i0", "BMI_i0", "FamilyHistory_HeartDiseaseORStroke_i0",
  "Baseline_AntiHypertensive_i0", "mean_SBP_0to5y_i0"
)

outcomes <- data.frame(
  label = c("MI", "Coronary revascularization", "Ischemic stroke",
            "PAD", "MALE", "All-cause mortality"),
  suffix = c("MyocardialInfarction", "CoronaryRevascularization",
             "IschaemicStroke", "PeripheralArterialDisease",
             "MajorAdverseLimbEvent", NA_character_)
)

outcome_columns <- unlist(lapply(c("Baseline_", "Incident_", "Time_"),
  function(prefix) paste0(prefix, outcomes$suffix[1:5], "_i0")))

cohort <- fst::read_fst(file.path(data_dir, "cohort_olink.fst"),
  columns = c("ID", "SBPV_VIM_z", "DBPV_VIM_z", clinical_covars,
              outcome_columns, "Death_Any_AllCause_i0",
              "Time_Death_Summary_i0"), as.data.table = FALSE)

groups <- c("Low/Low", "High/Low", "Low/High", "High/High")
comparison_groups <- groups[-1L]

colours <- c("Low/Low" = "#4D4D4D", "High/Low" = "#B36A4C",
             "Low/High" = "#4C78A8", "High/High" = "#8A4F7D")

# Four evenly spaced positions per outcome: LL, HL, LH and HH from top to bottom.
offsets <- c("Low/Low" = 0.30, "High/Low" = 0.10,
             "Low/High" = -0.10, "High/High" = -0.30)

limit_days <- 16 * 365.25

# Split two continuous traits at their medians among participants with both.
make_groups <- function(d, first, second) {
  eligible <- complete.cases(d[, c(first, second)])
  a <- median(d[[first]][eligible])
  b <- median(d[[second]][eligible])
  label <- ifelse(d[[first]] <= a,
    ifelse(d[[second]] <= b, "Low/Low", "Low/High"),
    ifelse(d[[second]] <= b, "High/Low", "High/High"))
  label[!eligible] <- NA_character_
  factor(label, levels = groups)
}

# Convert the cohort's Yes/No or 1/0 outcome fields to survival event status.
event01 <- function(x) {
  value <- as.character(x)
  ifelse(value %in% c("Yes", "1", "TRUE"), 1L,
         ifelse(value %in% c("No", "0", "FALSE"), 0L, NA_integer_))
}

# Prepare the original outcome-specific risk set; death uses the full cohort.
outcome_data <- function(d, i) {
  if (i == 6L) {
    event <- "Death_Any_AllCause_i0"
    time <- "Time_Death_Summary_i0"
    eligible <- rep(TRUE, nrow(d))
  } else {
    stem <- outcomes$suffix[i]
    event <- paste0("Incident_", stem, "_i0")
    time <- paste0("Time_", stem, "_i0")
    eligible <- as.character(d[[paste0("Baseline_", stem, "_i0")]]) == "No"
  }
  out <- d[which(eligible), unique(c("group", clinical_covars, event, time)),
           drop = FALSE]
  # Convert the cohort's Yes/No or 1/0 outcome fields to survival event status.
  out$.event <- event01(out[[event]])
  out$.time <- as.numeric(as.character(out[[time]]))
  out <- out[!is.na(out$group) & !is.na(out$.event) &
             is.finite(out$.time) & out$.time > 0, , drop = FALSE]
  out$group <- droplevels(out$group)
  out
}

# Measured SBPV + DBPV: median-defined four groups.
local({

# KM curves: measured_bpv_four_group_km.pdf
local({

data <- cohort

# Draw one survival panel with the original log-rank test and shared group colours.
plot_panel <- function(d, i, legend_labels) {
  fit <- survival::survfit(Surv(.time, .event) ~ group, data = d)
  test <- survival::survdiff(Surv(.time, .event) ~ group, data = d)
  p <- pchisq(test$chisq, df = length(test$n) - 1L, lower.tail = FALSE)
  p_text <- if (p < 0.001) "Log-rank P < 0.001" else
    sprintf("Log-rank P = %.3f", p)
  annotation <- sprintf("%s\nN = %s; events = %s", p_text,
                        format(nrow(d), big.mark = ",", scientific = FALSE),
                        format(sum(d$.event), big.mark = ",", scientific = FALSE))

  steps <- summary(fit, censored = FALSE)
  curve <- data.frame(time = steps$time, survival = steps$surv,
                      group = sub("^group=", "", steps$strata))
  end <- summary(fit, times = limit_days, extend = TRUE)
  curve <- rbind(data.frame(time = 0, survival = 1, group = groups), curve,
                 data.frame(time = end$time, survival = end$surv,
                            group = sub("^group=", "", end$strata)))
  curve <- curve[curve$time <= limit_days, , drop = FALSE]
  curve$group <- factor(curve$group, levels = groups)
  curve <- curve[order(curve$group, curve$time), , drop = FALSE]
  y_min <- max(0, floor((min(curve$survival) - 0.02) * 20) / 20)

  ggplot(curve, aes(time / 365.25, survival, colour = group)) +
    geom_step(linewidth = 0.65, direction = "hv") +
    annotate("text", x = 0.45, y = y_min + 0.11 * (1 - y_min),
             label = annotation, hjust = 0, size = 4.5, lineheight = 1.15,
             family = "sans") +
    scale_colour_manual(values = colours, limits = groups,
                        labels = legend_labels) +
    scale_x_continuous(breaks = seq(0, 16, 4), limits = c(0, 16), expand = c(0, 0)) +
    scale_y_continuous(breaks = pretty(c(y_min, 1), n = 4),
                       limits = c(y_min, 1), expand = c(0, 0)) +
    labs(title = outcomes$label[i],
         x = "Time (years)", y = "Survival probability", colour = NULL) +
    guides(colour = guide_legend(nrow = 2, byrow = TRUE,
                                  override.aes = list(linewidth = 1.3))) +
    theme_classic(base_family = "sans", base_size = 13) +
    theme(plot.title = element_text(face = "bold", size = 14),
          axis.title = element_text(size = 13),
          axis.text = element_text(size = 12, colour = "black"),
          panel.grid.major.y = element_line(colour = "#ECECEC", linewidth = 0.25),
          legend.position = "bottom",
          legend.key.width = grid::unit(20, "mm"),
          legend.text = element_text(size = 12),
          plot.margin = margin(8, 8, 5, 8))
}

# Assemble six outcomes in reading order and export one page for each grouping.
save_six_panels <- function(first, second, file_stem, legend_labels) {
  data$group <- make_groups(data, first, second)
  panels <- lapply(seq_len(nrow(outcomes)), function(i) {
    plot_panel(outcome_data(data, i), i, legend_labels)
  })
  figure <- wrap_plots(panels, ncol = 3, guides = "collect") &
    theme(legend.position = "bottom")
  ggsave(file.path(figure_dir, paste0(file_stem, ".pdf")), figure,
         width = 12.5, height = 9.0, device = cairo_pdf)
}

measured_labels <- c(
  "Low SBPV & Low DBPV", "High SBPV & Low DBPV",
  "Low SBPV & High DBPV", "High SBPV & High DBPV"
)
save_six_panels("SBPV_VIM_z", "DBPV_VIM_z",
                "measured_bpv_four_group_km", measured_labels)

})

# Adjusted HR forest: measured_bpv_four_group_hr_forest.pdf
local({

data <- cohort

# Keep baseline-free participants and complete covariates for each outcome.
# Translate Yes/No and 1/0 fields into model-ready outcome indicators.
prepare_outcome <- function(d, i) {
  if (i == 6L) {
    event <- "Death_Any_AllCause_i0"
    time <- "Time_Death_Summary_i0"
    eligible <- rep(TRUE, nrow(d))
  } else {
    stem <- outcomes$suffix[i]
    event <- paste0("Incident_", stem, "_i0")
    time <- paste0("Time_", stem, "_i0")
    eligible <- as.character(d[[paste0("Baseline_", stem, "_i0")]]) == "No"
  }
  keep <- c("group", clinical_covars, event, time)
  out <- d[which(eligible), unique(keep), drop = FALSE]
  out$.event <- event01(out[[event]])
  out$.time <- as.numeric(as.character(out[[time]]))
  required <- c("group", clinical_covars, ".event", ".time")
  out <- out[complete.cases(out[, required, drop = FALSE]) &
               out$.time > 0, , drop = FALSE]
  out$group <- droplevels(out$group)
  out
}

# Fit one adjusted Cox model after omitting groups with fewer than five events.
fit_outcome <- function(d, i) {
  event_counts <- table(factor(d$group[d$.event == 1L], levels = groups))
  eligible_groups <- names(event_counts)[event_counts >= 5L]
  if (!"Low/Low" %in% eligible_groups || length(eligible_groups) < 2L) {
    return(data.frame(outcome = outcomes$label[i], group = comparison_groups,
      estimate = NA_real_, low = NA_real_, high = NA_real_,
      y = nrow(outcomes) + 1L - i + offsets[comparison_groups]))
  }
  d <- droplevels(d[d$group %in% eligible_groups, , drop = FALSE])
  fit <- coxph(reformulate(c("group", clinical_covars),
                          response = "Surv(.time, .event)"), data = d)
  terms <- paste0("group", comparison_groups)
  beta <- coef(fit)[terms]
  se <- sqrt(diag(vcov(fit))[terms])
  result <- data.frame(outcome = outcomes$label[i], group = comparison_groups,
                       estimate = exp(beta), low = exp(beta - 1.96 * se),
                       high = exp(beta + 1.96 * se),
                       y = nrow(outcomes) + 1L - i + offsets[comparison_groups])
  result
}

# Draw the LL reference square (HR = 1, no CI) and three estimated contrasts.
plot_forest <- function(result) {
  result <- result[is.finite(result$estimate) & is.finite(result$low) &
                     is.finite(result$high), , drop = FALSE]
  result$group <- factor(result$group, levels = comparison_groups)
  x_range <- range(c(result$low, result$high), finite = TRUE)
  x_limits <- c(min(0.6, x_range[1] * 0.92), max(2.6, x_range[2] * 1.08))
  breaks <- c(0.5, 1, 2, 4, 8)
  breaks <- breaks[breaks >= x_limits[1] & breaks <= x_limits[2]]
  reference <- data.frame(
    estimate = 1, y = rev(seq_len(nrow(outcomes))) + offsets["Low/Low"]
  )

  ggplot(result, aes(y = y, colour = group)) +
    geom_hline(yintercept = seq(1.5, 5.5, 1), colour = "#EAEAEA",
               linewidth = 0.35) +
    geom_vline(xintercept = 1, colour = "#777777", linetype = "dashed",
               linewidth = 0.45) +
    geom_segment(aes(x = low, xend = high, yend = y), linewidth = 0.8) +
    geom_point(aes(x = estimate), size = 3.0) +
    geom_point(data = reference, aes(x = estimate, y = y),
               inherit.aes = FALSE, shape = 15, size = 2.5,
               colour = "black", show.legend = FALSE) +
    scale_colour_manual(values = colours) +
    scale_x_log10(breaks = breaks, labels = as.character(breaks),
                  limits = x_limits) +
    scale_y_continuous(breaks = 6:1,
                       labels = sub("Coronary revascularization",
                                    "Coronary\nrevascularization", outcomes$label),
                       limits = c(0.5, 6.5), expand = c(0, 0)) +
    labs(x = "Adjusted HR (log scale)", y = NULL) +
    theme_classic(base_family = "sans", base_size = 13) +
    theme(axis.text.y = element_text(size = 13, face = "bold",
                                       colour = "black"),
          axis.text.x = element_text(size = 13, colour = "black"),
          axis.title.x = element_text(size = 14),
          panel.grid.major.x = element_line(colour = "#EEEEEE",
                                             linewidth = 0.3),
          legend.position = "none",
          plot.margin = margin(8, 10, 8, 8))
}

# Export the adjusted HR forest for the measured SBPV/DBPV groups.
# Apply the same two-median, four-group definition as the KM figure.
save_forest <- function(first, second, file_stem) {
  data$group <- make_groups(data, first, second)
  result <- do.call(rbind, lapply(seq_len(nrow(outcomes)), function(i) {
    fit_outcome(prepare_outcome(data, i), i)
  }))
  figure <- plot_forest(result)
  ggsave(file.path(figure_dir, paste0(file_stem, ".pdf")), figure,
         width = 7.4, height = 7.1, device = cairo_pdf)
}

save_forest("SBPV_VIM_z", "DBPV_VIM_z",
  "measured_bpv_four_group_hr_forest")

})

})

# MPS-SBPV + MPS-DBPV: median-defined four groups.
local({

# Draw one survival panel with the original log-rank test and shared group colours.
plot_panel <- function(d, i, legend_labels) {
  fit <- survival::survfit(Surv(.time, .event) ~ group, data = d)
  test <- survival::survdiff(Surv(.time, .event) ~ group, data = d)
  p <- pchisq(test$chisq, df = length(test$n) - 1L, lower.tail = FALSE)
  p_text <- if (p < 0.001) "Log-rank P < 0.001" else
    sprintf("Log-rank P = %.3f", p)
  annotation <- sprintf("%s\nN = %s\nEvents = %s", p_text,
                        format(nrow(d), big.mark = ",", scientific = FALSE),
                        format(sum(d$.event), big.mark = ",", scientific = FALSE))

  steps <- summary(fit, censored = FALSE)
  curve <- data.frame(time = steps$time, survival = steps$surv,
                      group = sub("^group=", "", steps$strata))
  end <- summary(fit, times = limit_days, extend = TRUE)
  curve <- rbind(data.frame(time = 0, survival = 1, group = groups), curve,
                 data.frame(time = end$time, survival = end$surv,
                            group = sub("^group=", "", end$strata)))
  curve <- curve[curve$time <= limit_days, , drop = FALSE]
  curve$group <- factor(curve$group, levels = groups)
  curve <- curve[order(curve$group, curve$time), , drop = FALSE]
  y_min <- max(0, floor((min(curve$survival) - 0.02) * 20) / 20)

  ggplot(curve, aes(time / 365.25, survival, colour = group)) +
    geom_step(linewidth = 0.85, direction = "hv") +
    annotate("text", x = 0.45, y = y_min + 0.02 * (1 - y_min),
             label = annotation, hjust = 0, vjust = 0, size = 5.3,
             lineheight = 1.05,
             family = "sans") +
    scale_colour_manual(values = colours, limits = groups,
                        labels = legend_labels) +
    scale_x_continuous(breaks = seq(0, 16, 4), limits = c(0, 16), expand = c(0, 0)) +
    scale_y_continuous(breaks = pretty(c(y_min, 1), n = 4),
                       limits = c(y_min, 1), expand = c(0, 0)) +
    labs(title = outcomes$label[i],
         x = "Time (years)", y = "Survival probability", colour = NULL) +
    guides(colour = guide_legend(nrow = 2, byrow = TRUE,
                                  override.aes = list(linewidth = 1.3))) +
    theme_classic(base_family = "sans", base_size = 16) +
    theme(plot.title = element_text(face = "bold", size = 18),
          axis.title = element_text(size = 16),
          axis.text = element_text(size = 15, colour = "black"),
          panel.grid.major.y = element_line(colour = "#ECECEC", linewidth = 0.25),
          legend.position = "bottom",
          legend.key.width = grid::unit(18, "mm"),
          legend.text = element_text(size = 14),
          plot.margin = margin(8, 8, 5, 8))
}

# KM curves: mps_sbpv_dbpv_four_group_km.pdf
# Prepare the original outcome-specific risk set; death uses the full cohort.
save_km <- function(data, file_stem) {
  panels <- lapply(seq_len(nrow(outcomes)), function(i) {
    plot_panel(outcome_data(data, i), i, mps_labels)
  })
  figure <- wrap_plots(panels, ncol = 3, guides = "collect") &
    theme(legend.position = "bottom")
  ggsave(file.path(figure_dir, paste0(file_stem, "_km.pdf")), figure,
         width = 10.5, height = 8.4, device = cairo_pdf)
}

mps_labels <- c(
  "Low MPS-SBPV & Low MPS-DBPV", "High MPS-SBPV & Low MPS-DBPV",
  "Low MPS-SBPV & High MPS-DBPV", "High MPS-SBPV & High MPS-DBPV"
)
# Four evenly spaced positions per outcome: LL, HL, LH and HH from top to bottom.
# Fit one outcome-specific Cox model; exclude groups with fewer than five events.
fit_hr <- function(d, i) {
  d <- d[complete.cases(d[, clinical_covars, drop = FALSE]), , drop = FALSE]
  event_counts <- table(factor(d$group[d$.event == 1L], levels = groups))
  eligible_groups <- names(event_counts)[event_counts >= 5L]
  y <- nrow(outcomes) + 1L - i + offsets[comparison_groups]
  if (!"Low/Low" %in% eligible_groups || length(eligible_groups) < 2L) {
    return(data.frame(group = comparison_groups, estimate = NA_real_,
                      low = NA_real_, high = NA_real_, y = y))
  }
  d <- droplevels(d[d$group %in% eligible_groups, , drop = FALSE])
  fit <- coxph(reformulate(c("group", clinical_covars),
                          response = "Surv(.time, .event)"), data = d)
  terms <- paste0("group", comparison_groups)
  beta <- unname(coef(fit)[terms])
  se <- unname(sqrt(diag(vcov(fit))[terms]))
  data.frame(group = comparison_groups, estimate = exp(beta),
             low = exp(beta - 1.96 * se), high = exp(beta + 1.96 * se), y = y)
}

# Display the LL reference square and three estimated contrasts per outcome.
# LL is fixed at HR = 1 by definition and has no estimated confidence interval.
plot_hr <- function(result) {
  result <- result[is.finite(result$estimate) & is.finite(result$low) &
                     is.finite(result$high), , drop = FALSE]
  result$group <- factor(result$group, levels = comparison_groups)
  x_range <- range(c(result$low, result$high), finite = TRUE)
  x_limits <- c(min(0.6, x_range[1] * 0.92), max(2.6, x_range[2] * 1.08))
  # Use sparse, equally spaced doubling ticks on the logarithmic HR axis.
  breaks <- c(0.5, 1, 2, 4, 8)
  breaks <- breaks[breaks >= x_limits[1] & breaks <= x_limits[2]]
  reference <- data.frame(
    estimate = 1, y = rev(seq_len(nrow(outcomes))) + offsets["Low/Low"]
  )

  ggplot(result, aes(y = y, colour = group)) +
    geom_hline(yintercept = seq(1.5, 5.5, 1), colour = "#EAEAEA",
               linewidth = 0.35) +
    geom_vline(xintercept = 1, colour = "#777777", linetype = "dashed",
               linewidth = 0.45) +
    geom_segment(aes(x = low, xend = high, yend = y), linewidth = 1.1) +
    geom_point(aes(x = estimate), size = 4.0) +
    geom_point(data = reference, aes(x = estimate, y = y),
               inherit.aes = FALSE, shape = 15, size = 3.0,
               colour = "black", show.legend = FALSE) +
    scale_colour_manual(values = colours[comparison_groups]) +
    scale_x_log10(breaks = breaks, labels = as.character(breaks),
                  limits = x_limits) +
    scale_y_continuous(breaks = 6:1,
                       labels = sub("Coronary revascularization",
                                    "Coronary\nrevascularization", outcomes$label),
                       limits = c(0.5, 6.5), expand = c(0, 0)) +
    labs(x = "Adjusted HR (log scale)", y = NULL) +
    theme_classic(base_family = "sans", base_size = 16) +
    theme(axis.text.y = element_text(size = 18, face = "bold",
                                     colour = "black"),
          axis.text.x = element_text(size = 16, colour = "black"),
          axis.title.x = element_text(size = 18),
          panel.grid.major.x = element_line(colour = "#EEEEEE",
                                           linewidth = 0.3),
          legend.position = "none", plot.margin = margin(8, 10, 8, 8))
}

# Adjusted HR forest: mps_sbpv_dbpv_four_group_hr_forest.pdf
save_hr <- function(data, file_stem, figure_width) {
  result <- do.call(rbind, lapply(seq_len(nrow(outcomes)), function(i) {
    fit_hr(outcome_data(data, i), i)
  }))
  figure <- plot_hr(result)
  ggsave(file.path(figure_dir, paste0(file_stem, "_hr_forest.pdf")), figure,
         width = figure_width, height = 9.0, device = cairo_pdf)
}

# Use the saved elastic-net scores from script 08.
scores <- fst::read_fst(
  file.path(data_dir, "overall_mps", "mps_scores.fst"),
  columns = c("ID", "SBPV_MPS", "DBPV_MPS"), as.data.table = FALSE
)
scores <- scores[match(cohort$ID, scores$ID), , drop = FALSE]
data <- cbind(cohort, scores[, -1L, drop = FALSE])
# Split both MPSs at medians calculated among participants with both scores.
data$group <- make_groups(data, "SBPV_MPS", "DBPV_MPS")
save_km(data, "mps_sbpv_dbpv_four_group")
save_hr(data, "mps_sbpv_dbpv_four_group", 8.3)
})
