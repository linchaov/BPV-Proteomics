# ==============================================================================
# METADATA
# ==============================================================================
# Author: Chaoyang Lin, MD
# Goal: Draw cross-outcome protein heatmap, UpSet plot, and direction counts.
# Input: Saved analysis outputs in ./BPV_ASCVD.
# Outputs: PDFs and figure tables under ./BPV_ASCVD.
# ==============================================================================

# --- cross_outcome_protein_heatmap ---
local({

library(data.table)
library(dplyr)
library(ggplot2)

table_dir <- "./BPV_ASCVD/tables"
figure_dir <- "./BPV_ASCVD/figures"

dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

outcomes <- data.table(
  outcome = c(
    "Myocardial infarction", "Coronary revascularization",
    "Ischaemic stroke", "Peripheral arterial disease",
    "Major adverse limb event"
  ),
  display = c(
    "MI", "Coronary\nrevascularization", "Ischemic\nstroke", "PAD", "MALE"
  )
)

all_results <- fread(file.path(table_dir, "protein_outcome_cox_all.csv")) |>
  mutate(
    outcome = factor(outcome, levels = outcomes$outcome),
    display = factor(outcomes$display[match(outcome, outcomes$outcome)],
                     levels = outcomes$display),
    significant = !is.na(fdr) & fdr < 0.05
  )

protein_summary <- as.data.table(all_results)[
  , .(
    positive_significant_outcomes = sum(significant & log_hr > 0),
    negative_significant_outcomes = sum(significant & log_hr < 0),
    mean_log_hr_four = mean(log_hr[
      as.character(outcome) != "Major adverse limb event"
    ])
  ),
  by = protein
]

positive <- protein_summary[
  positive_significant_outcomes >= 2L &
    negative_significant_outcomes == 0L
]
setorder(positive, -mean_log_hr_four, -positive_significant_outcomes, protein)
positive <- head(positive, 20L)

negative <- protein_summary[
  negative_significant_outcomes >= 2L &
    positive_significant_outcomes == 0L
]
setorder(negative, mean_log_hr_four, -negative_significant_outcomes, protein)
negative <- head(negative, 12L)

protein_levels <- c(positive$protein, negative$protein)

heatmap_data <- as.data.table(all_results)[
  protein %in% protein_levels,
  .(protein, display = as.character(display), hr, log_hr, fdr, significant)
]
heatmap_data[, plotted_log_hr := fifelse(significant, log_hr, NA_real_)]
heatmap_data[, display := factor(display, levels = rev(outcomes$display))]
heatmap_data[, protein := factor(protein, levels = protein_levels)]

# The fill is log(HR), so HR = 1 sits at the neutral midpoint. Label the
# shared colour bar on the original HR scale; grey denotes FDR >= 0.05.
fill_limits <- range(heatmap_data$plotted_log_hr, na.rm = TRUE)
hr_breaks <- c(0.6, 0.8, 1, 1.5, 2, 2.5)
hr_breaks <- hr_breaks[log(hr_breaks) >= fill_limits[1] &
                         log(hr_breaks) <= fill_limits[2]]

combined_plot <- ggplot(
  heatmap_data, aes(protein, display, fill = plotted_log_hr)
) +
  geom_tile(colour = "white", linewidth = 0.8) +
  geom_vline(xintercept = nrow(positive) + 0.5, colour = "white", linewidth = 1.6) +
  scale_fill_gradientn(
    colours = c("#08306B", "#F7F7F7", "#8B1A1A"),
    values = scales::rescale(c(fill_limits[1], 0, fill_limits[2]),
                             from = fill_limits),
    limits = fill_limits,
    breaks = log(hr_breaks), labels = format(hr_breaks, trim = TRUE),
    na.value = "#ECECEC", name = "Hazard ratio",
    guide = guide_colourbar(
      title.position = "top", title.hjust = 0.5,
      barheight = grid::unit(55, "mm"),
      barwidth = grid::unit(4, "mm")
    )
  ) +
  scale_x_discrete(drop = FALSE) +
  scale_y_discrete(drop = FALSE) +
  coord_fixed() +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 14, base_family = "sans") +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(
      size = 13, face = "bold", colour = "#3A3A3A",
      angle = 55, hjust = 1, vjust = 1
    ),
    axis.text.y = element_text(size = 15, face = "bold", colour = "#4A4A4A"),
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 12),
    legend.position = "right",
    plot.margin = margin(6, 12, 8, 8)
  )

ggsave(
  file.path(figure_dir, "cross_outcome_consistent_proteins.pdf"),
  combined_plot,
  width = 15,
  height = 5.2,
  device = cairo_pdf,
  bg = "white"
)
})

# --- cross_outcome_protein_upset ---
local({

library(readxl)
library(data.table)
library(ggplot2)
library(patchwork)

input_file <- "./BPV_ASCVD/tables/supp_protein_outcomes.xlsx"
figure_file <- "./BPV_ASCVD/figures/cross_outcome_protein_upset.pdf"

dir.create(dirname(figure_file), recursive = TRUE, showWarnings = FALSE)

# Preserve the manuscript's outcome order while using compact display labels.
outcomes <- data.table(
  sheet = c(
    "Myocardial infarction",
    "Coronary revascularization",
    "Ischaemic stroke",
    "Peripheral arterial disease",
    "Major adverse limb event"
  ),
  set = c("MI", "Coronary revascularization", "Ischemic stroke", "PAD", "MALE")
)

protein_sets <- setNames(lapply(outcomes$sheet, function(sheet_name) {
  x <- readxl::read_excel(input_file, sheet = sheet_name)
  sort(unique(x[["Protein name"]][!is.na(x[["Protein name"]])]))
}), outcomes$set)

membership <- data.table(Protein = sort(unique(unlist(protein_sets, use.names = FALSE))))
for (set_name in outcomes$set) {
  membership[, (set_name) := Protein %in% protein_sets[[set_name]]]
}

membership[, degree := rowSums(.SD), .SDcols = outcomes$set]
membership[, intersection := apply(.SD, 1L, function(z) {
  paste(outcomes$set[as.logical(z)], collapse = " + ")
}), .SDcols = outcomes$set]

# Count mutually exclusive intersections. A protein contributes to exactly one bar.
intersections <- membership[
  , .(intersection_size = .N),
  by = c(outcomes$set, "degree", "intersection")
]
setorder(intersections, -intersection_size, -degree, intersection)
intersections[, intersection_id := seq_len(.N)]

set_sizes <- data.table(
  set = outcomes$set,
  set_size = vapply(protein_sets, length, integer(1))
)

matrix_data <- melt(
  intersections,
  id.vars = c("intersection_id", "intersection_size", "degree", "intersection"),
  measure.vars = outcomes$set,
  variable.name = "set", value.name = "present"
)[, .(intersection_id, set = as.character(set), present)]
setorder(matrix_data, intersection_id, set)

set_levels <- rev(outcomes$set)
matrix_data[, set := factor(set, levels = set_levels)]
set_sizes[, set := factor(set, levels = set_levels)]
intersections[, intersection_id_factor := factor(
  intersection_id, levels = intersection_id
)]
matrix_data[, intersection_id_factor := factor(
  intersection_id, levels = intersections$intersection_id
)]

# Connect the highest and lowest included outcomes within each intersection.
segments <- matrix_data[present == TRUE,
  .(
    y_min = set_levels[min(as.numeric(set))],
    y_max = set_levels[max(as.numeric(set))]
  ),
  by = intersection_id_factor
][y_min != y_max]
segments[, y_min := factor(y_min, levels = set_levels)]
segments[, y_max := factor(y_max, levels = set_levels)]

intersection_plot <- ggplot(
  intersections,
  aes(intersection_id_factor, intersection_size)
) +
  geom_col(width = 0.72, fill = "#5F6B73", colour = NA) +
  geom_text(
    aes(label = intersection_size),
    vjust = -0.35, size = 4.2, family = "sans"
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.12)),
    breaks = scales::breaks_pretty(n = 5)
  ) +
  labs(x = NULL, y = "Intersection size") +
  theme_classic(base_size = 14, base_family = "sans") +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),
    axis.title.y = element_text(size = 15),
    axis.text.y = element_text(size = 13, colour = "black"),
    plot.margin = margin(5, 8, 0, 0)
  )

matrix_plot <- ggplot(matrix_data, aes(intersection_id_factor, set)) +
  geom_segment(
    data = segments,
    aes(
      x = intersection_id_factor,
      xend = intersection_id_factor,
      y = y_min,
      yend = y_max
    ),
    inherit.aes = FALSE,
    colour = "#333333",
    linewidth = 0.55
  ) +
  geom_point(colour = "#D9D9D9", size = 2.6) +
  geom_point(
    data = matrix_data[present == TRUE],
    colour = "#333333", size = 3.0
  ) +
  scale_x_discrete(expand = expansion(add = 0.55)) +
  scale_y_discrete(drop = FALSE) +
  labs(x = NULL, y = NULL) +
  theme_void(base_size = 14, base_family = "sans") +
  theme(
    axis.text.y = element_blank(),
    plot.margin = margin(0, 8, 5, 0)
  )

set_size_plot <- ggplot(set_sizes, aes(set_size, set)) +
  geom_col(width = 0.62, fill = "#6F6F6F") +
  geom_text(
    aes(x = set_size / 2, label = set_size),
    hjust = 0.5, colour = "white", size = 4.3, family = "sans"
  ) +
  scale_x_reverse(
    expand = expansion(mult = c(0.08, 0)),
    breaks = scales::breaks_pretty(n = 3)
  ) +
  scale_y_discrete(drop = FALSE) +
  labs(x = "Set size", y = NULL) +
  theme_classic(base_size = 14, base_family = "sans") +
  theme(
    axis.title.x = element_text(size = 15),
    axis.text.x = element_text(size = 13, colour = "black"),
    axis.text.y = element_text(size = 15, colour = "black"),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank(),
    plot.margin = margin(0, 6, 5, 8)
  )

upset_plot <- wrap_plots(
  A = plot_spacer(),
  B = intersection_plot,
  C = set_size_plot,
  D = matrix_plot,
  design = "AABBBBBBBB\nAABBBBBBBB\nCCDDDDDDDD",
  heights = c(0.46, 0.46, 0.36)
)

ggsave(
  figure_file,
  upset_plot,
  width = 13.2,
  height = 7.6,
  device = cairo_pdf,
  bg = "white"
)
})

# --- bpv_protein_direction_counts ---
local({

library(data.table)
library(ggplot2)

root <- "./BPV_ASCVD"
table_dir <- file.path(root, "tables")
figure_dir <- file.path(root, "figures")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

bpv <- fread(file.path(table_dir, "supp_bpv_proteins.csv"))[
  Exposure %in% c("SBPV-VIM Z-Score", "DBPV-VIM Z-Score"), unique(`Protein name`)
]

associations <- fread(file.path(table_dir, "protein_outcome_cox_all.csv"))[
  protein %in% bpv & !is.na(fdr) & fdr < 0.05 & is.finite(log_hr)
]
associations[, Direction := ifelse(log_hr > 0, "Positive", "Negative")]

outcomes <- data.table(
  outcome = c("Myocardial infarction", "Coronary revascularization",
              "Ischaemic stroke", "Peripheral arterial disease",
              "Major adverse limb event"),
  label = c("MI", "Coronary\nrevascularization", "Ischemic\nstroke",
            "PAD", "MALE"),
  order = 1:5
)

counts <- associations[, .(Proteins = uniqueN(protein)), by = .(outcome, Direction)]
counts <- merge(
  CJ(outcome = outcomes$outcome, Direction = c("Positive", "Negative"), unique = TRUE),
  counts, by = c("outcome", "Direction"), all.x = TRUE
)
counts[is.na(Proteins), Proteins := 0L]
counts <- merge(counts, outcomes, by = "outcome")
setorder(counts, order, Direction)

bars <- dcast(counts, order + label ~ Direction, value.var = "Proteins")
bars[, Total := Positive + Negative]
positive_col <- "#C46A76"
negative_col <- "#6489AE"

plot <- ggplot(bars, aes(x = order)) +
  geom_col(aes(y = Total, fill = "Positive"), width = 0.48) +
  geom_col(aes(y = Negative, fill = "Negative"), width = 0.48) +
  geom_text(aes(y = Total + 6, label = Positive),
            size = 4.7, family = "Arial", colour = "#252525") +
  geom_segment(aes(xend = order, y = Negative, yend = Negative + 11),
               colour = negative_col, linewidth = 0.35) +
  geom_label(aes(y = Negative + 17, label = Negative),
             size = 4.5, family = "Arial", colour = negative_col,
             fill = "white", linewidth = 0, label.padding = unit(0.08, "lines")) +
  scale_fill_manual(values = c("Negative" = negative_col,
                               "Positive" = positive_col),
                    breaks = c("Negative", "Positive"), name = NULL) +
  scale_x_continuous(breaks = bars$order, labels = bars$label,
                     expand = expansion(add = 0.55)) +
  scale_y_continuous(breaks = scales::breaks_pretty(n = 5),
                     limits = c(0, max(bars$Total) * 1.15),
                     expand = expansion(mult = c(0, 0))) +
  labs(x = NULL, y = "Number of proteins") +
  theme_classic(base_size = 14, base_family = "Arial") +
  theme(
    panel.grid.major.y = element_line(colour = "#E5E5E5", linewidth = 0.3),
    axis.text.x = element_text(size = 12, colour = "black", lineheight = 0.9),
    axis.text.y = element_text(size = 13, colour = "black"),
    axis.title.y = element_text(size = 14, colour = "black"),
    legend.text = element_text(size = 13),
    legend.position = "top",
    legend.justification = "right",
    legend.box.just = "right",
    legend.background = element_rect(fill = "white", colour = NA),
    legend.key.size = unit(3.5, "mm"),
    plot.margin = margin(7, 9, 7, 5)
  )

ggsave(file.path(figure_dir, "bpv_union_protein_direction_counts.pdf"),
       plot, width = 6.2, height = 3.8, device = cairo_pdf, bg = "white")
})
