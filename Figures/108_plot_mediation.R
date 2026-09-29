# ==============================================================================
# METADATA
# ==============================================================================
# Author: Chaoyang Lin, MD
# Goal: Draw overall-MPS and representative-protein mediation figures.
# Input: Saved analysis outputs in ./BPV_ASCVD.
# Outputs: PDFs and figure tables under ./BPV_ASCVD.
# ==============================================================================

# --- mps_mediation_cox ---
local({
library(ggplot2)
library(patchwork)

root <- "./BPV_ASCVD"
results <- read.csv(file.path(root, "mediation/mps_cox/mediation_summary.csv"),
                    check.names = FALSE)
models <- data.frame(Trait = c("SBPV", "SBPV", "SBPV", "DBPV", "DBPV"),
                     Outcome = c("MI", "IS", "PAD", "IS", "PAD"))
colours <- c(SBPV = "#A6611A", DBPV = "#5875A4")

# A single neutral cell explains the diagram used for each result.
schematic <- ggplot() +
  annotate("text", x = 0.3, y = 9.55, hjust = 0,
           label = "Mediation model", size = 5.2, fontface = "bold") +
  annotate("text", x = 0.3, y = 8.75, hjust = 0,
           label = "HR scale (Per 1-SD increase)", size = 3.9) +
  annotate("segment", x = 1.9, y = 3.7, xend = 4.7, yend = 6.1,
           colour = "#444444", linewidth = 0.7,
           arrow = grid::arrow(length = grid::unit(2.3, "mm"))) +
  annotate("segment", x = 5.4, y = 6.1, xend = 8.1, yend = 3.7,
           colour = "#444444", linewidth = 0.7,
           arrow = grid::arrow(length = grid::unit(2.3, "mm"))) +
  annotate("segment", x = 2.55, y = 3.1, xend = 7.45, yend = 3.1,
           colour = "#444444", linewidth = 0.6,
           arrow = grid::arrow(length = grid::unit(2.3, "mm"))) +
  annotate("rect", xmin = 3.2, xmax = 6.8, ymin = 6.12, ymax = 7.06,
           fill = "white", colour = "#444444", linewidth = 0.4) +
  annotate("rect", xmin = 0.1, xmax = 2.4, ymin = 2.6, ymax = 3.65,
           fill = "white", colour = "#444444", linewidth = 0.4) +
  annotate("rect", xmin = 7.6, xmax = 9.9, ymin = 2.6, ymax = 3.65,
           fill = "white", colour = "#444444", linewidth = 0.4) +
  annotate("text", x = 5, y = 6.59, label = "Mediator", size = 4.2,
           fontface = "bold") +
  annotate("text", x = 1.25, y = 3.12, label = "Exposure", size = 3.7,
           fontface = "bold") +
  annotate("text", x = 8.75, y = 3.12, label = "Outcome", size = 3.7,
           fontface = "bold") +
  annotate("text", x = 5, y = 7.75, label = "IE: Indirect effect",
           size = 4.0, colour = "#444444") +
  annotate("text", x = 5, y = 3.55, label = "DE: Direct effect",
           size = 4.0, colour = "#444444") +
  annotate("text", x = 5, y = 2.24, label = "TE: Total effect",
           size = 3.9, colour = "#444444") +
  coord_cartesian(xlim = c(0, 10), ylim = c(1.55, 9.9), clip = "off") +
  theme_void(base_family = "sans") +
  theme(plot.background = element_rect(fill = "white", colour = NA),
        plot.margin = margin(3, 6, 3, 6))

# Draw one model per cell, using effect HRs and 95% CIs from script 10.
effect_panel <- function(trait, outcome) {
  d <- results[results$Trait == trait & results$Outcome == outcome, ]
  effect <- function(name, label) {
    x <- d[d$Effect == name, ]
    sprintf("%s: %.3f (%.3f–%.3f)", label,
            x$Estimate, x$CI_low, x$CI_high)
  }
  colour <- unname(colours[trait])
  box_fill <- if (trait == "SBPV") "#F3DFD2" else "#DFE8F4"

  ggplot() +
    annotate("text", x = 5, y = 8.58, label = effect("Indirect_HR", "IE"),
             size = 3.95, fontface = "bold", colour = "#252525") +
    annotate("text", x = 5, y = 7.85, label = sprintf("Proportion mediated: %.1f%%", d$Estimate[d$Effect == "PM_percent"]),
             size = 3.95, colour = "#444444") +
    annotate("segment", x = 1.9, y = 3.73, xend = 4.65, yend = 6.00,
             colour = colour, linewidth = 0.75,
             arrow = grid::arrow(length = grid::unit(2.3, "mm"))) +
    annotate("segment", x = 5.38, y = 6.00, xend = 8.1, yend = 3.73,
             colour = colour, linewidth = 0.75,
             arrow = grid::arrow(length = grid::unit(2.3, "mm"))) +
    annotate("segment", x = 2.35, y = 3.05, xend = 7.65, yend = 3.05,
             colour = "#555B62", linewidth = 0.6,
             arrow = grid::arrow(length = grid::unit(2.3, "mm"))) +
    annotate("rect", xmin = 3.50, xmax = 6.50, ymin = 6.00, ymax = 7.15,
             fill = box_fill, colour = NA) +
    annotate("rect", xmin = 0.15, xmax = 2.25, ymin = 2.52, ymax = 3.73,
             fill = box_fill, colour = NA) +
    annotate("rect", xmin = 7.75, xmax = 9.85, ymin = 2.52, ymax = 3.73,
             fill = box_fill, colour = NA) +
    annotate("text", x = 5, y = 6.59, label = paste0("MPS-", trait),
             size = 3.95, fontface = "bold") +
    annotate("text", x = 1.25, y = 3.12, label = trait,
             size = 4.15, fontface = "bold") +
    annotate("text", x = 8.80, y = 3.12, label = outcome,
             size = 4.15, fontface = "bold") +
    annotate("text", x = 5, y = 3.55, label = effect("Direct_HR", "DE"),
             size = 3.75, colour = "#252525") +
    annotate("text", x = 5, y = 2.24, label = effect("Total_HR", "TE"),
             size = 3.75, fontface = "bold", colour = "#252525") +
    coord_cartesian(xlim = c(0, 10), ylim = c(1.55, 9.9), clip = "off") +
    theme_void(base_family = "sans") +
    theme(plot.background = element_rect(fill = "white", colour = NA),
          plot.margin = margin(3, 6, 3, 6))
}

figure <- wrap_plots(
  c(list(schematic), Map(effect_panel, models$Trait, models$Outcome)), ncol = 2, nrow = 3
)
ggsave(file.path(root, "figures/overall_mps_mediation_regmedint_cox.pdf"), figure, width = 8.0, height = 8.1, device = cairo_pdf)
})

# --- representative_protein_mediation ---
local({

suppressPackageStartupMessages(library(ggplot2))

root <- "./BPV_ASCVD"
results <- read.csv(file.path(root, "mediation/functional_representative_proteins/results.csv"),
                    check.names = FALSE)

# Preserve the functional-group order.
proteins <- unique(results[, c("Protein", "belongs_to")])
proteins$group_number <- as.integer(sub("F", "", sub(",.*", "", proteins$belongs_to)))
proteins <- proteins[order(proteins$group_number, proteins$Protein), ]
proteins$y <- rev(seq_len(nrow(proteins)))

columns <- data.frame(
  Trait = c("SBPV", "SBPV", "SBPV", "DBPV", "DBPV"),
  Outcome = c("MI", "IS", "PAD", "IS", "PAD"),
  x = c(0.8, 1.6, 2.4, 3.45, 4.25)
)
plot_data <- merge(results, proteins[, c("Protein", "y")], by = "Protein")
plot_data <- merge(plot_data, columns, by = c("Trait", "Outcome"))
selected <- plot_data[is.finite(plot_data$PM_p) & plot_data$PM_p < 0.05, ]

# Empty cells retain the row and outcome positions, as in the original bubble map.
figure <- ggplot() +
  geom_rect(data = proteins[proteins$y %% 2L == 0L, ],
            aes(xmin = -1.50, xmax = 4.70, ymin = y - 0.48, ymax = y + 0.48),
            fill = "#F5F7F8", colour = NA) +
  annotate("rect", xmin = 0.45, xmax = 2.75,
           ymin = 0.5, ymax = 16.5,
           fill = NA, colour = "#73808A", linewidth = 0.45) +
  annotate("rect", xmin = 3.10, xmax = 4.60,
           ymin = 0.5, ymax = 16.5,
           fill = NA, colour = "#73808A", linewidth = 0.45) +
  geom_text(data = proteins, aes(x = -1.43, y = y, label = Protein),
            hjust = 0, family = "sans", size = 4.15, colour = "#26323A") +
  geom_text(data = proteins, aes(x = -0.33, y = y, label = belongs_to),
            hjust = 0, family = "sans", size = 4.15, colour = "#26323A") +
  geom_point(data = selected,
             aes(x = x, y = y, size = abs(PM_percent), fill = -log10(PM_p)),
             shape = 21, stroke = 0.35, colour = "#71575B") +
  annotate("text", x = -1.43, y = 17.15, hjust = 0,
           label = "Protein", fontface = "bold", size = 4.15) +
  annotate("text", x = -0.33, y = 17.15, hjust = 0,
           label = "Belongs to", fontface = "bold", size = 4.15) +
  annotate("text", x = 1.6, y = 18.10, label = "SBPV",
           fontface = "bold", size = 5.10) +
  annotate("text", x = 3.85, y = 18.10, label = "DBPV",
           fontface = "bold", size = 5.10) +
  annotate("text", x = columns$x, y = 17.15,
           label = c("MI", "IS", "PAD", "IS", "PAD"),
           fontface = "bold", size = 4.25, lineheight = 0.95) +
  scale_size_area(name = "Proportion mediated (%)",
                  max_size = 9, limits = c(0, max(20, abs(selected$PM_percent))),
                  breaks = c(5, 10, 15, 20),
                  guide = guide_legend(order = 1, title.position = "top",
                                       direction = "horizontal", nrow = 1,
                                       override.aes = list(fill = "#D8C3C4"))) +
  scale_fill_gradient(name = expression(-log[10](P[PM])),
                      low = "#F3DFDC", high = "#B35D65",
                      guide = guide_colourbar(order = 2,
                                              barwidth = grid::unit(42, "mm"),
                                              barheight = grid::unit(4, "mm"),
                                              title.position = "top")) +
  coord_cartesian(xlim = c(-1.50, 4.70), ylim = c(0.35, 18.3), clip = "off") +
  labs(x = NULL, y = NULL) +
  theme_void(base_family = "sans") +
  theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.spacing.x = grid::unit(10, "mm"),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 11.5),
    plot.margin = margin(2, 2, 2, 2),
    plot.background = element_rect(fill = "white", colour = NA)
  )

ggsave(file.path(root, "figures/functional_representative_protein_mediation.pdf"), figure, width = 7.5, height = 8.0,
       device = cairo_pdf, bg = "white")
})
