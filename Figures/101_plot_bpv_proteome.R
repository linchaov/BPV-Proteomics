# ==============================================================================
# METADATA
# ==============================================================================
# Author: Chaoyang Lin, MD
# Goal: Compare SBPV/DBPV protein effects and draw Manhattan plots.
# Input: Saved analysis outputs in ./BPV_ASCVD.
# Outputs: PDFs and figure tables under ./BPV_ASCVD.
# ==============================================================================

# --- sbpv_dbpv_protein_effects ---
local({

library(data.table)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(dplyr)

selected_proteins <- fread(
  "./BPV_ASCVD/data/bpv_protein_results.csv"
) |>
  filter(spec_id %in% c("SBPV_5y_VIM", "DBPV_5y_VIM"))
selected_sets <- split(selected_proteins$feature, selected_proteins$pressure)[c("SBPV", "DBPV")]

# The scatter uses all tested assays; colors and overlap depend on the prespecified protein screen
effect_comparison <- readRDS(
  "./BPV_ASCVD/data/bpv_protein_plot_data.rds"
) |>
  select(feature, pressure, beta, fdr, passes_protein_screen) |>
  pivot_wider(
    names_from = pressure,
    values_from = c(beta, fdr, passes_protein_screen),
    names_glue = "{.value}_{pressure}"
  ) |>
  filter(
    is.finite(beta_SBPV), is.finite(beta_DBPV),
    is.finite(fdr_SBPV), is.finite(fdr_DBPV),
    !is.na(passes_protein_screen_SBPV), !is.na(passes_protein_screen_DBPV)
  ) |>
  mutate(
    signature = case_when(
      passes_protein_screen_SBPV & passes_protein_screen_DBPV ~ "Both",
      passes_protein_screen_SBPV ~ "SBPV only",
      passes_protein_screen_DBPV ~ "DBPV only",
      TRUE ~ "Neither"
    ),
    signature = factor(signature, levels = c("Neither", "SBPV only", "DBPV only", "Both")),
    minimum_fdr = pmin(fdr_SBPV, fdr_DBPV)
  ) |>
  arrange(signature)

label_data <- bind_rows(
  effect_comparison |>
    filter(signature != "Neither", beta_SBPV > 0, beta_DBPV > 0) |>
    arrange(minimum_fdr, desc(beta_SBPV + beta_DBPV)) |>
    slice_head(n = 4),
  effect_comparison |>
    filter(signature != "Neither", beta_SBPV < 0, beta_DBPV < 0) |>
    arrange(minimum_fdr, desc(abs(beta_SBPV) + abs(beta_DBPV))) |>
    slice_head(n = 3)
) |>
  distinct(feature, .keep_all = TRUE)

correlation_label <- sprintf(
  "Pearson r = %.3f\nSpearman rho = %.3f\nN = %s",
  cor(effect_comparison$beta_SBPV, effect_comparison$beta_DBPV),
  cor(effect_comparison$beta_SBPV, effect_comparison$beta_DBPV, method = "spearman"),
  format(nrow(effect_comparison), big.mark = ",")
)
axis_limit <- 1.16 * max(abs(c(effect_comparison$beta_SBPV, effect_comparison$beta_DBPV)))

circle_angle <- seq(0, 2 * pi, length.out = 361)
overlap_shapes <- bind_rows(
  data.frame(
    set = "SBPV", x = -0.72 + 1.18 * cos(circle_angle),
    y = 1.18 * sin(circle_angle)
  ),
  data.frame(
    set = "DBPV", x = 0.72 + 1.18 * cos(circle_angle),
    y = 1.18 * sin(circle_angle)
  )
)
overlap_counts <- data.frame(
  x = c(-1.05, 0, 1.05), y = 0,
  label = c(
    length(setdiff(selected_sets$SBPV, selected_sets$DBPV)),
    length(intersect(selected_sets$SBPV, selected_sets$DBPV)),
    length(setdiff(selected_sets$DBPV, selected_sets$SBPV))
  )
)

overlap_plot <- ggplot(overlap_shapes, aes(x, y, group = set, fill = set)) +
  geom_polygon(alpha = 0.60, colour = "grey40", linewidth = 0.6) +
  geom_text(
    data = overlap_counts,
    aes(x, y, label = label),
    inherit.aes = FALSE, size = 8.5, family = "sans", fontface = "bold"
  ) +
  annotate("text", x = -1.15, y = 1.42, label = "SBPV",
           size = 7.5, family = "sans", fontface = "bold") +
  annotate("text", x = 1.15, y = 1.42, label = "DBPV",
           size = 7.5, family = "sans", fontface = "bold") +
  scale_fill_manual(values = c(SBPV = "#A6611A", DBPV = "#5875A4")) +
  coord_equal(xlim = c(-2.0, 2.0), ylim = c(-1.30, 1.55), clip = "off") +
  theme_void(base_size = 12) +
  theme(legend.position = "none", plot.margin = margin(3, 3, 3, 3))

effect_plot <- ggplot(
  effect_comparison,
  aes(beta_SBPV, beta_DBPV)
) +
  geom_hline(yintercept = 0, colour = "grey80", linewidth = 0.35) +
  geom_vline(xintercept = 0, colour = "grey80", linewidth = 0.35) +
  geom_abline(
    intercept = 0, slope = 1, colour = "grey55",
    linetype = "dashed", linewidth = 0.5
  ) +
  geom_point(
    data = filter(effect_comparison, signature == "Neither"),
    aes(colour = signature), size = 1.2, alpha = 0.48
  ) +
  geom_point(
    data = filter(effect_comparison, signature != "Neither"),
    aes(colour = signature), size = 1.6, alpha = 0.88
  ) +
  geom_text_repel(
    data = label_data,
    aes(label = feature),
    family = "sans", size = 5.2, colour = "black", seed = 20260713L,
    max.overlaps = Inf, min.segment.length = 0,
    box.padding = 0.48, point.padding = 0.25,
    force = 2, max.time = 2, max.iter = 20000,
    segment.colour = "grey55", segment.size = 0.3,
    show.legend = FALSE
  ) +
  annotate(
    "text", x = -0.95 * axis_limit, y = 0.95 * axis_limit,
    label = correlation_label, family = "sans", size = 5.4,
    hjust = 0, vjust = 1
  ) +
  scale_colour_manual(
    values = c(
      "Neither" = "grey75",
      "SBPV only" = "#A6611A",
      "DBPV only" = "#5875A4",
      "Both" = "#7B4F7C"
    ),
    breaks = c("Both", "SBPV only", "DBPV only", "Neither"),
    name = NULL,
    guide = guide_legend(ncol = 1, override.aes = list(size = 4, alpha = 1))
  ) +
  scale_x_continuous(
    limits = c(-axis_limit, axis_limit),
    breaks = scales::breaks_pretty(n = 6),
    expand = expansion(mult = 0)
  ) +
  scale_y_continuous(
    limits = c(-axis_limit, axis_limit),
    breaks = scales::breaks_pretty(n = 6),
    expand = expansion(mult = 0)
  ) +
  labs(
    x = expression("SBPV–protein association (" * beta * ")"),
    y = expression("DBPV–protein association (" * beta * ")")
  ) +
  coord_equal(clip = "off") +
  theme_bw(base_family = "sans") +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(colour = "#48545C", linewidth = 0.5),
    axis.title = element_text(size = 19, colour = "black"),
    axis.text = element_text(size = 15, colour = "black"),
    axis.ticks = element_line(colour = "#48545C", linewidth = 0.4),
    legend.text = element_text(size = 15),
    legend.position = c(0.98, 0.02),
    legend.justification = c(1, 0),
    legend.background = element_rect(fill = "white", colour = NA),
    legend.margin = margin(3, 3, 3, 3),
    plot.margin = margin(8, 10, 6, 8)
  )

ggsave(
  "./BPV_ASCVD/figures/bpv_protein_overlap.pdf",
  overlap_plot, width = 5.5, height = 4.2, device = cairo_pdf
)
ggsave(
  "./BPV_ASCVD/figures/bpv_protein_effects.pdf",
  effect_plot, width = 7.4, height = 6.4, device = cairo_pdf
)
})

# --- bpv_proteome_manhattan ---
local({

library(ggplot2)
library(ggrepel)
library(dplyr)

protein_associations <- readRDS(
  "./BPV_ASCVD/data/bpv_protein_plot_data.rds"
)
assay_genes <- readRDS(
  "./BPV_ASCVD/data/olink_gene_map.rds"
) |>
  filter(!is.na(entrez_id), entrez_id != "") |>
  distinct(feature, gene_symbol, entrez_id) |>
  group_by(feature) |>
  filter(n_distinct(entrez_id) == 1L) |>
  slice_head(n = 1L) |>
  ungroup()

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
gene_ranges <- suppressMessages(
  GenomicFeatures::genes(txdb, single.strand.genes.only = TRUE)
)
gene_coordinates <- as.data.frame(gene_ranges) |>
  mutate(entrez_id = as.character(names(gene_ranges))) |>
  transmute(
    entrez_id,
    chromosome = sub("^chr", "", as.character(seqnames)),
    gene_position = (as.numeric(start) + as.numeric(end)) / 2
  )

chr_levels <- c(as.character(1:22), "X", "Y")
chromosome_lengths <- GenomeInfoDb::seqlengths(txdb)
chromosome_layout <- data.frame(
  chromosome = sub("^chr", "", names(chromosome_lengths)),
  chromosome_length = as.numeric(chromosome_lengths)
) |>
  filter(chromosome %in% chr_levels, !is.na(chromosome_length)) |>
  mutate(chromosome_order = match(chromosome, chr_levels)) |>
  arrange(chromosome_order) |>
  mutate(
    offset = lag(cumsum(chromosome_length), default = 0),
    center = offset + chromosome_length / 2,
    xmax = offset + chromosome_length
  )

assay_coordinates <- assay_genes |>
  inner_join(gene_coordinates, by = "entrez_id") |>
  filter(chromosome %in% chr_levels) |>
  mutate(chromosome_order = match(chromosome, chr_levels)) |>
  arrange(chromosome_order, gene_position) |>
  distinct(feature, .keep_all = TRUE)

# Retain all mappable assays for the gray background, then highlight 
# proteins passing the complete screen. Coordinates are hg38 gene midpoints.
manhattan_data <- protein_associations |>
  filter(
    is.finite(p_value), is.finite(fdr)
  ) |>
  inner_join(assay_coordinates, by = "feature") |>
  left_join(
    chromosome_layout |> select(chromosome, offset),
    by = "chromosome"
  ) |>
  mutate(
    cumulative_position = offset + gene_position,
    neg_log10_p = -log10(pmax(p_value, .Machine$double.xmin))
  ) |>
  arrange(pressure, chromosome_order, gene_position, feature)

# Plot genome-ordered protein associations for one BPV trait.
plot_manhattan <- function(data, trait, sig_color) {
  data <- filter(data, pressure == trait)
  tested <- protein_associations |>
    filter(pressure == trait, is.finite(p_value), is.finite(fdr))
  fdr_hits <- sum(tested$fdr < 0.05)
  fdr_y <- if (fdr_hits > 0L) {
    -log10(0.05 * fdr_hits / nrow(tested))
  } else {
    NA_real_
  }
  label_data <- data |>
    filter(passes_protein_screen) |>
    arrange(p_value, desc(abs(beta))) |>
    slice_head(n = 20L)

  ggplot(data, aes(cumulative_position, neg_log10_p)) +
    geom_vline(
      data = filter(chromosome_layout, chromosome_order < max(chromosome_order)),
      aes(xintercept = xmax), inherit.aes = FALSE,
      colour = "grey75", linetype = "dashed", linewidth = 0.3
    ) +
    geom_hline(
      yintercept = -log10(0.05), colour = "grey55",
      linetype = "twodash", linewidth = 0.6
    ) +
    {if (is.finite(fdr_y)) geom_hline(
      yintercept = fdr_y, colour = sig_color,
      linetype = "dashed", linewidth = 0.7
    )} +
    geom_point(
      data = filter(data, !passes_protein_screen),
      colour = "grey65", size = 1.5, alpha = 0.85
    ) +
    geom_point(
      data = filter(data, passes_protein_screen),
      colour = sig_color, size = 2.1
    ) +
    geom_text_repel(
      data = label_data, aes(label = feature),
      family = "sans", size = 7, colour = sig_color,
      seed = 20260713L, max.overlaps = Inf,
      min.segment.length = 0, box.padding = 0.35,
      point.padding = 0.20, segment.colour = sig_color,
      show.legend = FALSE
    ) +
    scale_x_continuous(
      breaks = chromosome_layout$center,
      labels = chromosome_layout$chromosome,
      limits = c(0, max(chromosome_layout$xmax) * 1.01),
      expand = expansion(mult = c(0.005, 0.005))
    ) +
    scale_y_continuous(
      breaks = scales::breaks_pretty(n = 7),
      expand = expansion(mult = c(0.02, 0.16))
    ) +
    labs(
      title = trait, x = "Chromosome",
      y = expression(-log[10](italic(P)))
    ) +
    coord_cartesian(clip = "off") +
    theme_bw(base_family = "sans") +
    theme(
      panel.grid = element_blank(),
      plot.title = element_text(face = "bold", size = 20, hjust = 0.5),
      axis.title = element_text(size = 20),
      axis.text = element_text(size = 18, colour = "black"),
      axis.text.x = element_text(size = 15),
      plot.margin = margin(12, 20, 8, 12)
    )
}

ggsave(
  "./BPV_ASCVD/figures/sbpv_manhattan.pdf",
  plot_manhattan(manhattan_data, "SBPV", "#A6611A"),
  width = 14.5, height = 5.5, device = cairo_pdf
)
ggsave(
  "./BPV_ASCVD/figures/dbpv_manhattan.pdf",
  plot_manhattan(manhattan_data, "DBPV", "#5875A4"),
  width = 14.5, height = 5.5, device = cairo_pdf
)
})

