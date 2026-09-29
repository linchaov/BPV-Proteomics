# ==============================================================================
# METADATA
# ==============================================================================
# Author: Chaoyang Lin, MD
#
# Goal: Draw Metascape enrichment and functional network figures.
# Input: Saved analysis outputs in ./BPV_ASCVD.
# Outputs: PDFs and figure tables under ./BPV_ASCVD.
# ==============================================================================

# --- metascape_enrichment ---
local({
library(data.table)
library(dplyr)
library(pheatmap)

source_number <- c(
  "Hallmark Gene Sets" = 1,
  "KEGG Pathway" = 2,
  "GO Biological Processes" = 3,
  "Reactome Gene Sets" = 4
)

heatmap_file <- "./BPV_ASCVD/data/metascape/SBPV+DBPV/all.tpgra8hic/Enrichment_heatmap/HeatmapSelectedGO"

terms <- fread(paste0(heatmap_file, ".csv")) |>
  left_join(
    fread("./BPV_ASCVD/data/metascape/SBPV+DBPV/all.tpgra8hic/Enrichment_GO/_FINAL_GO.csv") |>
      distinct(GO, Category),
    by = "GO"
  ) |>
  mutate(label = paste0("(", unname(source_number[Category]), ") ", Description))

cdt <- fread(paste0(heatmap_file, ".cdt")) |>
  filter(grepl("^GENE", GID)) |>
  left_join(terms |> select(GO, label), by = c("GENE" = "GO"))

# Convert Metascape's JTreeView tree nodes into an R dendrogram without reclustering.
read_tree <- function(file, leaf_ids, labels) {
  x <- fread(file, header = FALSE, col.names = c("node", "left", "right", "similarity"))
  # Resolve internal node IDs and leaf IDs in Metascape's exported tree.
  index <- function(z) if (startsWith(z, "NODE")) match(z, x$node) else -match(z, leaf_ids)
  merge <- cbind(vapply(x$left, index, integer(1)), vapply(x$right, index, integer(1)))
  leaves <- function(i) unlist(lapply(merge[i, ], function(j) if (j < 0) -j else leaves(j)))
  structure(
    list(merge = merge, height = 1 - x$similarity, order = leaves(nrow(x)),
         labels = labels, method = "Metascape export", call = match.call(), dist.method = "export"),
    class = "hclust"
  )
}

mat <- as.matrix(terms[match(cdt$GENE, terms$GO), c("_LogP_DBPV", "_LogP_SBPV")])
storage.mode(mat) <- "numeric"
# A zero in the export denotes no enrichment; it is displayed in gray.
mat[mat == 0] <- NA_real_
mat <- -mat
dimnames(mat) <- list(cdt$label, c("DBPV", "SBPV"))

row_tree <- read_tree(paste0(heatmap_file, ".gtr"), cdt$GID, cdt$label)
row_tree <- as.hclust(rev(as.dendrogram(row_tree)))
col_tree <- read_tree(paste0(heatmap_file, ".atr"), c("ARRY1X", "ARRY2X"), colnames(mat))

metascape_colours <- scales::gradient_n_pal(
  c("#FFFFFF", "#FFFED1", "#FEE492", "#FEC44F", "#FE9829", "#D9600E", "#993404"),
  values = c(0, 2, 3, 4, 6, 10, 20) / 20
)(seq(0, 1, length.out = 200))

heatmap_plot <- pheatmap(
  mat, cluster_rows = row_tree, cluster_cols = col_tree,
  color = metascape_colours, breaks = seq(0, 20, length.out = 201), na_col = "#D9D9D9",
  legend = FALSE,
  treeheight_row = 85, treeheight_col = 45, cellwidth = 78, cellheight = 23,
  fontsize_row = 12, fontsize_col = 13, angle_col = 270, border_color = NA,
  silent = TRUE
)

grDevices::cairo_pdf(
  "./BPV_ASCVD/figures/metascape_enrichment.pdf",
  width = 14.2, height = 9
)
grid::grid.newpage()
grid::grid.draw(heatmap_plot$gtable)
grid::grid.raster(
  matrix(metascape_colours, nrow = 1), x = 0.50, y = 0.91,
  width = grid::unit(0.20, "npc"), height = grid::unit(0.028, "npc"), interpolate = TRUE
)
grid::grid.text(expression(-log[10](P)), x = 0.50, y = 0.952, gp = grid::gpar(fontsize = 14))
for (value in c(0, 2, 3, 4, 6, 10, 20)) {
  grid::grid.text(value, x = 0.40 + 0.20 * value / 20, y = 0.875,
                  gp = grid::gpar(fontsize = 10))
}
grid::grid.text(
  "(1) MSigDB Hallmark   (2) KEGG   (3) GO Biological Processes   (4) Reactome",
  x = 0.10, y = 0.035, just = "left", gp = grid::gpar(fontsize = 11)
)
grDevices::dev.off()

})

# --- metascape_functional_network ---
local({

library(xml2)
library(jsonlite)
library(dplyr)
library(ggplot2)
library(ggforce)

network <- read_xml(
  "./BPV_ASCVD/data/metascape/SBPV+DBPV/all.tpgra8hic/Enrichment_GO/GONetwork.xgmml"
)
node_xml <- xml_find_all(network, ".//*[local-name()='node']")
edge_xml <- xml_find_all(network, ".//*[local-name()='edge']")
# Extract one named Cytoscape node or edge attribute from the Metascape network.
attribute <- function(x, name) xml_attr(xml_find_first(x, paste0("./*[@name='", name, "']")), "value")

nodes <- data.frame(
  id = xml_attr(node_xml, "id"),
  group = as.integer(vapply(node_xml, attribute, character(1), "GROUP_ID")),
  proteins = as.integer(vapply(node_xml, attribute, character(1), "#GeneInGOAndHitList")),
  dbpv = as.integer(vapply(node_xml, attribute, character(1), "_MEMBER_DBPV")),
  sbpv = as.integer(vapply(node_xml, attribute, character(1), "_MEMBER_SBPV"))
) |>
  mutate(
    cluster = factor(group),
    enrichment = case_when(
      dbpv == 1L & sbpv == 1L ~ "Both",
      dbpv == 1L ~ "DBPV only",
      TRUE ~ "SBPV only"
    )
  )

# Keep Metascape's group IDs: these are the F01-F20 labels reused by Scripts.

layout_json <- fromJSON(
  sub(";$", "", sub("^var networks=", "", readLines(
    "./BPV_ASCVD/data/metascape/SBPV+DBPV/all.tpgra8hic/Enrichment_GO/GONetwork.js",
    warn = FALSE
  ))), simplifyVector = FALSE
)$GONetwork$elements$nodes
layout <- data.frame(
  id = vapply(layout_json, function(z) z$data$id_original, character(1)),
  x = vapply(layout_json, function(z) z$position$x, numeric(1)),
  y = -vapply(layout_json, function(z) z$position$y, numeric(1))
) |>
  mutate(
    x = if_else(x > 450, 450 + 0.65 * (x - 450), x),
    y = if_else(y < -600, -600 + 0.55 * (y + 600), y)
  )
nodes <- left_join(nodes, layout, by = "id")

edges <- data.frame(
  source = xml_attr(edge_xml, "source"),
  target = xml_attr(edge_xml, "target"),
  score = as.numeric(vapply(edge_xml, attribute, character(1), "SCORE"))
) |>
  mutate(
    x = nodes$x[match(source, nodes$id)], y = nodes$y[match(source, nodes$id)],
    xend = nodes$x[match(target, nodes$id)], yend = nodes$y[match(target, nodes$id)]
  )

# Metascape representative terms, with source prefixes removed and line breaks added.
cluster_labels <- c(
  "Vasculature\ndevelopment", "Antimicrobial humoral\nresponse",
  "Cytokine–cytokine receptor\ninteraction", "Locomotion",
  "Enzyme-linked receptor\nprotein signaling pathway",
  "Regulation of\ncell activation",
  "Regulation of cellular\nresponse to growth factor\nstimulus",
  "Animal organ\nmorphogenesis", "Renal system\nvasculature\ndevelopment",
  "Regulation of Insulin-like Growth Factor\n(IGF) transport and uptake by Insulin-like\nGrowth Factor Binding Proteins (IGFBPs)",
  "Epithelial–mesenchymal\ntransition", "IL6 JAK STAT3\nsignaling",
  "Activation of NF-kappaB-\ninducing kinase activity",
  "Positive regulation\nof MAPK cascade",
  "Regulation of T-cell\nreceptor signaling pathway",
  "Positive regulation of\nalpha-beta T-cell\nproliferation",
  "Extracellular matrix\norganization", "Alpha-defensins",
  "Inflammatory response",
  "Regulation of vascular\nendothelial growth factor\nsignaling pathway"
)
# Label positions and leader lines are fixed for the final layout.
cluster_centres <- nodes |>
  group_by(group) |>
  summarise(x = mean(x), y = mean(y), .groups = "drop") |>
  mutate(cluster_label = paste0(sprintf("F%02d.", group), cluster_labels[group])) |>
  left_join(data.frame(
    group = 1:20,
    label_x = c(-5, 120, 1410, -430, -405, 760, 255, 180, -455, 75,
                1200, 1540, 540, -435, 370, 760, 1040, -455, -60, 155),
    label_y = c(-135, -850, -300, -625, -210, 90, -255, -700, -375, -1210,
                135, 95, -640, -510, -390, -310, -185, -910, -1000, 135)
  ), by = "group") |>
  mutate(
    label_x = if_else(label_x > 450, 450 + 0.65 * (label_x - 450), label_x),
    label_y = if_else(label_y < -600, -600 + 0.55 * (label_y + 600), label_y)
  )

# Metascape hues identify functional term groups; BPV-list identity instead
# uses the Manhattan/scatter palette on node outlines.
cluster_colours <- c(
  "#D96D66", "#5C8DBA", "#6EA878", "#8E78AA", "#D99A53",
  "#DBBD58", "#AA8068", "#D492B5", "#9BABB5", "#72B9AF",
  "#D5CA8D", "#A9A6CB", "#E49A8E", "#7EAAC7", "#DBAC7C",
  "#A6BE74", "#DCAFC4", "#B6B6BB", "#A884B3", "#99C2A3"
)
names(cluster_colours) <- as.character(seq_along(cluster_colours))

# Singleton terms retain their nodes and edges but do not need an ellipse.
functional_network <- ggplot() +
  geom_mark_ellipse(
    data = nodes |> add_count(cluster) |> filter(n > 1),
    aes(x, y, group = cluster, fill = cluster),
    alpha = 0.055, colour = "#777777", linewidth = 0.3,
    expand = grid::unit(2.5, "mm"), show.legend = FALSE
  ) +
  geom_segment(
    data = edges, aes(x, y, xend = xend, yend = yend, linewidth = score),
    colour = "#694C9B", alpha = 0.24, show.legend = FALSE
  ) +
  geom_point(
    data = nodes,
    aes(x, y, fill = cluster, shape = enrichment, colour = enrichment, size = proteins),
    stroke = 0.8
  ) +
  geom_segment(
    data = filter(cluster_centres, group %in% c(4L, 8L, 10L, 11L, 15L)),
    aes(x, y, xend = label_x, yend = label_y),
    colour = "#777777", linewidth = 0.3
  ) +
  geom_label(
    data = cluster_centres, aes(label_x, label_y, label = cluster_label),
    size = 3.5, family = "sans", lineheight = 0.9,
    fill = "white", linewidth = 0, label.padding = grid::unit(0, "mm")
  ) +
  scale_fill_manual(values = cluster_colours, guide = "none") +
  scale_colour_manual(
    values = c("SBPV only" = "#A6611A", "DBPV only" = "#5875A4", "Both" = "#7B4F7C"),
    guide = "none"
  ) +
  scale_shape_manual(
    name = "Term enriched in", breaks = c("SBPV only", "DBPV only", "Both"),
    labels = c("SBPV", "DBPV", "Both"),
    values = c("SBPV only" = 22, "DBPV only" = 24, "Both" = 21),
    guide = guide_legend(override.aes = list(
      fill = "white", size = 4.5,
      colour = c("#A6611A", "#5875A4", "#7B4F7C")
    ))
  ) +
  scale_size_continuous(range = c(2.8, 6), guide = "none") +
  scale_linewidth(range = c(0.25, 0.9)) +
  coord_equal(xlim = c(-500, 1240), ylim = c(-980, 180), clip = "off") +
  labs(x = NULL, y = NULL) +
  theme_void(base_size = 14, base_family = "sans") +
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.88, 0.09),
    legend.background = element_rect(fill = "white", colour = NA),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9),
    legend.key.size = grid::unit(5.5, "mm"),
    plot.margin = margin(8, 8, 8, 8)
  )

ggsave(
  "./BPV_ASCVD/figures/metascape_network.pdf",
  functional_network, width = 11.5, height = 7.7, device = cairo_pdf, bg = "white"
)

})

