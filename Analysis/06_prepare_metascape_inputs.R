# ==============================================================================
# METADATA
# ==============================================================================
# Author: Chaoyang Lin, MD
# Contact: chaoyanglint@163.com / Chaoyang Lin
#
# Dependencies:
#   R >= 4.1.0
#   Required packages: dplyr, data.table
#
# Goal:
#   Prepare SBPV and DBPV foreground gene lists with the retained Olink background.
#
# Required input data:
#   - ./BPV_ASCVD/data/bpv_protein_results.csv
#   - ./BPV_ASCVD/data/olink_gene_map.rds
#   - ./BPV_ASCVD/data/olink_proteins.rds
#
# Main outputs:
#   - ./BPV_ASCVD/data/metascape/metascape_bpv.txt
# ==============================================================================

library(data.table)
library(dplyr)

selected_proteins <- fread(
  "./BPV_ASCVD/data/bpv_protein_results.csv"
) |>
  filter(spec_id %in% c("SBPV_5y_VIM", "DBPV_5y_VIM"))
protein_gene_map <- readRDS("./BPV_ASCVD/data/olink_gene_map.rds")

dir.create(
  "./BPV_ASCVD/data/metascape",
  recursive = TRUE,
  showWarnings = FALSE
)

gene_list <- function(features) {
  protein_gene_map |>
    filter(feature %in% features, !is.na(gene_symbol), gene_symbol != "") |>
    pull(gene_symbol) |>
    unique() |>
    sort() |>
    paste(collapse = ",")
}

# Export the two trait lists in one Metascape analysis; the measured Olink assay universe is the enrichment background.
writeLines(
  c(
    "#Name\tGenes",
    paste0("SBPV\t", gene_list(selected_proteins$feature[selected_proteins$pressure == "SBPV"])),
    paste0("DBPV\t", gene_list(selected_proteins$feature[selected_proteins$pressure == "DBPV"])),
    paste0("_BACKGROUND\t", gene_list(readRDS("./BPV_ASCVD/data/olink_proteins.rds")))
  ),
  "./BPV_ASCVD/data/metascape/metascape_bpv.txt"
)

