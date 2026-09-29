# ==============================================================================
# METADATA
# ==============================================================================
# Author: Chaoyang Lin, MD
# Dependencies: R >= 4.1.0; xml2, data.table, dplyr, tidyr, fst,
#               survival, broom
# Goal: Select one representative protein per Metascape functional group and
#       estimate its association with each incident cardiovascular outcome.
# Inputs: Metascape GONetwork.xgmml; olink_gene_map.rds;
#         bpv_protein_results.csv; cohort_olink.fst.
# Output: tables/supp_functional_proteins.csv.
# Notes: A protein can represent multiple F groups but is fitted only once
#        per outcome. Run after the combined SBPV-DBPV Metascape export.
# ==============================================================================

library(xml2)
library(data.table)
library(dplyr)
library(tidyr)

# For each Metascape functional group, choose one protein from the union of
# selected SBPV and DBPV proteins. Rank by term frequency, then BPV-PWAS FDR.
network <- read_xml(
  "./BPV_ASCVD/data/metascape/SBPV+DBPV/all.tpgra8hic/Enrichment_GO/GONetwork.xgmml"
)
node_xml <- xml_find_all(network, ".//*[local-name()='node']")
# Extract one named Cytoscape attribute from the Metascape network.
attribute <- function(x, name) xml_attr(xml_find_first(x, paste0("./*[@name='", name, "']")), "value")

terms <- data.frame(
  term_id = xml_attr(node_xml, "id"),
  group = as.integer(vapply(node_xml, attribute, character(1), "GROUP_ID")),
  hits = vapply(node_xml, attribute, character(1), "Hits"),
  SBPV = as.integer(vapply(node_xml, attribute, character(1), "_MEMBER_SBPV")),
  DBPV = as.integer(vapply(node_xml, attribute, character(1), "_MEMBER_DBPV"))
) |>
  filter(SBPV == 1L | DBPV == 1L) |>
  select(term_id, group, hits) |>
  distinct()

selected_proteins <- fread(
  "./BPV_ASCVD/data/bpv_protein_results.csv"
) |>
  filter(spec_id %in% c("SBPV_5y_VIM", "DBPV_5y_VIM")) |>
  group_by(feature) |>
  summarise(
    bpv_fdr = min(fdr),
    membership = if (n_distinct(pressure) == 2L) "Both" else paste0(first(pressure), " only"),
    .groups = "drop"
  )

protein_map <- readRDS(
  "./BPV_ASCVD/data/olink_gene_map.rds"
) |>
  select(feature, gene_symbol)

representatives <- terms |>
  separate_rows(hits, sep = "\\|") |>
  inner_join(protein_map, by = c("hits" = "gene_symbol"), relationship = "many-to-many") |>
  inner_join(selected_proteins, by = "feature") |>
  distinct(group, term_id, feature, hits, bpv_fdr, membership) |>
  count(group, feature, hits, bpv_fdr, membership, name = "terms_with_protein") |>
  arrange(group, desc(terms_with_protein), bpv_fdr, feature) |>
  group_by(group) |>
  slice_head(n = 1) |>
  ungroup() |>
  arrange(group) |>
  group_by(feature) |>
  summarise(
    belongs_to = paste0(sprintf("F%02d", group), collapse = ", "),
    membership = first(membership),
    .groups = "drop"
  )

fwrite(representatives, "./BPV_ASCVD/data/representative_proteins.csv")
