# ==============================================================================
# METADATA
# ==============================================================================
# Author: Chaoyang Lin, MD
# Contact: chaoyanglint@163.com / Chaoyang Lin
#
# Dependencies: R >= 4.1.0
# Required packages: fst, ggplot2, ggcorrplot
# Required input: ./BPV_ASCVD/data/cohort_olink.fst
# Main output: ./BPV_ASCVD/figures/bpv_age_spearman.pdf
# ==============================================================================

library(fst)
library(ggplot2)
library(ggcorrplot)

variables <- c(
  "Age" = "Age_i0",
  "SBPV VIM (5 y)" = "vim_SBP_0to5y_i0",
  "SBPV SD (5 y)" = "sd_SBP_0to5y_i0",
  "SBPV CV (5 y)" = "cv_SBP_0to5y_i0",
  "SBPV ARV (5 y)" = "arv_SBP_0to5y_i0",
  "SBPV VIM (10 y)" = "vim_SBP_0to10y_i0",
  "DBPV VIM (5 y)" = "vim_DBP_0to5y_i0",
  "DBPV SD (5 y)" = "sd_DBP_0to5y_i0",
  "DBPV CV (5 y)" = "cv_DBP_0to5y_i0",
  "DBPV ARV (5 y)" = "arv_DBP_0to5y_i0",
  "DBPV VIM (10 y)" = "vim_DBP_0to10y_i0"
)

cohort <- fst::read_fst(
  "./BPV_ASCVD/data/cohort_olink.fst",
  columns = unname(variables)
)
names(cohort) <- names(variables)

# Pairwise complete observations retain available information for each pair.
p <- ggcorrplot(
  cor(cohort, method = "spearman", use = "pairwise.complete.obs"),
  hc.order = FALSE, type = "lower", show.diag = TRUE,
  lab = TRUE, digits = 3, nsmall = 3, lab_size = 3.2,
  colors = c("#6698B6", "white", "#D48377"), outline.color = "white",
  legend.title = "Spearman rho", tl.cex = 10, tl.col = "black"
) +
  theme(
    text = element_text(family = "sans"),
    panel.grid = element_blank(),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9),
    plot.margin = margin(14, 14, 14, 14)
  )

ggsave(
  "./BPV_ASCVD/figures/bpv_age_spearman.pdf",
  p, width = 11, height = 10.5, device = cairo_pdf
)
