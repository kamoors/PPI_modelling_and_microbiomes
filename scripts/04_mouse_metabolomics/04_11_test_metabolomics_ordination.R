# SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0
# PCA/Procrustes statistics

# Primary manuscript test only; broader exploratory statistics are optional.

source("scripts/helpers/setup.R")
invisible(run_analysis(
  section = "04_mouse_metabolomics",
  directory = "metabolomics",
  script = "04_11_test_metabolomics_ordination.R",
  restricted = FALSE,
  code = {
suppressPackageStartupMessages({
    library(dplyr)
    library(tibble)
    library(readr)
    library(ggplot2)
    library(vegan)
    library(permute)
})
set.seed(20260714)
procrustes_permutations <- 9999L
bundle <- readRDS(project_input(project_path("metabolomics/Metabolomics_PCA/PPI_metabolomics_PCA_bundle.rds")))
pca_objects <- bundle$pca_objects
analysis_out_dir <- project_path("metabolomics/Metabolomics_PCA/PCA_statistical_analysis")
run_exvivo_procrustes <- function(nmr_obj, hilic_obj) {
    nmr_scores <- nmr_obj$scores %>% select(Match_ID, Donor, Group, Day_num, starts_with("PC"))
    hilic_scores <- hilic_obj$scores %>% select(Match_ID, starts_with("PC"))
    matched <- inner_join(nmr_scores, hilic_scores, by = "Match_ID", suffix = c("_NMR", "_HILIC")) %>%
        arrange(Match_ID)
    if (nrow(matched) < 8L) {
        return(list(summary = tibble(note = "Insufficient matched samples")))
    }
    X <- as.matrix(matched[, c("PC1_NMR", "PC2_NMR")])
    Y <- as.matrix(matched[, c("PC1_HILIC", "PC2_HILIC")])
    proc <- vegan::procrustes(X, Y, symmetric = TRUE)
    ctrl <- permute::how(nperm = procrustes_permutations)
    permute::setBlocks(ctrl) <- factor(matched$Donor)
    test <- vegan::protest(X, Y, permutations = ctrl)
    rotated <- as.data.frame(proc$Yrot) %>% setNames(c("HILIC_rotated_1", "HILIC_rotated_2")) %>% bind_cols(matched %>%
        select(Match_ID, Donor, Group, Day_num), .)
    nmr_plot <- as.data.frame(proc$X) %>% setNames(c("NMR_1", "NMR_2")) %>% bind_cols(matched %>% select(Match_ID),
        .)
    plot_data <- left_join(rotated, nmr_plot, by = "Match_ID")
    summary <- tibble(n_matched_samples = nrow(matched), procrustes_correlation = as.numeric(test$t0),
        permutation_p = as.numeric(test$signif), permutations = procrustes_permutations, restriction = "Permutations restricted within donor")
    p <- ggplot(plot_data, aes(colour = Group)) + geom_segment(aes(x = NMR_1, y = NMR_2, xend = HILIC_rotated_1,
        yend = HILIC_rotated_2), linewidth = 0.45, alpha = 0.55, arrow = grid::arrow(length = grid::unit(0.06,
        "inches"))) + geom_point(aes(NMR_1, NMR_2), shape = 16, size = 2.6) + geom_point(aes(HILIC_rotated_1,
        HILIC_rotated_2), shape = 1, size = 2.6) + scale_colour_manual(values = c(ctrl = "#2878B5", ppi = "#C73E3A")) +
        labs(title = "Concordance of baseline-centred ex vivo PCA configurations", subtitle = paste0("Filled: NMR; open: rotated HILIC; Procrustes r = ",
            formatC(test$t0, digits = 3, format = "f"), ", permutation P = ", format.pval(test$signif,
                digits = 3)), x = "Procrustes axis 1", y = "Procrustes axis 2", colour = "Inoculum") +
        theme_classic(base_size = 11)
    list(summary = summary, plot_data = plot_data, plot = p)
}
stopifnot(all(c("exvivo_nmr_delta", "exvivo_hilic_delta") %in% names(pca_objects)))
procrustes_result <- run_exvivo_procrustes(pca_objects$exvivo_nmr_delta, pca_objects$exvivo_hilic_delta)
write_csv(procrustes_result$summary, project_output(file.path(analysis_out_dir,
    "PCA_exvivo_NMR_HILIC_procrustes.csv"), name = "ex_vivo_nmr_hilic_procrustes_statistics.csv"))
write_csv(procrustes_result$plot_data, project_output(file.path(analysis_out_dir,
    "Procrustes_coordinates.csv"), name = "ex_vivo_nmr_hilic_procrustes_coordinates.csv"))
ggsave(project_output(file.path(analysis_out_dir, "Exvivo_NMR_HILIC_Procrustes.png"),
    name = "ex_vivo_nmr_hilic_procrustes.png"), procrustes_result$plot,
    width = 7, height = 5.8, dpi = 450, bg = "white")
print(procrustes_result$summary)
}))
