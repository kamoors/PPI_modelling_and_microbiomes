# SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0
# Sample diagnostics

source("scripts/helpers/setup.R")
invisible(run_analysis(
  section = "05_mouse_transcriptomics",
  directory = "rnaseq",
  script = "05_03_check_rna_samples.R",
  restricted = FALSE,
  code = {
library(DESeq2)
library(dplyr)
library(tibble)
project_source(project_code("scripts/05_mouse_transcriptomics/05_02_load_primary_deseq2_cache.R"))
run_outlier_diagnostics <- function(dds, organ_name, n_pcs_use = 5) {
    message("Running outlier diagnostics for: ", organ_name)
    vsd <- DESeq2::vst(dds, blind = TRUE)
    expr <- SummarizedExperiment::assay(vsd)
    sample_names <- colnames(expr)
    meta_df <- as.data.frame(SummarizedExperiment::colData(dds))
    cors <- stats::cor(expr, use = "pairwise.complete.obs", method = "pearson")
    median_cor <- apply(cors, 1, median, na.rm = TRUE)
    ctrl_samples <- rownames(subset(meta_df, group == "Ctrl"))
    ppi_samples <- rownames(subset(meta_df, group == "PPI"))
    if (length(ctrl_samples) < 2) {
        warning("Fewer than 2 Ctrl samples in ", organ_name)
    }
    if (length(ppi_samples) < 2) {
        warning("Fewer than 2 PPI samples in ", organ_name)
    }
    ctrl_centroid <- rowMeans(expr[, ctrl_samples, drop = FALSE], na.rm = TRUE)
    ppi_centroid <- rowMeans(expr[, ppi_samples, drop = FALSE], na.rm = TRUE)
    centroid_cor <- apply(expr, 2, function(x) {
        c(cor_to_ctrl = stats::cor(x, ctrl_centroid, use = "pairwise.complete.obs"), cor_to_ppi = stats::cor(x,
            ppi_centroid, use = "pairwise.complete.obs"))
    })
    centroid_cor <- t(centroid_cor)
    cooks <- SummarizedExperiment::assays(dds)[["cooks"]]
    if (is.null(cooks)) {
        message("  Cook's distances not found; running DESeq() for ", organ_name)
        dds_tmp <- DESeq2::DESeq(dds, quiet = TRUE)
        cooks <- SummarizedExperiment::assays(dds_tmp)[["cooks"]]
    }
    median_cooks <- apply(cooks, 2, median, na.rm = TRUE)
    pca <- stats::prcomp(t(expr), scale. = FALSE)
    n_pcs <- min(n_pcs_use, ncol(pca$x))
    pcs <- pca$x[, seq_len(n_pcs), drop = FALSE]
    md <- stats::mahalanobis(pcs, center = colMeans(pcs), cov = stats::cov(pcs))
    md_p <- stats::pchisq(md, df = ncol(pcs), lower.tail = FALSE)
    out <- data.frame(sample = sample_names, organ = organ_name, median_correlation = median_cor[sample_names],
        cor_to_ctrl = centroid_cor[sample_names, "cor_to_ctrl"], cor_to_ppi = centroid_cor[sample_names,
            "cor_to_ppi"], median_cooks = median_cooks[sample_names], mahalanobis = md[sample_names],
        mahalanobis_p = md_p[sample_names], stringsAsFactors = FALSE)
    out <- out %>% dplyr::mutate(rank_corr = rank(median_correlation, ties.method = "average"), rank_cooks = rank(-median_cooks,
        ties.method = "average"), rank_mahal = rank(-mahalanobis, ties.method = "average"), outlier_score = rank_corr +
        rank_cooks + rank_mahal) %>% dplyr::arrange(outlier_score)
    return(out)
}
outlier_results <- dplyr::bind_rows(lapply(names(dds_list), function(org) {
    run_outlier_diagnostics(dds = dds_list[[org]], organ_name = org, n_pcs_use = 5)
}))
top_outliers_by_organ <- outlier_results %>% dplyr::group_by(organ) %>% dplyr::slice_min(outlier_score,
    n = 5, with_ties = FALSE) %>% dplyr::ungroup() %>% dplyr::arrange(organ, outlier_score)
print(top_outliers_by_organ)
outlier_results_with_fdr <- outlier_results %>% dplyr::mutate(mahalanobis_fdr = p.adjust(mahalanobis_p,
    method = "BH")) %>% dplyr::arrange(mahalanobis_p)
print(outlier_results_with_fdr)
write.csv(outlier_results, project_output("RNAseq_outlier_diagnostics_all_organs.csv", name = "rn_aseq_outlier_diagnostics_all_organs.csv"),
    row.names = FALSE)
write.csv(top_outliers_by_organ, project_output("RNAseq_top_outlier_candidates_by_organ.csv", name = "rn_aseq_top_outlier_candidates_by_organ.csv"),
    row.names = FALSE)
write.csv(outlier_results_with_fdr, project_output("RNAseq_outlier_diagnostics_all_organs_with_mahalanobis_fdr.csv",
    name = "rn_aseq_outlier_diagnostics_all_organs_with_mahalanobis_fdr.csv"), row.names = FALSE)
}))
