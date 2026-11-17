# SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0
# Primary outlier-excluded DESeq2 contrasts; enrichment retained under scripts/optional/

source("scripts/helpers/setup.R")
invisible(run_analysis(
  section = "05_mouse_transcriptomics",
  directory = "rnaseq",
  script = "05_04_primary_deseq2_and_gsea.R",
  restricted = FALSE,
  code = {
library(DESeq2)
library(SummarizedExperiment)
library(dplyr)
library(tibble)
library(ggplot2)
load_rnaseq_script <- file.path(project_code("scripts/05_mouse_transcriptomics/05_02_load_primary_deseq2_cache.R"))
project_source(load_rnaseq_script)
out_dir <- file.path(getwd(), "RNAseq_without_C02_G04")
dir.create(project_output(out_dir, directory = TRUE), recursive = TRUE, showWarnings = FALSE)
exclude_samples <- c("25Jul1866-C02", "25Jul1866-G04")
run_deseq2_without_outliers <- function(dds, organ_name) {
    message("Running DESeq2 without outliers for: ", organ_name)
    keep_samples <- !colnames(dds) %in% exclude_samples
    dds_sub <- dds[, keep_samples]
    dds_sub$group <- droplevels(factor(dds_sub$group))
    dds_sub$group <- relevel(dds_sub$group, ref = "Ctrl")
    design(dds_sub) <- ~group
    keep_genes <- rowSums(DESeq2::counts(dds_sub) >= 10) >= ceiling(0.2 * ncol(dds_sub))
    dds_sub <- dds_sub[keep_genes, ]
    dds_sub <- DESeq2::DESeq(dds_sub)
    res <- DESeq2::results(dds_sub, contrast = c("group", "PPI", "Ctrl"))
    res_df <- as.data.frame(res) %>% tibble::rownames_to_column("gene") %>% dplyr::mutate(organ = organ_name,
        sig_0.1 = !is.na(padj) & padj < 0.1, sig_0.05 = !is.na(padj) & padj < 0.05) %>% dplyr::arrange(padj)
    list(dds = dds_sub, res = res, res_df = res_df)
}
deseq_no_outliers <- lapply(names(dds_list), function(org) {
    run_deseq2_without_outliers(dds = dds_list[[org]], organ_name = org)
})
names(deseq_no_outliers) <- names(dds_list)
res_no_outliers_list <- lapply(deseq_no_outliers, `[[`, "res_df")
for (org in names(res_no_outliers_list)) {
    write.csv(res_no_outliers_list[[org]], project_output(file.path(out_dir, paste0("DESeq2_", org, "_without_C02_G04.csv"))),
        row.names = FALSE)
}
res_no_outliers_all <- dplyr::bind_rows(res_no_outliers_list)
write.csv(res_no_outliers_all, project_output(file.path(out_dir, "DESeq2_all_organs_without_C02_G04.csv"),
    name = "de_seq2_all_organs_without_c02_g04.csv"), row.names = FALSE)
compare_to_original <- dplyr::bind_rows(lapply(names(deseq_no_outliers), function(org) {
    old <- as.data.frame(DESeq2::results(dds_list[[org]])) %>% tibble::rownames_to_column("gene") %>%
        dplyr::select(gene, log2FC_original = log2FoldChange, padj_original = padj)
    new <- deseq_no_outliers[[org]]$res_df %>% dplyr::select(gene, log2FC_without_outliers = log2FoldChange,
        padj_without_outliers = padj)
    old %>% dplyr::left_join(new, by = "gene") %>% dplyr::mutate(organ = org, sig_original = !is.na(padj_original) &
        padj_original < 0.1, sig_without_outliers = !is.na(padj_without_outliers) & padj_without_outliers <
        0.1, status = dplyr::case_when(sig_original & sig_without_outliers ~ "Significant_both", sig_original &
        !sig_without_outliers ~ "Only_original", !sig_original & sig_without_outliers ~ "Only_without_outliers",
        TRUE ~ "Not_significant"), delta_log2FC = log2FC_without_outliers - log2FC_original, abs_delta_log2FC = abs(delta_log2FC))
}))
write.csv(compare_to_original, project_output(file.path(out_dir, "DESeq2_original_vs_without_C02_G04.csv"),
    name = "de_seq2_original_vs_without_c02_g04.csv"), row.names = FALSE)
summary_tbl <- compare_to_original %>% dplyr::group_by(organ) %>% dplyr::summarise(n_sig_original = sum(sig_original,
    na.rm = TRUE), n_sig_without_outliers = sum(sig_without_outliers, na.rm = TRUE), overlap = sum(sig_original &
    sig_without_outliers, na.rm = TRUE), only_original = sum(sig_original & !sig_without_outliers, na.rm = TRUE),
    only_without_outliers = sum(!sig_original & sig_without_outliers, na.rm = TRUE), cor_log2FC = cor(log2FC_original,
        log2FC_without_outliers, use = "pairwise.complete.obs"), .groups = "drop")
write.csv(summary_tbl, project_output(file.path(out_dir, "DESeq2_summary_without_C02_G04.csv"), name = "de_seq2_summary_without_c02_g04.csv"),
    row.names = FALSE)
print(summary_tbl)
# Save fitted objects locally for optional enrichment and sensitivity analyses.
saveRDS(list(deseq_no_outliers = deseq_no_outliers,
    res_no_outliers_all = res_no_outliers_all, compare_to_original = compare_to_original,
    summary = summary_tbl), project_output(file.path(out_dir,
    "RNAseq_without_C02_G04_results.rds"), name = "primary_rna_without_outliers_results.rds"))
message("Primary DESeq2 fits complete. Optional enrichment is a separate script.")
}))
