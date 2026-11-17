# SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0
# Create outlier-excluded RNA cache

source("scripts/helpers/setup.R")
invisible(run_analysis(
  section = "05_mouse_transcriptomics",
  directory = "rnaseq/cache",
  script = "05_09_prepare_integration_deseq2_cache.R",
  restricted = FALSE,
  code = {
library(dplyr)
library(readxl)
library(DESeq2)
# tximport object is already generated; DESeq2 imports its matrices below.
library(readr)
base_dir <- project_path("rnaseq")
project_setwd(base_dir)
cache_dir <- file.path(base_dir, "cache")
dir.create(project_output(cache_dir, directory = TRUE), showWarnings = FALSE)
exclude_samples <- c("25Jul1866-C02", "25Jul1866-G04")
meta <- readxl::read_excel(project_input("Metadata.xlsx")) %>% dplyr::mutate(organ = sub(".*_", "", external_name),
    donor = sub("(_m.*)", "", external_name), group = ifelse(grepl("^Ctrl", donor), "Ctrl", "PPI")) %>%
    dplyr::filter(!name %in% exclude_samples)
txi <- readRDS(project_input("Salmon_TXI_object.rds"))
clean_names <- sub(".*(25Jul1866-[A-H][0-9]{2}).*", "\\1", colnames(txi$counts))
colnames(txi$counts) <- clean_names
colnames(txi$abundance) <- clean_names
colnames(txi$length) <- clean_names
txi$counts <- txi$counts[, match(meta$name, colnames(txi$counts))]
txi$abundance <- txi$abundance[, match(meta$name, colnames(txi$abundance))]
txi$length <- txi$length[, match(meta$name, colnames(txi$length))]
stopifnot(all(colnames(txi$counts) == meta$name))
organs <- unique(meta$organ)
prevalence_filter <- function(dds, min_count = 10, min_samps = ceiling(0.2 * ncol(dds))) {
    keep <- rowSums(DESeq2::counts(dds) >= min_count) >= min_samps
    dds[keep, ]
}
dds_list <- lapply(organs, function(org) {
    samples <- meta %>% filter(organ == org)
    idx <- match(samples$name, colnames(txi$counts))
    txi_sub <- lapply(txi, function(x) {
        if (is.matrix(x)) {
            x[, idx, drop = FALSE]
        }
        else {
            x
        }
    })
    dds <- DESeq2::DESeqDataSetFromTximport(txi = txi_sub, colData = samples, design = ~group)
    dds$group <- droplevels(factor(dds$group))
    dds$group <- relevel(dds$group, ref = "Ctrl")
    dds <- prevalence_filter(dds)
    dds <- DESeq2::DESeq(dds)
    dds
})
names(dds_list) <- organs
res_list <- lapply(names(dds_list), function(org) {
    res <- DESeq2::results(dds_list[[org]], contrast = c("group", "PPI", "Ctrl"))
    res <- as.data.frame(res)
    res$gene <- rownames(res)
    res$organ <- org
    res <- res %>% dplyr::mutate(sig_0.1 = !is.na(padj) & padj < 0.1, sig_0.05 = !is.na(padj) & padj <
        0.05)
    res
})
names(res_list) <- organs
gene_annotation <- readRDS(project_input(project_path("rnaseq/figure_coordinates_and_annotations.rds")))$gene_annotation
map_ens_sym <- function(res) {
    ind <- match(res$gene, gene_annotation$gene)
    res$geneSymbol <- gene_annotation$geneSymbol[ind]
    res$geneName <- gene_annotation$geneName[ind]
    res
}
res_list <- lapply(res_list, function(x) {
    map_ens_sym(x) %>% dplyr::arrange(padj)
})
res_all <- dplyr::bind_rows(res_list)
meta_brain <- as.data.frame(SummarizedExperiment::colData(dds_list$Brain))
meta_colon <- as.data.frame(SummarizedExperiment::colData(dds_list$Colon))
meta_liver <- as.data.frame(SummarizedExperiment::colData(dds_list$Liver))
counts_list <- lapply(dds_list, function(dds) {
    DESeq2::counts(dds, normalized = FALSE)
})
summary_tbl <- res_all %>% dplyr::group_by(organ) %>% dplyr::summarise(n_sig_0.1 = sum(sig_0.1, na.rm = TRUE),
    n_sig_0.05 = sum(sig_0.05, na.rm = TRUE), .groups = "drop")
print(summary_tbl)
rnaseq_cache_without_C02_G04 <- list(meta = meta, organs = organs, txi = txi, dds_list = dds_list, res_list = res_list,
    res_all = res_all, counts_list = counts_list, meta_brain = meta_brain, meta_colon = meta_colon, meta_liver = meta_liver,
    exclude_samples = exclude_samples, summary_tbl = summary_tbl)
saveRDS(rnaseq_cache_without_C02_G04, project_output(file.path(cache_dir, "rnaseq_deseq2_cache_without_C02_G04.rds"),
    name = "integration_deseq2_cache.rds"))
write.csv(res_all, project_output(file.path(cache_dir, "DESeq2_results_by_organ_without_C02_G04.csv"), name = "de_seq2_results_by_organ_without_c02_g04.csv"),
    row.names = FALSE)
write.csv(res_list$Liver, project_output(file.path(cache_dir, "DESeq2_results_liver_without_C02_G04.csv"),
    name = "de_seq2_results_liver_without_c02_g04.csv"), row.names = FALSE)
write.csv(res_list$Colon, project_output(file.path(cache_dir, "DESeq2_results_colon_without_C02_G04.csv"),
    name = "de_seq2_results_colon_without_c02_g04.csv"), row.names = FALSE)
write.csv(res_list$Brain, project_output(file.path(cache_dir, "DESeq2_results_brain_without_C02_G04.csv"),
    name = "de_seq2_results_brain_without_c02_g04.csv"), row.names = FALSE)
write.csv(summary_tbl, project_output(file.path(cache_dir, "DESeq2_summary_without_C02_G04.csv"), name = "de_seq2_summary_without_c02_g04.csv"),
    row.names = FALSE)
message("RNA-seq DESeq2 cache without C02/G04 saved.")
}))
