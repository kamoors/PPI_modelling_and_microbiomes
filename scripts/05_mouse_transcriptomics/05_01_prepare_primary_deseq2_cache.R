# SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0
# DESeq2 cache creation
# Run through run.R so input and output paths resolve within the package.
# Inputs: see data/input_manifest.csv and docs/RUN_ORDER.md.
# Writes generated files under results/; bundled inputs are read-only.

library(dplyr)
library(readxl)
library(DESeq2)
# tximport object is already generated; DESeq2 imports its matrices below.
library(readr)
base_dir <- ppi_path("rnaseq")
ppi_setwd(base_dir)
cache_dir <- file.path(base_dir, "cache")
dir.create(ppi_output(cache_dir, directory = TRUE), showWarnings = FALSE)
meta <- readxl::read_excel(ppi_input("Metadata.xlsx")) %>% dplyr::mutate(organ = sub(".*_", "", external_name), 
    donor = sub("(_m.*)", "", external_name), group = ifelse(grepl("^Ctrl", donor), "Ctrl", "PPI"))
txi <- readRDS(ppi_input("Salmon_TXI_object.rds"))
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
    res
})
names(res_list) <- organs
# Frozen annotations preserve study labels without a live database dependency.
gene_annotation <- readRDS(ppi_input(ppi_path("rnaseq/figure_coordinates_and_annotations.rds")))$gene_annotation
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
rnaseq_cache <- list(meta = meta, organs = organs, txi = txi, dds_list = dds_list, res_list = res_list, 
    res_all = res_all, counts_list = counts_list, meta_brain = meta_brain, meta_colon = meta_colon, meta_liver = meta_liver)
saveRDS(rnaseq_cache, ppi_output(file.path(cache_dir, "rnaseq_deseq2_cache.rds"), name = "primary_deseq2_cache.rds"))
write.csv(res_all, ppi_output(file.path(cache_dir, "DESeq2_results_by_organ.csv"), name = "de_seq2_results_by_organ.csv"), 
    row.names = FALSE)
write.csv(res_list$Liver, ppi_output(file.path(cache_dir, "DESeq2_results_liver.csv"), name = "de_seq2_results_liver.csv"), 
    row.names = FALSE)
write.csv(res_list$Colon, ppi_output(file.path(cache_dir, "DESeq2_results_colon.csv"), name = "de_seq2_results_colon.csv"), 
    row.names = FALSE)
write.csv(res_list$Brain, ppi_output(file.path(cache_dir, "DESeq2_results_brain.csv"), name = "de_seq2_results_brain.csv"), 
    row.names = FALSE)
message("RNA-seq DESeq2 cache saved.")

