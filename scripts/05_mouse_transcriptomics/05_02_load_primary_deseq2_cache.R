# SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0
# Load cached original RNA-seq objects

library(DESeq2)
library(dplyr)
base_dir <- project_path("rnaseq")
cache_file <- file.path(base_dir, "cache", "rnaseq_deseq2_cache.rds")
if (!project_file_exists(cache_file)) {
    stop("Cache file not found. Run 05_01_prepare_primary_deseq2_cache.R first.")
}
rnaseq_cache <- readRDS(project_input(cache_file))
meta <- rnaseq_cache$meta
organs <- rnaseq_cache$organs
txi <- rnaseq_cache$txi
dds_list <- rnaseq_cache$dds_list
res_list <- rnaseq_cache$res_list
res_all <- rnaseq_cache$res_all
counts_list <- rnaseq_cache$counts_list
meta_brain <- rnaseq_cache$meta_brain
meta_colon <- rnaseq_cache$meta_colon
meta_liver <- rnaseq_cache$meta_liver
message("Loaded cached RNA-seq objects:")
message("  dds_list")
message("  res_list")
message("  counts_list")
message("  meta / meta_brain / meta_colon / meta_liver")
