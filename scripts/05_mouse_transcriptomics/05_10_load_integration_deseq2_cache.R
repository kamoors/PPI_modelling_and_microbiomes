# SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0
# Load cached outlier-excluded RNA-seq objects

library(DESeq2)
library(dplyr)
base_dir <- project_path("rnaseq")
cache_file <- file.path(base_dir, "cache", "rnaseq_deseq2_cache_without_C02_G04.rds")
if (!project_file_exists(cache_file)) {
    stop("Cache file not found. Run 01b_cache_rnaseq_deseq2_objects_without_C02_G04.R first.")
}
rnaseq_cache_without_C02_G04 <- readRDS(project_input(cache_file))
meta <- rnaseq_cache_without_C02_G04$meta
organs <- rnaseq_cache_without_C02_G04$organs
txi <- rnaseq_cache_without_C02_G04$txi
dds_list <- rnaseq_cache_without_C02_G04$dds_list
res_list <- rnaseq_cache_without_C02_G04$res_list
res_all <- rnaseq_cache_without_C02_G04$res_all
counts_list <- rnaseq_cache_without_C02_G04$counts_list
meta_brain <- rnaseq_cache_without_C02_G04$meta_brain
meta_colon <- rnaseq_cache_without_C02_G04$meta_colon
meta_liver <- rnaseq_cache_without_C02_G04$meta_liver
exclude_samples <- rnaseq_cache_without_C02_G04$exclude_samples
summary_tbl <- rnaseq_cache_without_C02_G04$summary_tbl
message("Loaded cached RNA-seq objects without C02/G04:")
message("  dds_list")
message("  res_list")
message("  counts_list")
message("  meta / meta_brain / meta_colon / meta_liver")
message("  exclude_samples")
message("  summary_tbl")
message("Excluded samples:")
message(paste(exclude_samples, collapse = ", "))
