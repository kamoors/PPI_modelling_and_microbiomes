# SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0
# Corrected liver GO export
# Run through run.R so input and output paths resolve within the package.
# Inputs: see data/input_manifest.csv and docs/RUN_ORDER.md.
# Writes generated files under results/; bundled inputs are read-only.

required <- c("readxl", "dplyr", "tibble", "clusterProfiler", "org.Mm.eg.db", "AnnotationDbi")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0) {
    stop("Missing packages: ", paste(missing, collapse = ", "), "\nInstall Bioconductor packages with:\n", 
        "BiocManager::install(c('clusterProfiler','org.Mm.eg.db','AnnotationDbi'))")
}
library(readxl)
library(dplyr)
library(tibble)
library(clusterProfiler)
library(org.Mm.eg.db)
library(AnnotationDbi)
project_root <- ppi_path("")
# These primary and sensitivity result tables are analysis outputs, not raw input.
# Refitting the DESeq2 models is available through --workflow rna-primary.
tables <- readRDS(ppi_input(ppi_path("reference/supplementary_data_3_tables.rds")))
out_dir <- ppi_path("rnaseq/Results/GO_Liver_DESeq2_vs_limma")
deseq_sheet <- grep("DESeq2_Primary$", names(tables), ignore.case=TRUE, value=TRUE)
limma_sheet <- grep("limma_Primary$", names(tables), ignore.case=TRUE, value=TRUE)
stopifnot(length(deseq_sheet)==1L, length(limma_sheet)==1L)
deseq_all <- tables[[deseq_sheet]]
limma_all <- tables[[limma_sheet]]
required_cols <- c("gene", "organ", "log2FoldChange", "padj")
if (!all(required_cols %in% names(deseq_all))) {
    stop("DESeq2 table is missing required columns: ", paste(setdiff(required_cols, names(deseq_all)), 
        collapse = ", "))
}
if (!all(required_cols %in% names(limma_all))) {
    stop("limma table is missing required columns: ", paste(setdiff(required_cols, names(limma_all)), 
        collapse = ", "))
}
deseq_liver <- deseq_all %>% filter(organ == "Liver")
limma_liver <- limma_all %>% filter(organ == "Liver")
message("DESeq2 liver genes: ", nrow(deseq_liver))
message("limma liver genes: ", nrow(limma_liver))
strip_ensembl_version <- function(x) {
    sub("[.]\\d+$", "", as.character(x))
}
map_to_entrez <- function(ids) {
    ids <- unique(strip_ensembl_version(ids))
    ids <- ids[!is.na(ids) & ids != ""]
    AnnotationDbi::select(org.Mm.eg.db, keys = ids, keytype = "ENSEMBL", columns = c("ENTREZID", "SYMBOL")) %>% 
        filter(!is.na(ENTREZID)) %>% distinct(ENSEMBL, .keep_all = TRUE)
}
run_one_go <- function(df, method_name, direction, alpha = 0.1) {
    stopifnot(direction %in% c("Up", "Down"))
    if ("pvalue" %in% names(df)) {
        bg_ids <- df %>% filter(!is.na(pvalue)) %>% pull(gene)
    }
    else {
        bg_ids <- df$gene
    }
    sig <- df %>% filter(!is.na(padj), padj < alpha, !is.na(log2FoldChange), if (direction == "Up") 
        log2FoldChange > 0
    else log2FoldChange < 0)
    message(method_name, " ", direction, ": ", nrow(sig), " DEGs before ID mapping")
    sig_map <- map_to_entrez(sig$gene)
    bg_map <- map_to_entrez(bg_ids)
    sig_entrez <- unique(sig_map$ENTREZID)
    bg_entrez <- unique(bg_map$ENTREZID)
    message(method_name, " ", direction, ": ", length(sig_entrez), " mapped DEGs; ", length(bg_entrez), 
        " mapped background genes")
    if (length(sig_entrez) < 5) {
        return(tibble(method = method_name, direction = direction, ID = NA_character_, Description = NA_character_, 
            GeneRatio = NA_character_, BgRatio = NA_character_, pvalue = NA_real_, p.adjust = NA_real_, 
            qvalue = NA_real_, geneID = NA_character_, Count = NA_integer_, n_input_DEGs = nrow(sig), 
            n_mapped_DEGs = length(sig_entrez), n_background = length(bg_entrez), note = "Too few mapped DEGs for GO ORA"))
    }
    ego <- clusterProfiler::enrichGO(gene = sig_entrez, universe = bg_entrez, OrgDb = org.Mm.eg.db, keyType = "ENTREZID", 
        ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1, readable = TRUE)
    out <- as.data.frame(ego)
    if (nrow(out) == 0) {
        return(tibble(method = method_name, direction = direction, ID = NA_character_, Description = NA_character_, 
            GeneRatio = NA_character_, BgRatio = NA_character_, pvalue = NA_real_, p.adjust = NA_real_, 
            qvalue = NA_real_, geneID = NA_character_, Count = NA_integer_, n_input_DEGs = nrow(sig), 
            n_mapped_DEGs = length(sig_entrez), n_background = length(bg_entrez), note = "No GO BP terms returned"))
    }
    out %>% mutate(method = method_name, direction = direction, n_input_DEGs = nrow(sig), n_mapped_DEGs = length(sig_entrez), 
        n_background = length(bg_entrez), note = NA_character_, .before = 1)
}
go_all <- bind_rows(run_one_go(deseq_liver, "DESeq2", "Up"), run_one_go(deseq_liver, "DESeq2", "Down"), 
    run_one_go(limma_liver, "limma_duplicateCorrelation", "Up"), run_one_go(limma_liver, "limma_duplicateCorrelation", 
        "Down")) %>% arrange(method, direction, p.adjust, pvalue)
full_csv <- file.path(out_dir, "GO_BP_liver_DESeq2_vs_limma_all_results.csv")
sig_csv <- file.path(out_dir, "GO_BP_liver_DESeq2_vs_limma_FDR0p1.csv")
write.csv(go_all, ppi_output(full_csv), row.names = FALSE)
write.csv(go_all %>% filter(!is.na(p.adjust), p.adjust < 0.1), ppi_output(sig_csv), row.names = FALSE)
message("\nSaved:\n", full_csv)
message("\nSaved significant subset:\n", sig_csv)
summary_tbl <- go_all %>% dplyr::group_by(method, direction) %>% dplyr::summarise(n_input_DEGs = .data$n_input_DEGs[1], 
    n_mapped_DEGs = .data$n_mapped_DEGs[1], n_background = .data$n_background[1], n_GO_terms = sum(!is.na(.data$ID)), 
    n_GO_FDR0p1 = sum(!is.na(.data$p.adjust) & .data$p.adjust < 0.1), .groups = "drop")
print(summary_tbl)
write.csv(summary_tbl, ppi_output(file.path(out_dir, "GO_BP_liver_DESeq2_vs_limma_summary.csv"), name = "go_bp_liver_de_seq2_vs_limma_summary.csv"), 
    row.names = FALSE)
