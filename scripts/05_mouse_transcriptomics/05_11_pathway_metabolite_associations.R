# SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0
# Create pathway scores read by cross-layer preparation

source("scripts/helpers/setup.R")
invisible(run_analysis(
  section = "05_mouse_transcriptomics",
  directory = "rnaseq",
  script = "05_11_pathway_metabolite_associations.R",
  restricted = FALSE,
  code = {
library(dplyr)
library(tidyr)
library(purrr)
library(stringr)
library(GSVA)
library(msigdbr)
library(AnnotationDbi)
library(org.Mm.eg.db)
library(SummarizedExperiment)
library(writexl)
load_rnaseq_script <- file.path(project_code("scripts/05_mouse_transcriptomics/05_10_load_integration_deseq2_cache.R"))
load_metab_script <- file.path(project_code("scripts/00_shared_helpers/00_05_load_mouse_metabolomics.R"))
make_bridge_script <- file.path(project_code("scripts/00_shared_helpers/00_06_match_metabolomics_rna_samples.R"))
project_source(load_rnaseq_script)
project_source(load_metab_script)
project_source(make_bridge_script)
excluded_samples <- c("25Jul1866-C02", "25Jul1866-G04")
stopifnot(!any(excluded_samples %in% meta$name))
stopifnot(!any(excluded_samples %in% unlist(lapply(dds_list, colnames))))
message("Confirmed excluded samples are absent from RNA-seq objects.")
out_dir <- project_path("rnaseq/Results/hallmark_reactome_metabolite_correlations_without_C02_G04")
dir.create(project_output(out_dir, directory = TRUE), recursive = TRUE, showWarnings = FALSE)
min_n <- 6
fdr_cutoff <- 0.1
# Read the same saved result tables without requiring a fitted-object bundle.
gsea_no_outliers <- setNames(lapply(c("Liver", "Colon", "Brain"), function(org)
    setNames(lapply(c("HALLMARK", "REACTOME"), function(collection)
        read.csv(project_input(project_path(paste0("rnaseq/Results/RNAseq_without_C02_G04/GSEA_",
            org, "_", collection, "_without_C02_G04.csv"))), check.names=FALSE)),
        c("HALLMARK", "REACTOME"))), c("Liver", "Colon", "Brain"))
cp_gsea_sig <- purrr::imap_dfr(gsea_no_outliers, function(org_res, org) {
    dplyr::bind_rows(org_res$HALLMARK %>% dplyr::mutate(organ = org, library = "HALLMARK", ID = pathway,
        Description = pathway), org_res$REACTOME %>% dplyr::mutate(organ = org, library = "REACTOME",
        ID = pathway, Description = pathway))
}) %>% dplyr::filter(!is.na(padj), padj < fdr_cutoff) %>% dplyr::mutate(p.adjust = padj)
message("Significant no-outlier Hallmark/Reactome pathways:")
print(table(cp_gsea_sig$organ, cp_gsea_sig$library))
msig_h <- msigdbr::msigdbr(species = "Mus musculus", collection = "H")
msig_reactome <- msigdbr::msigdbr(species = "Mus musculus", collection = "C2", subcollection = "CP:REACTOME")
pathways_hr <- dplyr::bind_rows(msig_h %>% dplyr::select(gs_name, gene_symbol) %>% dplyr::mutate(library = "HALLMARK"),
    msig_reactome %>% dplyr::select(gs_name, gene_symbol) %>% dplyr::mutate(library = "REACTOME"))
pathway_sets_hr <- split(pathways_hr$gene_symbol, pathways_hr$gs_name)
pathway_sets_hr <- lapply(pathway_sets_hr, unique)
pathway_scores_for_organ_hr <- function(dds, organ_name, pathway_sets) {
    message("Calculating ssGSEA scores for ", organ_name)
    vsd <- DESeq2::vst(dds, blind = FALSE)
    expr <- SummarizedExperiment::assay(vsd)
    rownames(expr) <- sub("\\.\\d+$", "", rownames(expr))
    gene_map <- AnnotationDbi::select(org.Mm.eg.db, keys = rownames(expr), keytype = "ENSEMBL", columns = "SYMBOL") %>%
        dplyr::filter(!is.na(SYMBOL), SYMBOL != "") %>% dplyr::mutate(ENSEMBL = sub("\\.\\d+$", "", ENSEMBL)) %>%
        dplyr::distinct(ENSEMBL, .keep_all = TRUE)
    symbols <- gene_map$SYMBOL[match(rownames(expr), gene_map$ENSEMBL)]
    keep <- !is.na(symbols) & symbols != ""
    expr <- expr[keep, , drop = FALSE]
    rownames(expr) <- symbols[keep]
    expr <- rowsum(expr, group = rownames(expr))/as.vector(table(rownames(expr)))
    pathway_sets_use <- lapply(pathway_sets, function(g) intersect(g, rownames(expr)))
    pathway_sets_use <- pathway_sets_use[lengths(pathway_sets_use) >= 10 & lengths(pathway_sets_use) <=
        500]
    param <- GSVA::ssgseaParam(exprData = expr, geneSets = pathway_sets_use)
    scores <- GSVA::gsva(param)
    scores_df <- as.data.frame(t(scores), check.names = FALSE)
    scores_df$sample <- rownames(scores_df)
    scores_df$organ <- organ_name
    scores_df
}
pathway_scores_hr_list <- lapply(names(dds_list), function(org) {
    pathway_scores_for_organ_hr(dds = dds_list[[org]], organ_name = org, pathway_sets = pathway_sets_hr)
})
names(pathway_scores_hr_list) <- names(dds_list)
pathway_scores_hr_all <- dplyr::bind_rows(pathway_scores_hr_list)
meta_clean <- meta %>% dplyr::select(sample = name, group, donor)
pathway_scores_hr_all <- pathway_scores_hr_all %>% dplyr::left_join(meta_clean, by = "sample")
write.csv(pathway_scores_hr_all, project_output(file.path(out_dir, "Hallmark_Reactome_ssGSEA_scores_all_organs.csv"),
    name = "hallmark_reactome_ss_gsea_scores_all_organs.csv"), row.names = FALSE)
sig_hr_pathways_by_org <- cp_gsea_sig %>% dplyr::filter(library %in% c("HALLMARK", "REACTOME"), p.adjust <
    fdr_cutoff) %>% dplyr::mutate(pathway = ID) %>% dplyr::group_by(organ) %>% dplyr::summarise(pathways = list(unique(pathway)),
    .groups = "drop") %>% tibble::deframe()
matched_hr <- lapply(names(sig_hr_pathways_by_org), function(org) {
    data.frame(organ = org, pathway = sig_hr_pathways_by_org[[org]], present_in_scores = sig_hr_pathways_by_org[[org]] %in%
        colnames(pathway_scores_hr_all))
}) %>% dplyr::bind_rows()
write.csv(matched_hr, project_output(file.path(out_dir, "Hallmark_Reactome_pathway_matching_check.csv"),
    name = "hallmark_reactome_pathway_matching_check.csv"), row.names = FALSE)
print(table(matched_hr$organ, matched_hr$present_in_scores))
if (any(!matched_hr$present_in_scores)) {
    warning("Some significant pathways are not present in the ssGSEA score matrix. ", "Check pathway naming before interpreting correlation results.")
}
safe_cor <- function(x, y, method = "spearman", min_n = 6) {
    keep <- is.finite(x) & is.finite(y)
    if (sum(keep) < min_n) {
        return(c(rho = NA_real_, pvalue = NA_real_, n = sum(keep)))
    }
    test <- suppressWarnings(cor.test(x[keep], y[keep], method = method, exact = FALSE))
    c(rho = unname(test$estimate), pvalue = test$p.value, n = sum(keep))
}
make_gut_M_block <- function(dds, compartment, bridges, mice_combined) {
    organ_name <- as.character(unique(SummarizedExperiment::colData(dds)$organ))[1]
    br <- bridges[[organ_name]][[compartment]] %>% dplyr::filter(matched) %>% dplyr::select(Sample_Name,
        rna_sample)
    meta_gut <- pos_mice$metadata
    X <- as.data.frame(mice_combined, check.names = FALSE)
    X$Sample_Name <- meta_gut$Sample_Name
    M <- X %>% dplyr::inner_join(br, by = "Sample_Name") %>% dplyr::relocate(rna_sample) %>% dplyr::select(-Sample_Name) %>%
        dplyr::group_by(rna_sample) %>% dplyr::summarise(dplyr::across(where(is.numeric), ~mean(.x, na.rm = TRUE)),
        .groups = "drop")
    rnaseq_order <- colnames(dds)
    M2 <- matrix(NA_real_, nrow = length(rnaseq_order), ncol = ncol(M) - 1, dimnames = list(rnaseq_order,
        colnames(M)[-1]))
    M2[match(M$rna_sample, rnaseq_order), ] <- as.matrix(M[, -1, drop = FALSE])
    M2
}
make_serum_M_block <- function(dds, serum_combined, serum_meta, serum_bridges) {
    organ_name <- as.character(unique(SummarizedExperiment::colData(dds)$organ))[1]
    br <- serum_bridges[[organ_name]] %>% dplyr::filter(matched) %>% dplyr::select(Sample_Name, rna_sample)
    X <- as.data.frame(serum_combined, check.names = FALSE)
    X$Sample_Name <- serum_meta$Sample_Name
    M <- X %>% dplyr::inner_join(br, by = "Sample_Name") %>% dplyr::relocate(rna_sample) %>% dplyr::select(-Sample_Name) %>%
        dplyr::group_by(rna_sample) %>% dplyr::summarise(dplyr::across(where(is.numeric), ~mean(.x, na.rm = TRUE)),
        .groups = "drop")
    rnaseq_order <- colnames(dds)
    M2 <- matrix(NA_real_, nrow = length(rnaseq_order), ncol = ncol(M) - 1, dimnames = list(rnaseq_order,
        colnames(M)[-1]))
    M2[match(M$rna_sample, rnaseq_order), ] <- as.matrix(M[, -1, drop = FALSE])
    M2
}
correlate_gut_hr_pathways_metabolites <- function(organ_name, compartment) {
    dds <- dds_list[[organ_name]]
    M <- make_gut_M_block(dds = dds, compartment = compartment, bridges = bridges, mice_combined = mice_combined)
    score_df <- pathway_scores_hr_all %>% dplyr::filter(organ == organ_name) %>% dplyr::select(sample,
        group, dplyr::any_of(sig_hr_pathways_by_org[[organ_name]]))
    pathway_cols <- intersect(sig_hr_pathways_by_org[[organ_name]], colnames(score_df))
    common <- intersect(score_df$sample, rownames(M))
    if (length(pathway_cols) == 0 || length(common) < min_n) {
        return(tibble::tibble())
    }
    score_df <- score_df %>% dplyr::filter(sample %in% common) %>% dplyr::arrange(match(sample, common))
    M <- M[score_df$sample, , drop = FALSE]
    out <- expand.grid(pathway = pathway_cols, metabolite = colnames(M), stringsAsFactors = FALSE)
    purrr::pmap_dfr(out, function(pathway, metabolite) {
        test <- safe_cor(x = score_df[[pathway]], y = M[, metabolite], min_n = min_n)
        tibble::tibble(organ = organ_name, compartment = compartment, source = "gut", pathway = pathway,
            metabolite = metabolite, rho = as.numeric(test["rho"]), pvalue = as.numeric(test["pvalue"]),
            n = as.numeric(test["n"]))
    })
}
gut_hr_cor_results <- list()
for (org in names(dds_list)) {
    for (comp in c("cecum", "feces")) {
        message("Correlating gut: ", org, " / ", comp)
        gut_hr_cor_results[[paste(org, comp, sep = "_")]] <- correlate_gut_hr_pathways_metabolites(organ_name = org,
            compartment = comp)
    }
}
gut_hr_pathway_metabolite_cor <- dplyr::bind_rows(gut_hr_cor_results) %>% dplyr::group_by(organ, compartment) %>%
    dplyr::mutate(padj = p.adjust(pvalue, method = "BH")) %>% dplyr::ungroup() %>% dplyr::arrange(padj)
correlate_serum_hr_pathways_metabolites <- function(organ_name) {
    dds <- dds_list[[organ_name]]
    M <- make_serum_M_block(dds = dds, serum_combined = serum_combined, serum_meta = pos_serum$metadata,
        serum_bridges = serum_bridges)
    score_df <- pathway_scores_hr_all %>% dplyr::filter(organ == organ_name) %>% dplyr::select(sample,
        group, dplyr::any_of(sig_hr_pathways_by_org[[organ_name]]))
    pathway_cols <- intersect(sig_hr_pathways_by_org[[organ_name]], colnames(score_df))
    common <- intersect(score_df$sample, rownames(M))
    if (length(pathway_cols) == 0 || length(common) < min_n) {
        return(tibble::tibble())
    }
    score_df <- score_df %>% dplyr::filter(sample %in% common) %>% dplyr::arrange(match(sample, common))
    M <- M[score_df$sample, , drop = FALSE]
    out <- expand.grid(pathway = pathway_cols, metabolite = colnames(M), stringsAsFactors = FALSE)
    purrr::pmap_dfr(out, function(pathway, metabolite) {
        test <- safe_cor(x = score_df[[pathway]], y = M[, metabolite], min_n = min_n)
        tibble::tibble(organ = organ_name, compartment = "serum", source = "serum", pathway = pathway,
            metabolite = metabolite, rho = as.numeric(test["rho"]), pvalue = as.numeric(test["pvalue"]),
            n = as.numeric(test["n"]))
    })
}
serum_hr_cor_results <- list()
for (org in names(dds_list)) {
    message("Correlating serum: ", org)
    serum_hr_cor_results[[org]] <- correlate_serum_hr_pathways_metabolites(organ_name = org)
}
serum_hr_pathway_metabolite_cor <- dplyr::bind_rows(serum_hr_cor_results) %>% dplyr::group_by(organ) %>%
    dplyr::mutate(padj = p.adjust(pvalue, method = "BH")) %>% dplyr::ungroup() %>% dplyr::arrange(padj)
hr_pathway_metabolite_cor_all <- dplyr::bind_rows(gut_hr_pathway_metabolite_cor, serum_hr_pathway_metabolite_cor) %>%
    dplyr::arrange(padj)
hallmark_pathway_metabolite_cor_all <- hr_pathway_metabolite_cor_all %>% dplyr::filter(grepl("^HALLMARK_",
    pathway)) %>% dplyr::group_by(organ, compartment) %>% dplyr::mutate(padj_hallmark_only = p.adjust(pvalue,
    method = "BH")) %>% dplyr::ungroup() %>% dplyr::arrange(padj_hallmark_only)
hallmark_pathway_metabolite_cor_sig <- hallmark_pathway_metabolite_cor_all %>% dplyr::filter(!is.na(padj_hallmark_only),
    padj_hallmark_only < 0.1)
write.csv(hallmark_pathway_metabolite_cor_all, project_output(file.path(out_dir, "ALL_Hallmark_ONLY_pathway_metabolite_correlations_recalculated_FDR.csv"),
    name = "all_hallmark_only_pathway_metabolite_correlations_recalculated_fdr.csv"), row.names = FALSE)
write.csv(hallmark_pathway_metabolite_cor_sig, project_output(file.path(out_dir, "SIGNIFICANT_Hallmark_ONLY_pathway_metabolite_correlations_FDR0p1.csv"),
    name = "significant_hallmark_only_pathway_metabolite_correlations_fdr0p1.csv"), row.names = FALSE)
message("Hallmark-only significant correlations, recalculated FDR < 0.1:")
print(hallmark_pathway_metabolite_cor_sig %>% dplyr::count(organ, compartment, sort = TRUE))
write.csv(gut_hr_pathway_metabolite_cor, project_output(file.path(out_dir, "gut_Hallmark_Reactome_pathway_metabolite_correlations.csv"),
    name = "gut_hallmark_reactome_pathway_metabolite_correlations.csv"), row.names = FALSE)
write.csv(serum_hr_pathway_metabolite_cor, project_output(file.path(out_dir, "serum_Hallmark_Reactome_pathway_metabolite_correlations.csv"),
    name = "serum_hallmark_reactome_pathway_metabolite_correlations.csv"), row.names = FALSE)
write.csv(hr_pathway_metabolite_cor_all, project_output(file.path(out_dir, "ALL_Hallmark_Reactome_pathway_metabolite_correlations.csv"),
    name = "all_hallmark_reactome_pathway_metabolite_correlations.csv"), row.names = FALSE)
writexl::write_xlsx(list(all_correlations = hr_pathway_metabolite_cor_all, significant_FDR_0p1 = hr_pathway_metabolite_cor_all %>%
    dplyr::filter(!is.na(padj), padj < 0.1), hallmark_only_all = hallmark_pathway_metabolite_cor_all,
    hallmark_only_sig_FDR_0p1 = hallmark_pathway_metabolite_cor_sig, gut = gut_hr_pathway_metabolite_cor,
    serum = serum_hr_pathway_metabolite_cor, pathway_matching_check = matched_hr), project_output(file.path(out_dir,
    "Hallmark_Reactome_pathway_metabolite_correlations.xlsx"), name = "hallmark_reactome_pathway_metabolite_correlations.xlsx"))
saveRDS(list(pathway_scores_hr_all = pathway_scores_hr_all, gut_cor = gut_hr_pathway_metabolite_cor,
    serum_cor = serum_hr_pathway_metabolite_cor, all_cor = hr_pathway_metabolite_cor_all, hallmark_only_all = hallmark_pathway_metabolite_cor_all,
    hallmark_only_sig = hallmark_pathway_metabolite_cor_sig, matched_hr = matched_hr, cp_gsea_sig = cp_gsea_sig,
    sig_hr_pathways_by_org = sig_hr_pathways_by_org, excluded_samples = excluded_samples), project_output(file.path(out_dir,
    "Hallmark_Reactome_pathway_metabolite_correlations.rds"), name = "hallmark_reactome_pathway_metabolite_correlations.rds"))
message("Done.")
}))
