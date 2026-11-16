# SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0
# Outlier-excluded DESeq2 and GSEA
# Run through run.R so input and output paths resolve within the package.
# Inputs: see data/input_manifest.csv and docs/RUN_ORDER.md.
# Writes generated files under results/; bundled inputs are read-only.

library(DESeq2)
library(SummarizedExperiment)
library(dplyr)
library(tibble)
library(ggplot2)
library(pheatmap)
library(clusterProfiler)
library(org.Mm.eg.db)
library(AnnotationDbi)
library(writexl)
library(msigdbr)
library(fgsea)
load_rnaseq_script <- file.path(ppi_code("scripts/05_mouse_transcriptomics/05_02_load_primary_deseq2_cache.R"))
ppi_source(load_rnaseq_script)
out_dir <- file.path(getwd(), "RNAseq_without_C02_G04")
dir.create(ppi_output(out_dir, directory = TRUE), recursive = TRUE, showWarnings = FALSE)
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
    write.csv(res_no_outliers_list[[org]], ppi_output(file.path(out_dir, paste0("DESeq2_", org, "_without_C02_G04.csv"))), 
        row.names = FALSE)
}
res_no_outliers_all <- dplyr::bind_rows(res_no_outliers_list)
write.csv(res_no_outliers_all, ppi_output(file.path(out_dir, "DESeq2_all_organs_without_C02_G04.csv"), 
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
write.csv(compare_to_original, ppi_output(file.path(out_dir, "DESeq2_original_vs_without_C02_G04.csv"), 
    name = "de_seq2_original_vs_without_c02_g04.csv"), row.names = FALSE)
summary_tbl <- compare_to_original %>% dplyr::group_by(organ) %>% dplyr::summarise(n_sig_original = sum(sig_original, 
    na.rm = TRUE), n_sig_without_outliers = sum(sig_without_outliers, na.rm = TRUE), overlap = sum(sig_original & 
    sig_without_outliers, na.rm = TRUE), only_original = sum(sig_original & !sig_without_outliers, na.rm = TRUE), 
    only_without_outliers = sum(!sig_original & sig_without_outliers, na.rm = TRUE), cor_log2FC = cor(log2FC_original, 
        log2FC_without_outliers, use = "pairwise.complete.obs"), .groups = "drop")
write.csv(summary_tbl, ppi_output(file.path(out_dir, "DESeq2_summary_without_C02_G04.csv"), name = "de_seq2_summary_without_c02_g04.csv"), 
    row.names = FALSE)
print(summary_tbl)
for (org in names(deseq_no_outliers)) {
    dds_sub <- deseq_no_outliers[[org]]$dds
    vsd <- DESeq2::vst(dds_sub, blind = TRUE)
    pca_df <- DESeq2::plotPCA(vsd, intgroup = "group", returnData = TRUE)
    pca_df$sample <- rownames(pca_df)
    percentVar <- round(100 * attr(pca_df, "percentVar"))
    p <- ggplot2::ggplot(pca_df, ggplot2::aes(PC1, PC2, color = group, label = sample)) + ggplot2::geom_point(size = 4) + 
        ggplot2::geom_text(vjust = -0.5, size = 3, check_overlap = TRUE) + ggplot2::theme_bw() + ggplot2::labs(title = paste0(org, 
        " PCA without C02/G04"), x = paste0("PC1: ", percentVar[1], "%"), y = paste0("PC2: ", percentVar[2], 
        "%"))
    ggplot2::ggsave(ppi_output(file.path(out_dir, paste0("PCA_", org, "_without_C02_G04.pdf"))), p, width = 7, 
        height = 5)
}
library(DESeq2)
library(ggplot2)
library(cowplot)
library(ggrepel)
library(dplyr)
library(tibble)
library(vegan)
library(AnnotationDbi)
library(org.Mm.eg.db)
out_dir <- file.path(getwd(), "RNAseq_without_C02_G04")
organs <- c("Liver", "Colon", "Brain")
group_cols <- c(Ctrl = "#FF6E00", PPI = "#08306B")
direction_cols <- c(Down = "dodgerblue", NS = "grey70", Up = "firebrick")
permanova_from_vsd <- function(vsd, group_var = "group", permutations = 9999) {
    X <- t(SummarizedExperiment::assay(vsd))
    d <- dist(X, method = "euclidean")
    meta_df <- as.data.frame(SummarizedExperiment::colData(vsd))
    form <- as.formula(paste("d ~", group_var))
    fit <- vegan::adonis2(formula = form, data = meta_df, permutations = permutations)
    list(r2 = unname(fit$R2[1]), p = unname(fit$`Pr(>F)`[1]), fit = fit)
}
pca_for_organ <- function(dds, organ_name, intgroup = "group", permutations = 9999, label_pos = c("tl", 
    "tr", "bl", "br")) {
    label_pos <- match.arg(label_pos)
    vsd <- DESeq2::vst(dds, blind = FALSE)
    pdat <- DESeq2::plotPCA(vsd, intgroup = intgroup, returnData = TRUE, ntop = min(20000, nrow(vsd)))
    percentVar <- round(100 * attr(pdat, "percentVar"))
    pv <- permanova_from_vsd(vsd = vsd, group_var = intgroup, permutations = permutations)
    lab <- paste0("PERMANOVA\n", "R² = ", sprintf("%.3f", pv$r2), "\np = ", format.pval(pv$p, digits = 3, 
        eps = 0.001))
    pos <- switch(label_pos, tl = list(x = -Inf, y = Inf, hjust = -0.05, vjust = 1.05), tr = list(x = Inf, 
        y = Inf, hjust = 1.05, vjust = 1.05), bl = list(x = -Inf, y = -Inf, hjust = -0.05, vjust = -0.05), 
        br = list(x = Inf, y = -Inf, hjust = 1.05, vjust = -0.05))
    ggplot(pdat, aes(PC1, PC2, color = !!rlang::sym(intgroup))) + geom_point(size = 3) + annotate("text", 
        x = pos$x, y = pos$y, label = lab, hjust = pos$hjust, vjust = pos$vjust, size = 3) + scale_color_manual(name = "Group", 
        values = group_cols) + xlab(paste0("PC1: ", percentVar[1], "%")) + ylab(paste0("PC2: ", percentVar[2], 
        "%")) + ggtitle(organ_name) + theme_classic() + theme(legend.position = "bottom", plot.title = element_text(hjust = 0.5, 
        face = "bold"))
}
add_gene_symbols <- function(res_df) {
    map <- AnnotationDbi::select(org.Mm.eg.db, keys = res_df$gene, keytype = "ENSEMBL", columns = c("SYMBOL", 
        "GENENAME")) %>% dplyr::filter(!is.na(SYMBOL), SYMBOL != "") %>% dplyr::distinct(ENSEMBL, .keep_all = TRUE)
    res_df %>% dplyr::left_join(map, by = c(gene = "ENSEMBL")) %>% dplyr::rename(geneSymbol = SYMBOL, 
        geneName = GENENAME)
}
volcano_for_organ <- function(res_df, organ_name, alpha = 0.1, n_labels = 5) {
    df <- res_df %>% dplyr::mutate(signif = dplyr::case_when(!is.na(padj) & padj < alpha & log2FoldChange > 
        0 ~ "Up", !is.na(padj) & padj < alpha & log2FoldChange < 0 ~ "Down", TRUE ~ "NS"), label = dplyr::case_when(!is.na(geneSymbol) & 
        geneSymbol != "" ~ geneSymbol, TRUE ~ gene), neg_log10_padj = -log10(padj))
    topGenes <- df %>% dplyr::filter(!is.na(padj), is.finite(neg_log10_padj), !is.na(log2FoldChange)) %>% 
        dplyr::arrange(padj) %>% dplyr::slice_head(n = n_labels)
    ggplot(df, aes(x = log2FoldChange, y = neg_log10_padj, color = signif)) + geom_point(alpha = 0.6, 
        size = 1.2) + geom_hline(yintercept = -log10(alpha), linetype = "dashed", color = "black") + 
        ggrepel::geom_text_repel(data = topGenes, aes(label = label), size = 2.7, max.overlaps = Inf) + 
        scale_color_manual(name = "Direction", values = direction_cols) + labs(x = "log2FC", y = "-log10(FDR)") + 
        theme_classic() + theme(legend.position = "bottom")
}
gsea_dotplot_for_organ <- function(fgsea_df, organ_name, padj_cutoff = 0.1, show_n = 10) {
    df <- fgsea_df %>% dplyr::mutate(direction = ifelse(NES > 0, "Positive", "Negative"), pathway_clean = gsub("^HALLMARK_", 
        "", pathway), pathway_clean = gsub("_", " ", pathway_clean), pathway_clean = stringr::str_to_title(pathway_clean), 
        neg_log10_padj = -log10(padj)) %>% dplyr::filter(!is.na(padj), is.finite(neg_log10_padj)) %>% 
        dplyr::arrange(padj)
    pos <- df %>% dplyr::filter(NES > 0) %>% dplyr::slice_min(padj, n = ceiling(show_n/2))
    neg <- df %>% dplyr::filter(NES < 0) %>% dplyr::slice_min(padj, n = floor(show_n/2))
    plot_df <- dplyr::bind_rows(pos, neg) %>% dplyr::arrange(NES) %>% dplyr::mutate(pathway_clean = factor(pathway_clean, 
        levels = pathway_clean))
    ggplot(plot_df, aes(x = NES, y = pathway_clean, size = neg_log10_padj, color = direction)) + geom_vline(xintercept = 0, 
        linetype = "dashed", color = "grey50") + geom_point(alpha = 0.85) + scale_color_manual(name = "Direction", 
        values = c(Positive = "firebrick", Negative = "dodgerblue")) + labs(x = "NES", y = NULL, size = "-log10(FDR)") + 
        theme_classic() + theme(legend.position = "bottom", axis.text.y = element_text(size = 8))
}
pca_plots <- lapply(organs, function(org) {
    pca_for_organ(dds = deseq_no_outliers[[org]]$dds, organ_name = org, label_pos = "br")
})
volcano_plots <- lapply(organs, function(org) {
    res_df <- deseq_no_outliers[[org]]$res_df %>% add_gene_symbols()
    volcano_for_organ(res_df = res_df, organ_name = org, alpha = 0.1, n_labels = 5)
})
gsea_plots <- lapply(organs, function(org) {
    gsea_dotplot_for_organ(fgsea_df = gsea_no_outliers[[org]]$HALLMARK, organ_name = org, show_n = 20)
})
library(cowplot)
pca_row_plots <- plot_grid(plotlist = lapply(pca_plots, function(x) x + theme(legend.position = "none")), 
    ncol = 3, labels = c("A", "B", "C"), label_size = 16, label_fontface = "bold")
pca_row <- plot_grid(pca_row_plots, pca_legend, ncol = 1, rel_heights = c(1, 0.12))
volcano_row_plots <- plot_grid(plotlist = lapply(volcano_plots, function(x) x + theme(legend.position = "none")), 
    ncol = 3, labels = c("D", "E", "F"), label_size = 16, label_fontface = "bold")
volcano_row <- plot_grid(volcano_row_plots, volcano_legend, ncol = 1, rel_heights = c(1, 0.12))
gsea_row_plots <- plot_grid(plotlist = lapply(gsea_plots, function(x) x + theme(legend.position = "none")), 
    ncol = 3, labels = c("G", "H", "I"), label_size = 16, label_fontface = "bold")
gsea_row <- plot_grid(gsea_row_plots, gsea_legend, ncol = 1, rel_heights = c(1, 0.15))
final_panel <- plot_grid(pca_row, volcano_row, gsea_row, ncol = 1, rel_heights = c(1, 1, 1.8))
print(final_panel)
ggsave(ppi_output(file.path(out_dir, "PCA_volcano_HallmarkGSEA_panel_without_C02_G04.pdf"), name = "pca_volcano_hallmark_gsea_panel_without_c02_g04.pdf"), 
    final_panel, width = 5.2 * length(organs), height = 12)
ggsave(ppi_output(file.path(out_dir, "PCA_volcano_HallmarkGSEA_panel_without_C02_G04.png"), name = "pca_volcano_hallmark_gsea_panel_without_c02_g04.png"), 
    final_panel, width = 5.2 * length(organs), height = 12, dpi = 300)
map_ens_to_entrez <- function(genes) {
    AnnotationDbi::select(org.Mm.eg.db, keys = genes, keytype = "ENSEMBL", columns = c("ENTREZID", "SYMBOL")) %>% 
        dplyr::filter(!is.na(ENTREZID)) %>% dplyr::distinct(ENSEMBL, .keep_all = TRUE)
}
go_results <- list()
for (org in names(deseq_no_outliers)) {
    res_df <- deseq_no_outliers[[org]]$res_df
    degs <- res_df %>% dplyr::filter(sig_0.1) %>% dplyr::pull(gene)
    bg <- res_df$gene
    if (length(degs) < 10) {
        message("Skipping GO for ", org, ": too few DEGs")
        next
    }
    deg_map <- map_ens_to_entrez(degs)
    bg_map <- map_ens_to_entrez(bg)
    ego <- clusterProfiler::enrichGO(gene = unique(deg_map$ENTREZID), universe = unique(bg_map$ENTREZID), 
        OrgDb = org.Mm.eg.db, keyType = "ENTREZID", ont = "BP", pAdjustMethod = "BH", readable = TRUE)
    go_results[[org]] <- ego
    write.csv(as.data.frame(ego), ppi_output(file.path(out_dir, paste0("GO_BP_", org, "_without_C02_G04.csv"))), 
        row.names = FALSE)
    pdf(ppi_output(file.path(out_dir, paste0("GO_BP_", org, "_without_C02_G04_dotplot.pdf"))), width = 8, 
        height = 6)
    print(clusterProfiler::dotplot(ego, showCategory = 20) + ggplot2::ggtitle(paste(org, "GO BP")))
    dev.off()
}
saveRDS(list(deseq_no_outliers = deseq_no_outliers, res_no_outliers_all = res_no_outliers_all, compare_to_original = compare_to_original, 
    summary_tbl = summary_tbl, go_results = go_results, exclude_samples = exclude_samples), ppi_output(file.path(out_dir, 
    "RNAseq_without_C02_G04_full_objects.rds"), name = "rn_aseq_without_c02_g04_full_objects.rds"))
writexl::write_xlsx(list(summary = summary_tbl, comparison = compare_to_original, DESeq2_all_organs = res_no_outliers_all), 
    ppi_output(file.path(out_dir, "RNAseq_without_C02_G04_summary.xlsx"), name = "rn_aseq_without_c02_g04_summary.xlsx"))
message("DESeq2 and GO enrichment done. Results saved to: ", out_dir)
msig_hallmark <- msigdbr::msigdbr(species = "Mus musculus", category = "H") %>% dplyr::select(gs_name, 
    gene_symbol)
msig_reactome <- msigdbr::msigdbr(species = "Mus musculus", category = "C2", subcategory = "CP:REACTOME") %>% 
    dplyr::select(gs_name, gene_symbol)
hallmark_sets <- split(msig_hallmark$gene_symbol, msig_hallmark$gs_name)
hallmark_sets <- lapply(hallmark_sets, unique)
reactome_sets <- split(msig_reactome$gene_symbol, msig_reactome$gs_name)
reactome_sets <- lapply(reactome_sets, unique)
make_ranked_symbol_vector <- function(res_df) {
    map <- AnnotationDbi::select(org.Mm.eg.db, keys = res_df$gene, keytype = "ENSEMBL", columns = c("SYMBOL")) %>% 
        dplyr::filter(!is.na(SYMBOL), SYMBOL != "") %>% dplyr::distinct(ENSEMBL, .keep_all = TRUE)
    tmp <- res_df %>% dplyr::left_join(map, by = c(gene = "ENSEMBL")) %>% dplyr::filter(!is.na(SYMBOL), 
        !is.na(stat), is.finite(stat)) %>% dplyr::group_by(SYMBOL) %>% dplyr::slice_max(order_by = abs(stat), 
        n = 1, with_ties = FALSE) %>% dplyr::ungroup()
    ranks <- tmp$stat
    names(ranks) <- tmp$SYMBOL
    sort(ranks, decreasing = TRUE)
}
run_fgsea_one <- function(pathways, ranks, label, minSize = 15, maxSize = 500, eps = 1e-06) {
    pathways_use <- lapply(pathways, function(g) intersect(g, names(ranks)))
    pathways_use <- pathways_use[lengths(pathways_use) >= minSize & lengths(pathways_use) <= maxSize]
    fg <- fgsea::fgseaMultilevel(pathways = pathways_use, stats = ranks, minSize = minSize, maxSize = maxSize, 
        eps = eps, nproc = 1) %>% as.data.frame() %>% dplyr::arrange(padj) %>% dplyr::mutate(leadingEdge = vapply(leadingEdge, 
        function(x) paste(x, collapse = ";"), character(1)))
    write.csv(fg, ppi_output(file.path(out_dir, paste0("GSEA_", label, ".csv"))), row.names = FALSE)
    fg
}
gsea_no_outliers <- list()
ranks_no_outliers <- list()
for (org in names(deseq_no_outliers)) {
    message("Running outlier-excluded GSEA for: ", org)
    res_df <- deseq_no_outliers[[org]]$res_df
    ranks <- make_ranked_symbol_vector(res_df)
    ranks_no_outliers[[org]] <- ranks
    fg_h <- run_fgsea_one(pathways = hallmark_sets, ranks = ranks, label = paste0(org, "_HALLMARK_without_C02_G04"), 
        minSize = 15, maxSize = 300)
    fg_r <- run_fgsea_one(pathways = reactome_sets, ranks = ranks, label = paste0(org, "_REACTOME_without_C02_G04"), 
        minSize = 15, maxSize = 500)
    gsea_no_outliers[[org]] <- list(HALLMARK = fg_h, REACTOME = fg_r)
}
saveRDS(gsea_no_outliers, ppi_output(file.path(out_dir, "GSEA_without_C02_G04_results.rds"), name = "gsea_without_c02_g04_results.rds"))
saveRDS(ranks_no_outliers, ppi_output(file.path(out_dir, "ranked_gene_vectors_without_C02_G04.rds"), 
    name = "ranked_gene_vectors_without_c02_g04.rds"))
gsea_original <- list()
ranks_original <- list()
for (org in names(res_list)) {
    message("Running original GSEA for comparison: ", org)
    res_orig_df <- as.data.frame(res_list[[org]]) %>% tibble::rownames_to_column("gene")
    ranks_orig <- make_ranked_symbol_vector(res_orig_df)
    ranks_original[[org]] <- ranks_orig
    fg_h_orig <- run_fgsea_one(pathways = hallmark_sets, ranks = ranks_orig, label = paste0(org, "_HALLMARK_original"), 
        minSize = 15, maxSize = 300)
    fg_r_orig <- run_fgsea_one(pathways = reactome_sets, ranks = ranks_orig, label = paste0(org, "_REACTOME_original"), 
        minSize = 15, maxSize = 500)
    gsea_original[[org]] <- list(HALLMARK = fg_h_orig, REACTOME = fg_r_orig)
}
saveRDS(gsea_original, ppi_output(file.path(out_dir, "GSEA_original_results.rds"), name = "gsea_original_results.rds"))
compare_gsea_one <- function(fg_orig, fg_new, organ, library_name) {
    comp <- fg_orig %>% dplyr::select(pathway, NES_original = NES, pval_original = pval, padj_original = padj) %>% 
        dplyr::inner_join(fg_new %>% dplyr::select(pathway, NES_without_outliers = NES, pval_without_outliers = pval, 
            padj_without_outliers = padj), by = "pathway") %>% dplyr::mutate(organ = organ, library = library_name, 
        sig_original = padj_original < 0.1, sig_without_outliers = padj_without_outliers < 0.1, status = dplyr::case_when(sig_original & 
            sig_without_outliers ~ "Significant_both", sig_original & !sig_without_outliers ~ "Only_original", 
            !sig_original & sig_without_outliers ~ "Only_without_outliers", TRUE ~ "Not_significant"), 
        delta_NES = NES_without_outliers - NES_original, flipped_direction = sign(NES_original) != sign(NES_without_outliers)) %>% 
        dplyr::arrange(padj_without_outliers)
    write.csv(comp, ppi_output(file.path(out_dir, paste0("GSEA_comparison_", organ, "_", library_name, 
        "_original_vs_without_C02_G04.csv"))), row.names = FALSE)
    p <- ggplot2::ggplot(comp, ggplot2::aes(x = NES_original, y = NES_without_outliers, color = status)) + 
        ggplot2::geom_point(size = 2, alpha = 0.8) + ggplot2::geom_abline(slope = 1, intercept = 0, linetype = 2) + 
        ggplot2::theme_bw() + ggplot2::labs(title = paste(organ, library_name, "GSEA NES comparison"), 
        subtitle = paste0("NES correlation = ", round(cor(comp$NES_original, comp$NES_without_outliers, 
            use = "complete.obs"), 3)), x = "NES original", y = "NES without C02/G04")
    ggplot2::ggsave(ppi_output(file.path(out_dir, paste0("GSEA_NES_comparison_", organ, "_", library_name, 
        ".pdf"))), p, width = 6, height = 5)
    comp
}
gsea_comparisons <- list()
for (org in names(gsea_no_outliers)) {
    gsea_comparisons[[paste0(org, "_HALLMARK")]] <- compare_gsea_one(fg_orig = gsea_original[[org]]$HALLMARK, 
        fg_new = gsea_no_outliers[[org]]$HALLMARK, organ = org, library_name = "HALLMARK")
    gsea_comparisons[[paste0(org, "_REACTOME")]] <- compare_gsea_one(fg_orig = gsea_original[[org]]$REACTOME, 
        fg_new = gsea_no_outliers[[org]]$REACTOME, organ = org, library_name = "REACTOME")
}
gsea_comparison_all <- dplyr::bind_rows(gsea_comparisons)
write.csv(gsea_comparison_all, ppi_output(file.path(out_dir, "GSEA_comparison_all_organs_original_vs_without_C02_G04.csv"), 
    name = "gsea_comparison_all_organs_original_vs_without_c02_g04.csv"), row.names = FALSE)
gsea_summary <- gsea_comparison_all %>% dplyr::group_by(organ, library) %>% dplyr::summarise(n_sig_original = sum(sig_original, 
    na.rm = TRUE), n_sig_without_outliers = sum(sig_without_outliers, na.rm = TRUE), n_sig_both = sum(sig_original & 
    sig_without_outliers, na.rm = TRUE), n_only_original = sum(sig_original & !sig_without_outliers, 
    na.rm = TRUE), n_only_without_outliers = sum(!sig_original & sig_without_outliers, na.rm = TRUE), 
    n_flipped_direction = sum(flipped_direction, na.rm = TRUE), NES_correlation = cor(NES_original, NES_without_outliers, 
        use = "complete.obs"), .groups = "drop")
write.csv(gsea_summary, ppi_output(file.path(out_dir, "GSEA_summary_original_vs_without_C02_G04.csv"), 
    name = "gsea_summary_original_vs_without_c02_g04.csv"), row.names = FALSE)
print(gsea_summary)
gsea_flipped_important <- gsea_comparison_all %>% dplyr::filter(flipped_direction, sig_original | sig_without_outliers) %>% 
    dplyr::arrange(organ, library, padj_without_outliers)
write.csv(gsea_flipped_important, ppi_output(file.path(out_dir, "GSEA_direction_flips_significant_pathways.csv"), 
    name = "gsea_direction_flips_significant_pathways.csv"), row.names = FALSE)
gsea_robust <- gsea_comparison_all %>% dplyr::filter(sig_original, sig_without_outliers, sign(NES_original) == 
    sign(NES_without_outliers)) %>% dplyr::arrange(organ, library, padj_without_outliers)
write.csv(gsea_robust, ppi_output(file.path(out_dir, "GSEA_robust_significant_both_same_direction.csv"), 
    name = "gsea_robust_significant_both_same_direction.csv"), row.names = FALSE)
gsea_revealed <- gsea_comparison_all %>% dplyr::filter(!sig_original, sig_without_outliers) %>% dplyr::arrange(organ, 
    library, padj_without_outliers)
write.csv(gsea_revealed, ppi_output(file.path(out_dir, "GSEA_revealed_without_C02_G04.csv"), name = "gsea_revealed_without_c02_g04.csv"), 
    row.names = FALSE)
writexl::write_xlsx(list(DESeq2_summary = summary_tbl, GSEA_summary = gsea_summary, GSEA_robust = gsea_robust, 
    GSEA_flipped_important = gsea_flipped_important, GSEA_revealed = gsea_revealed), ppi_output(file.path(out_dir, 
    "RNAseq_without_C02_G04_GSEA_summary.xlsx"), name = "rn_aseq_without_c02_g04_gsea_summary.xlsx"))
saveRDS(list(deseq_no_outliers = deseq_no_outliers, compare_to_original = compare_to_original, summary_tbl = summary_tbl, 
    go_results = go_results, gsea_no_outliers = gsea_no_outliers, gsea_original = gsea_original, gsea_comparisons = gsea_comparisons, 
    gsea_summary = gsea_summary, gsea_robust = gsea_robust, gsea_flipped_important = gsea_flipped_important, 
    gsea_revealed = gsea_revealed, exclude_samples = exclude_samples), ppi_output(file.path(out_dir, 
    "RNAseq_without_C02_G04_GSEA_full_objects.rds"), name = "rn_aseq_without_c02_g04_gsea_full_objects.rds"))
message("RNA-seq sensitivity analysis complete. Results saved to: ", out_dir)
