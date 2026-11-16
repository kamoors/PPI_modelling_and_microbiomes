# SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0
# Figure 5: mouse transcriptomics.
# Run through run.R so input and output paths resolve within the package.
# Inputs: see data/input_manifest.csv and docs/RUN_ORDER.md.
# Writes generated files under results/; bundled inputs are read-only.

required_packages <- c("dplyr", "tibble", "ggplot2", "ggrepel", "cowplot", "stringr", "readr")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) {
    stop("Install required packages first: ", paste(missing_packages, collapse = ", "))
}
ppi_check_plot_versions(required_packages)
FIG_ROOT <- ppi_path("figures")
DATA_DIR <- file.path(FIG_ROOT, "source_data", "Mouse_Transcriptomics")
OUT_DIR <- file.path(FIG_ROOT, "final")
dir.create(ppi_output(DATA_DIR, directory = TRUE), recursive = TRUE, showWarnings = FALSE)
dir.create(ppi_output(OUT_DIR, directory = TRUE), recursive = TRUE, showWarnings = FALSE)
REFRESH_SOURCE_DATA <- FALSE
PPI_ROOT <- ppi_path("")
RNA_ROOT <- file.path(PPI_ROOT, "rnaseq")
stage_input <- function(filename, candidates) {
    target <- file.path(DATA_DIR, filename)
    if (REFRESH_SOURCE_DATA || !ppi_file_exists(target)) {
        candidates <- unique(candidates[nzchar(candidates)])
        hits <- candidates[ppi_file_exists(candidates)]
        if (!length(hits)) {
            stop("Could not stage ", filename, ".\nChecked:\n", paste(candidates, collapse = "\n"))
        }
        if (!file.copy(ppi_input(hits[[1]]), ppi_output(target), overwrite = TRUE)) 
            stop("Failed to copy ", hits[[1]])
        message("Staged: ", filename)
    }
    target
}
# Frozen plotting tables avoid distributing fitted DESeq2 models.
cache_file <- stage_input("rna_plot_inputs.rds", character())
candidate_result_dirs <- c(file.path(RNA_ROOT, "Results", "RNAseq_without_C02_G04"), file.path(RNA_ROOT, 
    "RNAseq_without_C02_G04"))
find_candidates <- function(filename) {
    direct <- file.path(candidate_result_dirs, filename)
    recursive <- ppi_list_files(RNA_ROOT, pattern = paste0("^", gsub("([.])", "\\\\.", filename), "$"), 
        recursive = TRUE, full.names = TRUE)
    c(direct, recursive)
}
organs <- c("Liver", "Colon", "Brain")
de_files <- setNames(vector("character", length(organs)), organs)
gsea_files <- setNames(vector("character", length(organs)), organs)
for (org in organs) {
    de_name <- paste0("DESeq2_", org, "_without_C02_G04.csv")
    gsea_name <- paste0("GSEA_", org, "_HALLMARK_without_C02_G04.csv")
    de_files[[org]] <- stage_input(de_name, find_candidates(de_name))
    gsea_files[[org]] <- stage_input(gsea_name, find_candidates(gsea_name))
}
readr::write_csv(data.frame(file = basename(c(cache_file, unname(de_files), unname(gsea_files))), centralized_path = c(cache_file, 
    unname(de_files), unname(gsea_files))), ppi_output(file.path(DATA_DIR, "Mouse_Transcriptomics_source_manifest.csv"), 
    name = "mouse_transcriptomics_source_manifest.csv"))
BASE_SIZE <- 18
FIG_W_IN <- 18
FIG_H_IN <- 14.5
DPI <- 600
group_cols <- c(Ctrl = "#FF6E00", PPI = "#08306B")
direction_cols <- c(Down = "dodgerblue", NS = "grey70", Up = "firebrick")
base_theme <- ggplot2::theme_classic(base_size = BASE_SIZE, base_family = "Arial") + ggplot2::theme(axis.text = ggplot2::element_text(size = 16, 
    colour = "#2F2F2F"), axis.title = ggplot2::element_text(size = 18, colour = "#2F2F2F"), plot.title = ggplot2::element_text(size = 19, 
    hjust = 0.5, face = "bold"), legend.title = ggplot2::element_text(size = 16), legend.text = ggplot2::element_text(size = 15), 
    plot.margin = ggplot2::margin(7, 7, 7, 7))
plot_inputs <- readRDS(ppi_input(cache_file))
pca_for_organ <- function(organ_name) {
    hit <- names(plot_inputs$pca)[tolower(names(plot_inputs$pca)) == tolower(organ_name)]
    if (length(hit) != 1L) stop("Missing PCA table for ", organ_name)
    pdat <- plot_inputs$pca[[hit]]
    percentVar <- round(100 * attr(pdat, "percentVar"))
    ggplot2::ggplot(pdat, ggplot2::aes(PC1, PC2, colour = group)) + ggplot2::geom_point(size = 3.4, alpha = 0.9) + 
        ggplot2::scale_colour_manual(name = "Group", values = group_cols) + ggplot2::labs(title = organ_name, 
        x = paste0("PC1: ", percentVar[1], "%"), y = paste0("PC2: ", percentVar[2], "%")) + base_theme + 
        ggplot2::theme(legend.position = "bottom")
}
add_gene_symbols <- function(res_df) {
    map <- plot_inputs$gene_annotation
    if (!"geneSymbol" %in% names(res_df))
        res_df$geneSymbol <- map$geneSymbol[match(res_df$gene, map$gene)]
    res_df
}
volcano_for_organ <- function(res_df, alpha = 0.1, n_labels = 5) {
    df <- dplyr::mutate(add_gene_symbols(res_df), signif = dplyr::case_when(!is.na(padj) & padj < alpha & 
        log2FoldChange > 0 ~ "Up", !is.na(padj) & padj < alpha & log2FoldChange < 0 ~ "Down", TRUE ~ 
        "NS"), label = ifelse(!is.na(geneSymbol) & geneSymbol != "", geneSymbol, gene), neg_log10_padj = -log10(pmax(padj, 
        .Machine$double.xmin)))
    topGenes <- dplyr::slice_head(dplyr::arrange(dplyr::filter(df, !is.na(padj), is.finite(neg_log10_padj), 
        !is.na(log2FoldChange)), padj), n = n_labels)
    ggplot2::ggplot(df, ggplot2::aes(log2FoldChange, neg_log10_padj, colour = signif)) + ggplot2::geom_point(alpha = 0.6, 
        size = 1.25) + ggplot2::geom_hline(yintercept = -log10(alpha), linetype = "dashed", colour = "black", 
        linewidth = 0.35) + ggrepel::geom_text_repel(data = topGenes, ggplot2::aes(label = label), size = 6, 
        max.overlaps = Inf, segment.colour = "#9E9E9E") + ggplot2::scale_colour_manual(name = "Direction", 
        values = direction_cols) + ggplot2::labs(x = "log2FC", y = expression(-log[10]("FDR"))) + base_theme + 
        ggplot2::theme(legend.position = "bottom")
}
gsea_dotplot_for_organ <- function(fgsea_df, show_n = 20) {
    df <- dplyr::arrange(dplyr::filter(dplyr::mutate(fgsea_df, direction = ifelse(NES > 0, "Positive", 
        "Negative"), pathway_clean = stringr::str_to_title(gsub("_", " ", gsub("^HALLMARK_", "", pathway))), 
        neg_log10_padj = -log10(pmax(padj, .Machine$double.xmin))), !is.na(padj), is.finite(neg_log10_padj)), 
        padj)
    pos <- dplyr::slice_min(dplyr::filter(df, NES > 0), padj, n = ceiling(show_n/2))
    neg <- dplyr::slice_min(dplyr::filter(df, NES < 0), padj, n = floor(show_n/2))
    plot_df <- dplyr::mutate(dplyr::arrange(dplyr::bind_rows(pos, neg), NES), pathway_clean = factor(pathway_clean, 
        levels = pathway_clean))
    ggplot2::ggplot(plot_df, ggplot2::aes(NES, pathway_clean, size = neg_log10_padj, colour = direction)) + 
        ggplot2::geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50", linewidth = 0.35) + 
        ggplot2::geom_point(alpha = 0.88) + ggplot2::scale_colour_manual(name = "Direction", values = c(Positive = "firebrick", 
        Negative = "dodgerblue")) + ggplot2::scale_size_continuous(range = c(2.2, 5.2), name = "-log10(FDR)") + 
        ggplot2::labs(x = "NES", y = NULL) + base_theme + ggplot2::theme(legend.position = "bottom", 
        axis.text.y = ggplot2::element_text(size = 15))
}
pca_plots <- lapply(organs, function(org) pca_for_organ(org))
volcano_plots <- lapply(organs, function(org) {
    res <- readr::read_csv(ppi_input(de_files[[org]]), show_col_types = FALSE)
    volcano_for_organ(res)
})
gsea_plots <- lapply(organs, function(org) {
    fg <- readr::read_csv(ppi_input(gsea_files[[org]]), show_col_types = FALSE)
    gsea_dotplot_for_organ(fg, show_n = 20)
})
pca_legend <- cowplot::get_legend(pca_plots[[1]] + ggplot2::theme(legend.position = "bottom"))
volcano_legend <- cowplot::get_legend(volcano_plots[[1]] + ggplot2::theme(legend.position = "bottom"))
gsea_legend <- cowplot::get_legend(gsea_plots[[1]] + ggplot2::theme(legend.position = "bottom"))
pca_row_plots <- cowplot::plot_grid(plotlist = lapply(pca_plots, function(x) x + ggplot2::theme(legend.position = "none")), 
    ncol = 3, labels = c("A", "B", "C"), label_size = 24, label_fontface = "bold")
pca_row <- cowplot::plot_grid(pca_row_plots, pca_legend, ncol = 1, rel_heights = c(1, 0.1))
volcano_row_plots <- cowplot::plot_grid(plotlist = lapply(volcano_plots, function(x) x + ggplot2::theme(legend.position = "none")), 
    ncol = 3, labels = c("D", "E", "F"), label_size = 24, label_fontface = "bold")
volcano_row <- cowplot::plot_grid(volcano_row_plots, volcano_legend, ncol = 1, rel_heights = c(1, 0.1))
gsea_row_plots <- cowplot::plot_grid(plotlist = lapply(gsea_plots, function(x) x + ggplot2::theme(legend.position = "none")), 
    ncol = 3, labels = c("G", "H", "I"), label_size = 24, label_fontface = "bold")
gsea_row <- cowplot::plot_grid(gsea_row_plots, gsea_legend, ncol = 1, rel_heights = c(1, 0.12))
final_figure <- cowplot::plot_grid(pca_row, volcano_row, gsea_row, ncol = 1, rel_heights = c(1, 1, 1.72))
print(final_figure)
pdf_file <- file.path(OUT_DIR, "Mouse_Transcriptomics.pdf")
png_file <- file.path(OUT_DIR, "Mouse_Transcriptomics.png")
grDevices::cairo_pdf(filename = ppi_output(pdf_file), width = FIG_W_IN, height = FIG_H_IN)
print(final_figure)
grDevices::dev.off()
if (requireNamespace("ragg", quietly = TRUE)) {
    ragg::agg_png(ppi_output(png_file), width = FIG_W_IN, height = FIG_H_IN, units = "in", res = DPI, 
        background = "white")
    print(final_figure)
    grDevices::dev.off()
} else {
    grDevices::png(ppi_output(png_file), width = FIG_W_IN, height = FIG_H_IN, units = "in", res = DPI, 
        bg = "white")
    print(final_figure)
    grDevices::dev.off()
}
message("Figure 5 written to: ", OUT_DIR)
