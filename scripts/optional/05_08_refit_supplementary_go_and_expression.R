# SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0
# Liver GO and selected expression

source("scripts/helpers/setup.R")
invisible(run_analysis(
  section = "05_mouse_transcriptomics",
  directory = "figures/scripts",
  script = "05_08_refit_supplementary_go_and_expression.R",
  restricted = FALSE,
  code = {
required_packages <- c("DESeq2", "SummarizedExperiment", "dplyr", "tidyr", "tibble", "readr", "ggplot2",
    "cowplot", "stringr", "AnnotationDbi", "org.Mm.eg.db", "clusterProfiler")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) {
    stop("Install required packages first: ", paste(missing_packages, collapse = ", "))
}
FIG_ROOT <- project_path("figures")
PPI_ROOT <- project_path("")
RNA_ROOT <- file.path(PPI_ROOT, "rnaseq")
DATA_DIR <- file.path(FIG_ROOT, "source_data", "Mouse_Transcriptomics_Supplementary")
OUT_DIR <- file.path(FIG_ROOT, "final")
dir.create(project_output(DATA_DIR, directory = TRUE), recursive = TRUE, showWarnings = FALSE)
dir.create(project_output(OUT_DIR, directory = TRUE), recursive = TRUE, showWarnings = FALSE)
REFRESH_SOURCE_DATA <- FALSE
stage_input <- function(filename, candidates) {
    target <- file.path(DATA_DIR, filename)
    if (REFRESH_SOURCE_DATA || !project_file_exists(target)) {
        candidates <- unique(candidates[nzchar(candidates)])
        hits <- candidates[project_file_exists(candidates)]
        if (!length(hits)) {
            stop("Could not stage ", filename, ".\nChecked:\n", paste(candidates, collapse = "\n"))
        }
        if (!file.copy(project_input(hits[[1]]), project_output(target), overwrite = TRUE)) {
            stop("Failed to copy ", hits[[1]])
        }
        message("Staged: ", filename)
    }
    target
}
cache_file <- stage_input("rnaseq_deseq2_cache.rds", file.path(RNA_ROOT, "cache", "rnaseq_deseq2_cache.rds"))
candidate_result_dirs <- c(file.path(RNA_ROOT, "Results", "RNAseq_without_C02_G04"), file.path(RNA_ROOT,
    "RNAseq_without_C02_G04"))
find_candidates <- function(filename) {
    direct <- file.path(candidate_result_dirs, filename)
    recursive <- project_list_files(RNA_ROOT, pattern = paste0("^", gsub("([.])", "\\\\.", filename), "$"),
        recursive = TRUE, full.names = TRUE)
    c(direct, recursive)
}
organs <- c("Liver", "Colon", "Brain")
de_files <- setNames(vector("character", length(organs)), organs)
for (org in organs) {
    de_name <- paste0("DESeq2_", org, "_without_C02_G04.csv")
    de_files[[org]] <- stage_input(de_name, find_candidates(de_name))
}
readr::write_csv(data.frame(file = basename(c(cache_file, unname(de_files))), centralized_path = c(cache_file,
    unname(de_files))), project_output(file.path(DATA_DIR, "Mouse_Transcriptomics_Supplementary_source_manifest.csv"),
    name = "mouse_transcriptomics_supplementary_source_manifest.csv"))
BASE_SIZE <- 18
FIG_W_IN <- 18
FIG_H_IN <- 14.5
DPI <- 600
ALPHA <- 0.1
TOP_GO_PER_DIRECTION <- 10
group_cols <- c(Ctrl = "#FF6E00", PPI = "#08306B")
direction_cols <- c(`Higher in PPI` = "firebrick", `Lower in PPI` = "dodgerblue")
base_theme <- ggplot2::theme_classic(base_size = BASE_SIZE, base_family = "Arial") + ggplot2::theme(axis.text = ggplot2::element_text(size = 16,
    colour = "#2F2F2F"), axis.title = ggplot2::element_text(size = 18, colour = "#2F2F2F"), plot.title = ggplot2::element_text(size = 19,
    hjust = 0.5, face = "bold"), legend.title = ggplot2::element_text(size = 16), legend.text = ggplot2::element_text(size = 15),
    plot.margin = ggplot2::margin(7, 7, 7, 7))
exclude_samples <- c("25Jul1866-C02", "25Jul1866-G04")
rnaseq_cache <- readRDS(project_input(cache_file))
if (!"dds_list" %in% names(rnaseq_cache)) {
    stop("Cache does not contain dds_list.")
}
dds_list <- rnaseq_cache$dds_list
prepare_dds <- function(dds) {
    dds <- dds[, !colnames(dds) %in% exclude_samples]
    dds$group <- droplevels(factor(dds$group))
    dds$group <- relevel(dds$group, ref = "Ctrl")
    DESeq2::design(dds) <- ~group
    keep <- rowSums(DESeq2::counts(dds) >= 10) >= ceiling(0.2 * ncol(dds))
    dds[keep, ]
}
dds_no_outliers <- lapply(dds_list, prepare_dds)
get_dds <- function(org) {
    hit <- names(dds_no_outliers)[tolower(names(dds_no_outliers)) == tolower(org)]
    if (!length(hit))
        stop("No dds found for organ: ", org)
    dds_no_outliers[[hit[[1]]]]
}
strip_ensembl_version <- function(x) sub("\\..*$", "", as.character(x))
add_gene_symbols <- function(res_df) {
    if (!all(c("gene", "log2FoldChange", "padj") %in% names(res_df))) {
        stop("DESeq2 result table must contain gene, log2FoldChange and padj.")
    }
    res_df <- dplyr::mutate(res_df, gene_base = strip_ensembl_version(gene))
    symbol_map <- AnnotationDbi::mapIds(org.Mm.eg.db, keys = unique(res_df$gene_base), keytype = "ENSEMBL",
        column = "SYMBOL", multiVals = "first")
    dplyr::mutate(res_df, geneSymbol = unname(symbol_map[gene_base]), gene_label = dplyr::if_else(!is.na(geneSymbol) &
        geneSymbol != "", geneSymbol, gene_base))
}
read_de <- function(org) {
    dplyr::mutate(add_gene_symbols(readr::read_csv(project_input(de_files[[org]]), show_col_types = FALSE)),
        organ = org)
}
de_results <- setNames(lapply(organs, read_de), organs)
liver_de <- de_results[["Liver"]]
tested_liver <- if ("pvalue" %in% names(liver_de)) {
    dplyr::filter(liver_de, !is.na(pvalue))
} else {
    liver_de
}
entrez_map <- AnnotationDbi::mapIds(org.Mm.eg.db, keys = unique(tested_liver$gene_base), keytype = "ENSEMBL",
    column = "ENTREZID", multiVals = "first")
universe_entrez <- unique(stats::na.omit(unname(entrez_map)))
run_directional_go <- function(de_df, direction_label) {
    selected <- dplyr::filter(de_df, !is.na(padj), padj < ALPHA)
    if (direction_label == "Higher in PPI") {
        selected <- dplyr::filter(selected, log2FoldChange > 0)
    }
    else {
        selected <- dplyr::filter(selected, log2FoldChange < 0)
    }
    selected_entrez <- unique(stats::na.omit(unname(entrez_map[selected$gene_base])))
    if (!length(selected_entrez)) {
        stop("No mapped significant genes for GO direction: ", direction_label)
    }
    ego <- clusterProfiler::enrichGO(gene = selected_entrez, universe = universe_entrez, OrgDb = org.Mm.eg.db,
        keyType = "ENTREZID", ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1, readable = TRUE)
    dplyr::mutate(dplyr::filter(tibble::as_tibble(as.data.frame(ego)), !is.na(p.adjust), p.adjust < ALPHA),
        direction = direction_label)
}
go_all <- dplyr::bind_rows(run_directional_go(liver_de, "Higher in PPI"), run_directional_go(liver_de,
    "Lower in PPI"))
if (!nrow(go_all)) {
    stop("No liver GO Biological Process terms passed FDR < ", ALPHA)
}
go_plot_df <- dplyr::arrange(dplyr::mutate(dplyr::ungroup(dplyr::slice_min(dplyr::group_by(go_all, direction),
    order_by = p.adjust, n = TOP_GO_PER_DIRECTION, with_ties = FALSE)), direction = factor(direction,
    levels = c("Higher in PPI", "Lower in PPI")), neg_log10_fdr = -log10(pmax(p.adjust, .Machine$double.xmin)),
    term_wrapped = stringr::str_wrap(Description, width = 48), term_key = paste(term_wrapped, direction,
        sep = "___")), direction, dplyr::desc(neg_log10_fdr))
go_plot_df$term_key <- factor(go_plot_df$term_key, levels = rev(unique(go_plot_df$term_key)))
p_go <- ggplot2::ggplot(go_plot_df, ggplot2::aes(x = neg_log10_fdr, y = term_key, size = Count, colour = direction)) +
    ggplot2::geom_vline(xintercept = -log10(ALPHA), linetype = "dashed", colour = "grey50", linewidth = 0.35) +
    ggplot2::geom_point(alpha = 0.88) + ggplot2::facet_grid(rows = ggplot2::vars(direction), scales = "free_y",
    space = "free_y") + ggplot2::scale_y_discrete(labels = function(x) sub("___.*$", "", x)) + ggplot2::scale_colour_manual(values = direction_cols,
    guide = "none") + ggplot2::scale_size_continuous(range = c(2.5, 6), name = "Gene count") + ggplot2::labs(x = expression(-log[10]("FDR")),
    y = NULL, title = "Liver GO biological processes") + base_theme + ggplot2::theme(axis.text.y = ggplot2::element_text(size = 13),
    strip.background = ggplot2::element_blank(), strip.text.y = ggplot2::element_text(size = 15, face = "bold"),
    legend.position = "bottom")
significant_genes <- function(org) {
    dplyr::arrange(dplyr::filter(de_results[[org]], !is.na(padj), padj < ALPHA), padj)
}
colon_sig <- significant_genes("Colon")
brain_sig <- significant_genes("Brain")
liver_sig <- significant_genes("Liver")
if (!nrow(colon_sig)) stop("No significant colon genes at FDR < ", ALPHA)
if (!nrow(brain_sig)) stop("No significant brain genes at FDR < ", ALPHA)
lcn2_brain <- dplyr::filter(brain_sig, geneSymbol == "Lcn2")
lcn2_liver <- dplyr::filter(liver_sig, geneSymbol == "Lcn2")
if (nrow(lcn2_brain) != 1L || nrow(lcn2_liver) != 1L) {
    stop("Lcn2 must map uniquely and pass FDR < ", ALPHA, " in both brain and liver before it is shown as a cross-tissue result.")
}
brain_other_sig <- dplyr::filter(brain_sig, geneSymbol != "Lcn2")
if (!nrow(brain_other_sig)) {
    stop("No significant brain genes remain after moving Lcn2 to panel D.")
}
extract_expression <- function(org, selected_df) {
    dds <- get_dds(org)
    norm_counts <- DESeq2::counts(dds, normalized = TRUE)
    count_key <- dplyr::distinct(tibble::tibble(count_row = rownames(norm_counts), gene_base = strip_ensembl_version(rownames(norm_counts))),
        gene_base, .keep_all = TRUE)
    selected_df <- dplyr::mutate(dplyr::distinct(selected_df, gene_base, .keep_all = TRUE), gene_order = dplyr::row_number())
    selected_rows <- dplyr::left_join(selected_df, count_key, by = "gene_base")
    if (anyNA(selected_rows$count_row)) {
        stop("Could not find selected ", org, " gene(s) in the count matrix: ", paste(selected_rows$gene_base[is.na(selected_rows$count_row)],
            collapse = ", "))
    }
    count_df <- tidyr::pivot_longer(tibble::rownames_to_column(as.data.frame(t(norm_counts[selected_rows$count_row,
        , drop = FALSE])), "sample"), cols = -sample, names_to = "count_row", values_to = "normalized_count")
    sample_meta <- dplyr::transmute(tibble::rownames_to_column(as.data.frame(SummarizedExperiment::colData(dds)),
        "sample"), sample, group = factor(as.character(group), levels = c("Ctrl", "PPI")))
    dplyr::mutate(dplyr::left_join(dplyr::left_join(count_df, dplyr::select(selected_rows, count_row,
        gene_base, gene_label, geneSymbol, log2FoldChange, padj, gene_order), by = "count_row"), sample_meta,
        by = "sample"), organ = org, log2_normalized_count = log2(normalized_count + 1), gene_label = factor(gene_label,
        levels = selected_rows$gene_label[order(selected_rows$gene_order)]))
}
colon_expression <- extract_expression("Colon", colon_sig)
brain_expression <- extract_expression("Brain", brain_other_sig)
lcn2_expression <- dplyr::mutate(dplyr::bind_rows(extract_expression("Brain", lcn2_brain), extract_expression("Liver",
    lcn2_liver)), organ = factor(organ, levels = c("Brain", "Liver")))
expression_facet_plot <- function(plot_df, plot_title) {
    ggplot2::ggplot(plot_df, ggplot2::aes(x = group, y = log2_normalized_count, colour = group)) + ggplot2::geom_boxplot(ggplot2::aes(fill = group),
        width = 0.56, outlier.shape = NA, alpha = 0.16, linewidth = 0.65) + ggplot2::geom_point(position = ggplot2::position_jitter(width = 0.1,
        height = 0), size = 2.5, alpha = 0.9) + ggplot2::facet_wrap(ggplot2::vars(gene_label), nrow = 1,
        scales = "free_y") + ggplot2::scale_colour_manual(name = "Group", values = group_cols) + ggplot2::scale_fill_manual(values = group_cols,
        guide = "none") + ggplot2::labs(x = NULL, y = expression(log[2]("normalized count" + 1)), title = plot_title) +
        base_theme + ggplot2::theme(strip.background = ggplot2::element_blank(), strip.text = ggplot2::element_text(size = 15,
        face = "italic"), legend.position = "bottom")
}
p_colon <- expression_facet_plot(colon_expression, "Colon")
p_brain <- expression_facet_plot(brain_expression, "Brain")
p_lcn2 <- ggplot2::ggplot(lcn2_expression, ggplot2::aes(x = organ, y = log2_normalized_count, colour = group,
    fill = group, group = interaction(organ, group))) + ggplot2::geom_boxplot(position = ggplot2::position_dodge(width = 0.72),
    width = 0.58, outlier.shape = NA, alpha = 0.16, linewidth = 0.65) + ggplot2::geom_point(position = ggplot2::position_jitterdodge(jitter.width = 0.1,
    jitter.height = 0, dodge.width = 0.72), size = 2.5, alpha = 0.9) + ggplot2::scale_colour_manual(name = "Group",
    values = group_cols) + ggplot2::scale_fill_manual(values = group_cols, guide = "none") + ggplot2::labs(x = NULL,
    y = expression(log[2]("normalized count" + 1)), title = expression(italic("Lcn2"))) + base_theme +
    ggplot2::theme(legend.position = "bottom")
expression_legend <- cowplot::get_legend(p_colon + ggplot2::theme(legend.position = "bottom"))
expression_row <- cowplot::plot_grid(p_colon + ggplot2::theme(legend.position = "none"), p_brain + ggplot2::theme(legend.position = "none"),
    p_lcn2 + ggplot2::theme(legend.position = "none"), ncol = 3, rel_widths = c(1.12, 1.12, 0.76), labels = c("B",
        "C", "D"), label_size = 24, label_fontface = "bold")
expression_row_with_legend <- cowplot::plot_grid(expression_row, expression_legend, ncol = 1, rel_heights = c(1,
    0.1))
final_figure <- cowplot::plot_grid(p_go, expression_row_with_legend, ncol = 1, rel_heights = c(1.55,
    0.95), labels = c("A", ""), label_size = 24, label_fontface = "bold")
print(final_figure)
readr::write_csv(dplyr::select(dplyr::mutate(go_all, display_key = paste(ID, direction, sep = "___"),
    displayed = display_key %in% paste(go_plot_df$ID, go_plot_df$direction, sep = "___")), -display_key),
    project_output(file.path(DATA_DIR, "Figure_Sx_liver_GO_source_data.csv"), name = "figure_sx_liver_go_source_data.csv"))
readr::write_csv(dplyr::mutate(dplyr::bind_rows(colon_expression, brain_expression, lcn2_expression),
    gene_label = as.character(gene_label)), project_output(file.path(DATA_DIR, "Figure_Sx_gene_expression_source_data.csv"),
    name = "figure_sx_gene_expression_source_data.csv"))
pdf_file <- file.path(OUT_DIR, "Mouse_Transcriptomics_Supplementary_GO_and_gene_expression.pdf")
png_file <- file.path(OUT_DIR, "Mouse_Transcriptomics_Supplementary_GO_and_gene_expression.png")
grDevices::cairo_pdf(filename = project_output(pdf_file), width = FIG_W_IN, height = FIG_H_IN)
print(final_figure)
grDevices::dev.off()
if (requireNamespace("ragg", quietly = TRUE)) {
    ragg::agg_png(project_output(png_file), width = FIG_W_IN, height = FIG_H_IN, units = "in", res = DPI,
        background = "white")
    print(final_figure)
    grDevices::dev.off()
} else {
    grDevices::png(project_output(png_file), width = FIG_W_IN, height = FIG_H_IN, units = "in", res = DPI,
        bg = "white")
    print(final_figure)
    grDevices::dev.off()
}
message("Supplementary transcriptomics figure written to: ", OUT_DIR)
}))
