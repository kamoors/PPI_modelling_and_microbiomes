# SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0
# Liver GO and selected expression
# Run through run.R so input and output paths resolve within the package.
# Inputs: see data/input_manifest.csv and docs/RUN_ORDER.md.
# Writes generated files under results/; bundled inputs are read-only.

required_packages <- c("dplyr", "tidyr", "tibble", "readr", "ggplot2", "cowplot", "stringr")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) {
    stop("Install required packages first: ", paste(missing_packages, collapse = ", "))
}
ppi_check_plot_versions(required_packages)
FIG_ROOT <- ppi_path("figures")
PPI_ROOT <- ppi_path("")
RNA_ROOT <- file.path(PPI_ROOT, "rnaseq")
DATA_DIR <- file.path(FIG_ROOT, "source_data", "Mouse_Transcriptomics_Supplementary")
OUT_DIR <- file.path(FIG_ROOT, "final")
dir.create(ppi_output(DATA_DIR, directory = TRUE), recursive = TRUE, showWarnings = FALSE)
dir.create(ppi_output(OUT_DIR, directory = TRUE), recursive = TRUE, showWarnings = FALSE)
# These are published plotting tables. GO refitting is retained separately.
go_file <- file.path(DATA_DIR, "Figure_Sx_liver_GO_source_data.csv")
expression_file <- file.path(DATA_DIR, "Figure_Sx_gene_expression_source_data.csv")
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
go_all <- readr::read_csv(ppi_input(go_file), show_col_types = FALSE)
expression <- readr::read_csv(ppi_input(expression_file), show_col_types = FALSE)
expression$group <- factor(expression$group, levels = c("Ctrl", "PPI"))
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
ordered_expression <- function(df) {
    labels <- unique(df$gene_label[order(df$gene_order)])
    df$gene_label <- factor(df$gene_label, levels = labels)
    df
}
colon_expression <- ordered_expression(dplyr::filter(expression, organ == "Colon"))
brain_expression <- ordered_expression(dplyr::filter(expression, organ == "Brain", is.na(geneSymbol) | geneSymbol != "Lcn2"))
lcn2_expression <- dplyr::mutate(dplyr::filter(expression, geneSymbol == "Lcn2"),
    organ = factor(organ, levels = c("Brain", "Liver")))
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
    ppi_output(file.path(DATA_DIR, "Figure_Sx_liver_GO_source_data.csv"), name = "figure_sx_liver_go_source_data.csv"))
readr::write_csv(dplyr::mutate(dplyr::bind_rows(colon_expression, brain_expression, lcn2_expression), 
    gene_label = as.character(gene_label)), ppi_output(file.path(DATA_DIR, "Figure_Sx_gene_expression_source_data.csv"), 
    name = "figure_sx_gene_expression_source_data.csv"))
pdf_file <- file.path(OUT_DIR, "Mouse_Transcriptomics_Supplementary_GO_and_gene_expression.pdf")
png_file <- file.path(OUT_DIR, "Mouse_Transcriptomics_Supplementary_GO_and_gene_expression.png")
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
message("Supplementary transcriptomics figure written to: ", OUT_DIR)
