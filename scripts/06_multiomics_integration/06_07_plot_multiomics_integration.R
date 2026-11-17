# SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0
# Figure 6: cross-layer associations.

source("scripts/helpers/setup.R")
invisible(run_analysis(
  section = "06_multiomics_integration",
  directory = "figures/scripts",
  script = "06_07_plot_multiomics_integration.R",
  restricted = FALSE,
  code = {
required_packages <- c("readr", "dplyr", "tidyr", "stringr", "forcats", "ggplot2", "cowplot", "tibble",
    "grid")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) {
    stop("Install required packages first: ", paste(missing_packages, collapse = ", "))
}
project_check_plot_versions(required_packages)
FIG_ROOT <- project_path("figures")
DATA_DIR <- file.path(FIG_ROOT, "source_data", "CrossLayer_Integration")
OUT_DIR <- file.path(FIG_ROOT, "final")
dir.create(project_output(DATA_DIR, directory = TRUE), recursive = TRUE, showWarnings = FALSE)
dir.create(project_output(OUT_DIR, directory = TRUE), recursive = TRUE, showWarnings = FALSE)
REFRESH_SOURCE_DATA <- FALSE
PPI_ROOT <- project_path("")
CROSS_ROOT <- file.path(PPI_ROOT, "integration")
LEGACY_FIG_DIR <- file.path(CROSS_ROOT, "06_cross_layer_massive_figure_focus_concordance_colon_liver_v3_no_titles")
GENE_DIR <- file.path(CROSS_ROOT, "05_gene_metabolite_DESeq2_Group_Metabolite_Hallmark_Reactome", "complete_gene_results")
stage_input <- function(filename, candidates) {
    target <- file.path(DATA_DIR, filename)
    if (REFRESH_SOURCE_DATA || !project_file_exists(target)) {
        hits <- candidates[project_file_exists(candidates)]
        if (!length(hits)) {
            stop("Could not stage ", filename, ".\nChecked:\n", paste(candidates, collapse = "\n"))
        }
        if (!file.copy(project_input(hits[[1]]), project_output(target), overwrite = TRUE))
            stop("Failed to copy ", hits[[1]])
        message("Staged: ", filename)
    }
    target
}
panel_A_file <- stage_input("panel_A_focused_effect_data.csv", c(file.path(LEGACY_FIG_DIR, "panel_A_focused_effect_data.csv")))
panel_C_file <- stage_input("panel_C_colon_DESeq2_GSEA_data.csv", c(file.path(LEGACY_FIG_DIR, "panel_C_colon_DESeq2_GSEA_data.csv")))
panel_D_file <- stage_input("panel_D_liver_DESeq2_GSEA_data.csv", c(file.path(LEGACY_FIG_DIR, "panel_D_liver_DESeq2_GSEA_data.csv")))
gene_specs <- dplyr::ungroup(dplyr::mutate(dplyr::rowwise(tibble::tribble(~dataset_id, ~Organ, ~filename,
    "serum_lcms", "Liver", "gene_results__serum_lcms__Liver.csv.gz", "serum_lcms", "Colon", "gene_results__serum_lcms__Colon.csv.gz",
    "cecum_lcms", "Liver", "gene_results__cecum_lcms__Liver.csv.gz", "cecum_lcms", "Colon", "gene_results__cecum_lcms__Colon.csv.gz",
    "cecum_nmr", "Liver", "gene_results__cecum_nmr__Liver.csv.gz", "cecum_nmr", "Colon", "gene_results__cecum_nmr__Colon.csv.gz")),
    path = stage_input(filename, c(file.path(GENE_DIR, filename)))))
readr::write_csv(data.frame(file = c(basename(panel_A_file), basename(panel_C_file), basename(panel_D_file),
    gene_specs$filename), centralized_path = c(panel_A_file, panel_C_file, panel_D_file, gene_specs$path)),
    project_output(file.path(DATA_DIR, "CrossLayer_Integration_source_manifest.csv"), name = "cross_layer_integration_source_manifest.csv"))
BASE_SIZE <- 18
FIG_W_IN <- 18
FIG_H_IN <- 15.2
DPI <- 600
model_label <- "counts ~ Group + Metabolite"
theme_levels <- c("Aromatic amino-acid metabolism", "Lipid metabolism", "Bile acid and haem metabolism")
theme_palette <- c(`Aromatic amino-acid metabolism` = "#4C78A8", `Lipid metabolism` = "#F58518", `Bile acid and haem metabolism` = "#54A24B")
manuscript_theme <- ggplot2::theme_bw(base_size = BASE_SIZE, base_family = "Arial") + ggplot2::theme(panel.grid = ggplot2::element_blank(),
    panel.border = ggplot2::element_rect(colour = "#3A3A3A", fill = NA, linewidth = 0.35), axis.line = ggplot2::element_line(colour = "#3A3A3A",
        linewidth = 0.3), axis.ticks = ggplot2::element_line(colour = "#3A3A3A", linewidth = 0.3), axis.text = ggplot2::element_text(size = 16,
        colour = "#2F2F2F"), axis.title = ggplot2::element_text(size = 18, colour = "#2F2F2F"), strip.background = ggplot2::element_rect(fill = "#D9D9D9",
        colour = "#3A3A3A", linewidth = 0.35), strip.text = ggplot2::element_text(face = "bold", size = 16,
        colour = "#2F2F2F"), legend.title = ggplot2::element_text(face = "bold", size = 16), legend.text = ggplot2::element_text(size = 15),
    plot.title = ggplot2::element_blank(), plot.subtitle = ggplot2::element_blank(), plot.margin = ggplot2::margin(7,
        7, 7, 7))
panel_A_labeller <- ggplot2::as_labeller(c(`Aromatic amino-acid metabolism` = "Aromatic amino-acid\nmetabolism",
    `Lipid metabolism` = "Lipid\nmetabolism", `Bile acid and haem metabolism` = "Bile acid and\nhaem metabolism"))
panel_A_data <- readr::read_csv(project_input(panel_A_file), show_col_types = FALSE)
if (!"term" %in% names(panel_A_data)) {
    panel_A_data$term <- paste0(panel_A_data$genus, " - ", panel_A_data$metabolite_short)
}
panel_A_data <- dplyr::arrange(dplyr::mutate(panel_A_data, theme = factor(theme, levels = theme_levels),
    term_key = paste0(theme_order, "__", plot_order, "__", term)), theme_order, plot_order)
panel_A_levels <- rev(panel_A_data$term_key)
panel_A_data <- dplyr::mutate(panel_A_data, term_key = factor(term_key, levels = panel_A_levels))
panel_A <- ggplot2::ggplot(panel_A_data, ggplot2::aes(beta, term_key)) + ggplot2::geom_vline(xintercept = 0,
    linetype = "dashed", linewidth = 0.35, colour = "#6F6F6F") + ggplot2::geom_segment(ggplot2::aes(x = conf_low,
    xend = conf_high, yend = term_key), linewidth = 0.45, colour = "#555555") + ggplot2::geom_point(ggplot2::aes(colour = theme),
    size = 2.8) + ggplot2::facet_grid(rows = ggplot2::vars(theme), scales = "free_y", space = "free_y",
    labeller = ggplot2::labeller(theme = panel_A_labeller)) + ggplot2::scale_colour_manual(values = theme_palette,
    guide = "none") + ggplot2::scale_y_discrete(labels = function(x) sub("^[0-9]+__[0-9]+__", "", x)) +
    ggplot2::labs(x = "Association coefficient (beta)", y = NULL) + manuscript_theme + ggplot2::theme(axis.text.y = ggplot2::element_text(size = 15.5),
    strip.text.y = ggplot2::element_text(angle = 0, size = 13.5))
hpp_names <- c(serum_lcms = "3-(3-Hydroxyphenyl)propanoic acid", cecum_lcms = "3-(3-Hydroxyphenyl)propionic acid",
    cecum_nmr = "3-hydroxyphenylpropionate")
read_hpp_gene_statistics <- function(dataset_id, organ, path) {
    target_metabolite <- unname(hpp_names[[dataset_id]])
    dat <- dplyr::select(dplyr::filter(readr::read_csv(project_input(path), col_types = readr::cols_only(metabolite = readr::col_character(),
        model = readr::col_character(), gene = readr::col_character(), statistic = readr::col_double()),
        progress = FALSE), model == model_label, metabolite == target_metabolite, !is.na(gene), is.finite(statistic)),
        gene, statistic)
    if (!nrow(dat))
        stop("No 3-HPP rows for ", dataset_id, " / ", organ)
    names(dat)[names(dat) == "statistic"] <- paste0("stat_", dataset_id)
    dat
}
gene_stats <- dplyr::ungroup(dplyr::mutate(dplyr::rowwise(gene_specs), hpp_data = list(read_hpp_gene_statistics(dataset_id,
    Organ, path))))
get_hpp <- function(dataset_id, organ) {
    idx <- which(gene_stats$dataset_id == dataset_id & gene_stats$Organ == organ)
    if (length(idx) != 1L)
        stop("Could not identify ", dataset_id, " / ", organ)
    gene_stats$hpp_data[[idx]]
}
comparison_specs <- tibble::tribble(~Organ, ~comparison_dataset, ~comparison_label, "Liver", "cecum_lcms",
    "Serum vs caecal LC-MS", "Liver", "cecum_nmr", "Serum LC-MS vs caecal NMR", "Colon", "cecum_lcms",
    "Serum vs caecal LC-MS", "Colon", "cecum_nmr", "Serum LC-MS vs caecal NMR")
build_concordance <- function(organ, comparison_dataset, comparison_label) {
    serum <- get_hpp("serum_lcms", organ)
    comp <- get_hpp(comparison_dataset, organ)
    merged <- dplyr::inner_join(serum, comp, by = "gene")
    comp_col <- paste0("stat_", comparison_dataset)
    rho <- stats::cor(merged$stat_serum_lcms, merged[[comp_col]], method = "spearman", use = "complete.obs")
    tibble::tibble(gene = merged$gene, serum_statistic = merged$stat_serum_lcms, comparison_statistic = merged[[comp_col]],
        Organ = organ, comparison_label = comparison_label, rho = rho)
}
panel_B_data <- dplyr::bind_rows(Map(build_concordance, comparison_specs$Organ, comparison_specs$comparison_dataset,
    comparison_specs$comparison_label))
panel_B_labels <- dplyr::mutate(dplyr::distinct(panel_B_data, Organ, comparison_label, rho), label = paste0("ρ = ",
    sprintf("%.3f", rho)))
panel_B <- ggplot2::ggplot(panel_B_data, ggplot2::aes(serum_statistic, comparison_statistic)) + ggplot2::geom_abline(slope = 1,
    intercept = 0, linetype = "dashed", linewidth = 0.35, colour = "#6F6F6F") + ggplot2::geom_point(size = 0.75,
    alpha = 0.16, colour = "#333333") + ggplot2::geom_text(data = panel_B_labels, ggplot2::aes(x = -Inf,
    y = Inf, label = label), inherit.aes = FALSE, hjust = -0.08, vjust = 1.2, size = 5) + ggplot2::facet_grid(rows = ggplot2::vars(Organ),
    cols = ggplot2::vars(comparison_label), scales = "free") + ggplot2::labs(x = "Serum 3-HPP Wald statistic",
    y = "Comparison-layer 3-HPP Wald statistic") + manuscript_theme + ggplot2::theme(axis.text = ggplot2::element_text(size = 15.5),
    axis.title = ggplot2::element_text(size = 17), strip.text = ggplot2::element_text(size = 15.5), panel.spacing = grid::unit(1.2,
        "mm"))
profile_levels <- c("3-HPP\nserum", "3-HPP\ncaecal LC-MS", "3-HPP\ncaecal NMR", "3-PPA\nserum", "N-cinnamoyl-\nglycine",
    "4-HPP\nserum", "FA 22:4\nserum", "Stercobilin\ncaecum", "Mesobilirubinogen\ncaecum", "Urobilinogen/\nurobilin")
pathway_levels <- c("E2F targets", "G2M checkpoint", "MYC targets", "Mitotic spindle", "mTORC1 signalling",
    "Protein secretion", "Cholesterol homeostasis", "Cholesterol biosynthesis", "Interferon-alpha response",
    "Interferon-gamma response", "TNF-NF-kB signalling", "Inflammatory response", "IL6-JAK-STAT3 signalling",
    "Oxidative phosphorylation")
prepare_gsea_grid <- function(path) {
    x <- readr::read_csv(project_input(path), show_col_types = FALSE)
    dplyr::mutate(x, profile_label = factor(as.character(profile_label), levels = profile_levels), pathway_label = factor(as.character(pathway_label),
        levels = rev(pathway_levels)), NES_plot = as.numeric(NES_plot))
}
colon_grid <- prepare_gsea_grid(panel_C_file)
liver_grid <- prepare_gsea_grid(panel_D_file)
make_gsea_panel <- function(grid) {
    ggplot2::ggplot(grid, ggplot2::aes(profile_label, pathway_label, fill = NES_plot)) + ggplot2::geom_tile(colour = "white",
        linewidth = 0.35) + ggplot2::scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
        midpoint = 0, na.value = "#F2F2F2", name = "NES") + ggplot2::labs(x = NULL, y = NULL) + manuscript_theme +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 40, hjust = 1, vjust = 1, size = 15),
            axis.text.y = ggplot2::element_text(size = 15.5), legend.position = "bottom", legend.box.spacing = grid::unit(0,
                "mm"), legend.margin = ggplot2::margin(t = -5, r = 0, b = 0, l = 0), plot.margin = ggplot2::margin(t = 5,
                r = 7, b = 1, l = 7))
}
panel_C <- make_gsea_panel(colon_grid)
panel_D <- make_gsea_panel(liver_grid)
row_1 <- cowplot::plot_grid(panel_A, panel_B, nrow = 1, rel_widths = c(1.25, 1.25), labels = c("A", "B"),
    label_size = 24, label_fontface = "bold", label_x = 0, label_y = 1, hjust = 0, vjust = 1)
row_2 <- cowplot::plot_grid(panel_C, nrow = 1, labels = "C", label_size = 24, label_fontface = "bold",
    label_x = 0, label_y = 1, hjust = 0, vjust = 1)
row_3 <- cowplot::plot_grid(panel_D, nrow = 1, labels = "D", label_size = 24, label_fontface = "bold",
    label_x = 0, label_y = 1, hjust = 0, vjust = 1)
final_figure <- cowplot::plot_grid(row_1, row_2, row_3, ncol = 1, rel_heights = c(1.08, 1.13, 1.13))
print(final_figure)
pdf_file <- file.path(OUT_DIR, "CrossLayer_Integration.pdf")
png_file <- file.path(OUT_DIR, "CrossLayer_Integration.png")
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
message("Figure 6 written to: ", OUT_DIR)
}))
