# SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0
# Figure 4: experimental metabolomics.

source("scripts/helpers/setup.R")
invisible(run_analysis(
  section = "04_mouse_metabolomics",
  directory = "figures/scripts",
  script = "04_13_plot_mouse_metabolomics.R",
  restricted = FALSE,
  code = {
required_packages <- c("dplyr", "tidyr", "stringr", "readr", "ggplot2", "cowplot")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) {
    stop("Install required packages first: ", paste(missing_packages, collapse = ", "))
}
project_check_plot_versions(required_packages)
FIG_ROOT <- project_path("figures")
DATA_DIR <- file.path(FIG_ROOT, "source_data", "Experimental_Metabolomics")
OUT_DIR <- file.path(FIG_ROOT, "final")
dir.create(project_output(DATA_DIR, directory = TRUE), recursive = TRUE, showWarnings = FALSE)
dir.create(project_output(OUT_DIR, directory = TRUE), recursive = TRUE, showWarnings = FALSE)
REFRESH_SOURCE_DATA <- FALSE
PPI_ROOT <- project_path("")
LEGACY_DIR <- file.path(PPI_ROOT, "metabolomics", "Figure3_remade_mouse_dupCor_combined_significance")
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
panel_files <- character()
for (letter in LETTERS[1:5]) {
    panel_files <- c(panel_files, stage_input(paste0("panel_", letter, "_metadata.csv"), c(file.path(LEGACY_DIR,
        paste0("panel_", letter, "_metadata.csv")))), stage_input(paste0("panel_", letter, "_plot_matrix.csv"),
        c(file.path(LEGACY_DIR, paste0("panel_", letter, "_plot_matrix.csv")))))
}
heatmap_F_file <- stage_input("panel_F_heatmap_long.csv", c(file.path(LEGACY_DIR, "panel_F_heatmap_long.csv")))
heatmap_G_file <- stage_input("panel_G_heatmap_long.csv", c(file.path(LEGACY_DIR, "panel_G_heatmap_long.csv")))
heatmap_H_file <- stage_input("panel_H_heatmap_long.csv", c(file.path(LEGACY_DIR, "panel_H_heatmap_long.csv")))
readr::write_csv(data.frame(file = basename(c(panel_files, heatmap_F_file, heatmap_G_file, heatmap_H_file)),
    centralized_path = c(panel_files, heatmap_F_file, heatmap_G_file, heatmap_H_file)), project_output(file.path(DATA_DIR,
    "Experimental_Metabolomics_source_manifest.csv"), name = "experimental_metabolomics_source_manifest.csv"))
BASE_SIZE <- 18
FIG_W_IN <- 18
FIG_H_IN <- 17
DPI <- 600
base_theme <- ggplot2::theme_bw(base_size = BASE_SIZE, base_family = "Arial") + ggplot2::theme(panel.grid = ggplot2::element_line(colour = "#E6E6E6",
    linewidth = 0.25), panel.border = ggplot2::element_rect(colour = "#4A4A4A", fill = NA, linewidth = 0.4),
    axis.ticks = ggplot2::element_line(colour = "#4A4A4A", linewidth = 0.3), axis.line = ggplot2::element_line(colour = "#4A4A4A",
        linewidth = 0.3), axis.text = ggplot2::element_text(size = 16, colour = "#2F2F2F"), axis.title = ggplot2::element_text(size = 18,
        colour = "#2F2F2F"), strip.text = ggplot2::element_text(size = 16), legend.title = ggplot2::element_text(size = 16),
    legend.text = ggplot2::element_text(size = 15), legend.background = ggplot2::element_blank(), legend.key = ggplot2::element_blank(),
    plot.margin = ggplot2::margin(7, 7, 7, 7))
load_panel_source <- function(letter) {
    meta <- readr::read_csv(project_input(file.path(DATA_DIR, paste0("panel_", letter, "_metadata.csv"))),
        show_col_types = FALSE)
    both <- readr::read_csv(project_input(file.path(DATA_DIR, paste0("panel_", letter, "_plot_matrix.csv"))),
        show_col_types = FALSE)
    feature_cols <- setdiff(names(both), names(meta))
    if (length(feature_cols) < 2) {
        stop("Could not identify >=2 feature columns for panel ", letter)
    }
    X <- as.data.frame(both[, feature_cols, drop = FALSE], check.names = FALSE)
    list(metadata = meta, matrix = X)
}
median_impute_cols <- function(df) {
    df <- as.data.frame(df, check.names = FALSE)
    for (j in seq_along(df)) {
        x <- as.numeric(df[[j]])
        med <- stats::median(x, na.rm = TRUE)
        if (!is.finite(med))
            med <- 0
        x[is.na(x)] <- med
        df[[j]] <- x
    }
    df
}
scale_features_for_pca <- function(df) {
    X <- median_impute_cols(df)
    keep <- vapply(X, function(x) stats::sd(x, na.rm = TRUE) > 0, logical(1))
    X <- X[, keep, drop = FALSE]
    if (ncol(X) < 2)
        stop("Fewer than two variable features remain for PCA.")
    scale(X)
}
find_donor_col <- function(meta) {
    hits <- c("donor", "Donor", "inoculum_donor", "Inoculum_Donor")
    hit <- hits[hits %in% names(meta)]
    if (!length(hit))
        stop("Could not identify donor column in metadata: ", paste(names(meta), collapse = ", "))
    hit[[1]]
}
find_group_col <- function(meta) {
    hits <- c("Group", "group", "Separation")
    hit <- hits[hits %in% names(meta)]
    if (!length(hit))
        stop("Could not identify group column in metadata.")
    hit[[1]]
}
normalize_group <- function(x) {
    z <- tolower(as.character(x))
    ifelse(z %in% c("ctrl", "control", "0"), "Ctrl", "PPI")
}
euclidean_group_r2 <- function(X_scaled, group) {
    X_scaled <- as.matrix(X_scaled)
    group <- factor(group)
    grand <- colMeans(X_scaled)
    ss_total <- sum((X_scaled - matrix(grand, nrow(X_scaled), ncol(X_scaled), byrow = TRUE))^2)
    ss_within <- 0
    for (g in levels(group)) {
        Xi <- X_scaled[group == g, , drop = FALSE]
        centre <- colMeans(Xi)
        ss_within <- ss_within + sum((Xi - matrix(centre, nrow(Xi), ncol(Xi), byrow = TRUE))^2)
    }
    1 - ss_within/ss_total
}
exact_donor_label_permanova <- function(X_scaled, metadata) {
    donor_col <- find_donor_col(metadata)
    group_col <- find_group_col(metadata)
    donor <- as.character(metadata[[donor_col]])
    group <- normalize_group(metadata[[group_col]])
    donor_group <- unique(data.frame(donor = donor, group = group, stringsAsFactors = FALSE))
    if (nrow(donor_group) != 6L || sum(donor_group$group == "PPI") != 3L || sum(donor_group$group ==
        "Ctrl") != 3L) {
        stop("Exact donor-label test expects six donors with 3 PPI and 3 Ctrl. Found:\n", paste(capture.output(print(donor_group)),
            collapse = "\n"))
    }
    observed_r2 <- euclidean_group_r2(X_scaled, group)
    combos <- utils::combn(donor_group$donor, 3, simplify = FALSE)
    perm_r2 <- vapply(combos, function(ppi_donors) {
        perm_group <- ifelse(donor %in% ppi_donors, "PPI", "Ctrl")
        euclidean_group_r2(X_scaled, perm_group)
    }, numeric(1))
    p_exact <- mean(perm_r2 >= observed_r2 - 1e-12)
    list(r2 = observed_r2, p = p_exact, n_permutations = length(perm_r2))
}
make_pca_panel <- function(source, panel_label, shape_col = NULL, annotation_position = c("tr", "br",
    "tl", "bl")) {
    annotation_position <- match.arg(annotation_position)
    X_scaled <- scale_features_for_pca(source$matrix)
    pca <- stats::prcomp(X_scaled, center = FALSE, scale. = FALSE)
    var_exp <- 100 * pca$sdev^2/sum(pca$sdev^2)
    scores <- as.data.frame(pca$x[, 1:2, drop = FALSE])
    meta <- source$metadata
    group_col <- find_group_col(meta)
    meta$Group_plot <- normalize_group(meta[[group_col]])
    scores <- dplyr::bind_cols(meta, scores)
    perm <- exact_donor_label_permanova(X_scaled, meta)
    lab <- panel_label
    p <- ggplot2::ggplot(scores, ggplot2::aes(PC1, PC2, colour = Group_plot))
    if (!is.null(shape_col) && shape_col %in% names(scores)) {
        p <- p + ggplot2::geom_point(ggplot2::aes(shape = .data[[shape_col]]), size = 2.1)
    }
    else {
        p <- p + ggplot2::geom_point(size = 2.1)
    }
    p <- p + ggplot2::scale_colour_manual(values = c(Ctrl = "#E67E22", PPI = "#1F4E8C"), name = "Group")
    pos <- switch(annotation_position, tr = list(x = Inf, y = Inf, hjust = 1.04, vjust = 1.08), br = list(x = Inf,
        y = -Inf, hjust = 1.04, vjust = -0.08), tl = list(x = -Inf, y = Inf, hjust = -0.04, vjust = 1.08),
        bl = list(x = -Inf, y = -Inf, hjust = -0.04, vjust = -0.08))
    p + ggplot2::annotate("text", x = pos$x, y = pos$y, label = lab, hjust = pos$hjust, vjust = pos$vjust,
        size = 4.3) + ggplot2::labs(x = paste0("PC1 (", sprintf("%.1f", var_exp[1]), "%)"), y = paste0("PC2 (",
        sprintf("%.1f", var_exp[2]), "%)")) + base_theme + ggplot2::theme(legend.position = "bottom",
        legend.box = "vertical", legend.margin = ggplot2::margin(3, 8, 3, 8),
        legend.box.spacing = grid::unit(1, "mm")) + ggplot2::guides(
        colour = ggplot2::guide_legend(order = 1, nrow = 1),
        shape = ggplot2::guide_legend(order = 2, nrow = 1))
}
make_sig_stars <- function(fdr) {
    dplyr::case_when(is.na(fdr) ~ "", fdr < 0.01 ~ "***", fdr < 0.05 ~ "**", fdr < 0.1 ~ "*", TRUE ~
        "")
}
build_heatmap_panel <- function(long_df, dataset_order, top_n_total = NULL, y_title = NULL, metabolite_font_size = 15) {
    long_df <- dplyr::ungroup(dplyr::mutate(dplyr::group_by(dplyr::mutate(long_df, Dataset = as.character(Dataset),
        Metabolite = as.character(Metabolite), effect = as.numeric(effect), FDR = as.numeric(FDR)), Dataset),
        scaled_effect = {
            m <- max(abs(effect), na.rm = TRUE)
            if (is.finite(m) && m > 0)
                effect/m
            else effect
        }, sig_label = make_sig_stars(FDR)))
    if (!is.null(top_n_total)) {
        keep <- dplyr::pull(dplyr::slice_head(dplyr::arrange(dplyr::summarise(dplyr::group_by(long_df,
            Metabolite), max_abs = max(abs(effect), na.rm = TRUE), .groups = "drop"), dplyr::desc(max_abs),
            Metabolite), n = top_n_total), Metabolite)
        long_df <- dplyr::filter(long_df, Metabolite %in% keep)
    }
    ord <- dplyr::arrange(dplyr::summarise(dplyr::group_by(long_df, Metabolite), order_effect = mean(scaled_effect,
        na.rm = TRUE), .groups = "drop"), dplyr::desc(order_effect), Metabolite)
    levels_y <- rev(ord$Metabolite)
    grid <- dplyr::mutate(dplyr::left_join(tidyr::expand_grid(Metabolite = levels_y, Dataset = dataset_order),
        long_df, by = c("Metabolite", "Dataset")), Dataset = factor(Dataset, levels = dataset_order),
        Metabolite = factor(Metabolite, levels = levels_y))
    ggplot2::ggplot(grid, ggplot2::aes(Dataset, Metabolite, fill = scaled_effect)) + ggplot2::geom_tile(width = 0.78,
        height = 0.88, colour = NA) + ggplot2::geom_text(ggplot2::aes(label = sig_label), size = 4, na.rm = TRUE) +
        ggplot2::scale_fill_gradient2(low = "#3B4CC0", mid = "#F3F3F3", high = "#D7191C", midpoint = 0,
            limits = c(-1, 1), breaks = c(-1, -0.5, 0, 0.5, 1),
            labels = c("-1", "-0.5", "0", "0.5", "1"),
            na.value = "#D9D9D9", name = "Scaled effect") + ggplot2::labs(x = NULL,
        y = y_title) + base_theme + ggplot2::theme(panel.grid = ggplot2::element_blank(), panel.border = ggplot2::element_blank(),
        axis.line = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(), axis.text.x = ggplot2::element_text(angle = 45,
            hjust = 1, vjust = 1, size = 15), axis.text.y = ggplot2::element_text(size = metabolite_font_size),
        axis.title.y = ggplot2::element_text(size = 17), legend.position = "none")
}
A <- load_panel_source("A")
B <- load_panel_source("B")
C <- load_panel_source("C")
D <- load_panel_source("D")
E <- load_panel_source("E")
panel_A <- make_pca_panel(A, "NMR ex vivo", shape_col = if ("Day" %in% names(A$metadata)) "Day" else NULL,
    annotation_position = "br")
panel_B <- make_pca_panel(B, "HILIC ex vivo", shape_col = if ("Day" %in% names(B$metadata)) "Day" else NULL,
    annotation_position = "tr")
panel_C <- make_pca_panel(C, "NMR caecum", annotation_position = "br")
panel_D <- make_pca_panel(D, "LC–MS caecum/faeces", shape_col = if ("Organ" %in% names(D$metadata)) "Organ" else NULL,
    annotation_position = "tl")
panel_E <- make_pca_panel(E, "LC–MS serum", annotation_position = "br")
long_F <- readr::read_csv(project_input(heatmap_F_file), show_col_types = FALSE)
long_G <- readr::read_csv(project_input(heatmap_G_file), show_col_types = FALSE)
long_H <- readr::read_csv(project_input(heatmap_H_file), show_col_types = FALSE)
panel_F <- build_heatmap_panel(long_F, dataset_order = c("HILIC (Ex vivo)", "NMR (Ex vivo)"), y_title = "Metabolite",
    metabolite_font_size = 15.5)
panel_G <- build_heatmap_panel(long_G, dataset_order = c("LCMS (Caecum)", "LCMS (Faeces)", "NMR (Caecum)"),
    y_title = NULL, metabolite_font_size = 15)
panel_H <- build_heatmap_panel(long_H, dataset_order = c("LCMS (Serum)"), top_n_total = 30, y_title = NULL,
    metabolite_font_size = 15.5)
legend_plot <- panel_H + ggplot2::theme(legend.position = "bottom",
    legend.margin = ggplot2::margin(6, 12, 8, 12)) + ggplot2::guides(fill = ggplot2::guide_colourbar(
    direction = "horizontal", barwidth = grid::unit(80, "mm"), barheight = grid::unit(4, "mm"),
    title.position = "top", title.hjust = 0.5, label.position = "bottom"))
heatmap_legend <- cowplot::get_legend(legend_plot)
pca_top_left <- cowplot::plot_grid(panel_A, panel_C, nrow = 1, labels = c("A", "C"), label_size = 24,
    label_fontface = "bold", label_x = 0, label_y = 1, hjust = 0, vjust = 1)
pca_bottom_left <- cowplot::plot_grid(panel_B, panel_D, nrow = 1, labels = c("B", "D"), label_size = 24,
    label_fontface = "bold", label_x = 0, label_y = 1, hjust = 0, vjust = 1)
pca_left_block <- cowplot::plot_grid(pca_top_left, pca_bottom_left, ncol = 1)
pca_block <- cowplot::plot_grid(pca_left_block, panel_E, nrow = 1, rel_widths = c(2, 1), labels = c("",
    "E"), label_size = 24, label_fontface = "bold", label_x = 0, label_y = 1, hjust = 0, vjust = 1)
panel_H_compact <- cowplot::plot_grid(NULL, panel_H, NULL, nrow = 1, rel_widths = c(0.18, 0.64, 0.18))
heatmap_row <- cowplot::plot_grid(panel_F, panel_G, panel_H_compact, nrow = 1, labels = c("F", "G", "H"),
    label_size = 24, label_fontface = "bold", label_x = 0, label_y = 1, hjust = 0, vjust = 1)
final_figure <- cowplot::plot_grid(pca_block, heatmap_row, heatmap_legend, ncol = 1, rel_heights = c(1.7,
    1.35, 0.18))
print(final_figure)
pdf_file <- file.path(OUT_DIR, "Experimental_Metabolomics.pdf")
png_file <- file.path(OUT_DIR, "Experimental_Metabolomics.png")
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
message("Figure 4 written to: ", OUT_DIR)
}))
