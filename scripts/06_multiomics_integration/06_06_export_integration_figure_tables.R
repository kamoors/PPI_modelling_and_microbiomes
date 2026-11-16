# SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0
# Export exact focused panel tables
# Run through run.R so input and output paths resolve within the package.
# Inputs: see data/input_manifest.csv and docs/RUN_ORDER.md.
# Writes generated files under results/; bundled inputs are read-only.

required_packages <- c("readr", "dplyr", "tidyr", "stringr", "forcats", "ggplot2", "cowplot", "tibble")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0L) {
    stop("Install the following packages before running this script: ", paste(missing_packages, collapse = ", "))
}
crosslayer_root <- ppi_path("integration")
pipeline02_dir <- file.path(crosslayer_root, "02_association_results")
pipeline05_dir <- file.path(crosslayer_root, "05_gene_metabolite_DESeq2_Group_Metabolite_Hallmark_Reactome")
gene_result_dir <- file.path(pipeline05_dir, "complete_gene_results")
output_dir <- file.path(crosslayer_root, "06_cross_layer_massive_figure_focus_concordance_colon_liver_v3_no_titles")
dir.create(ppi_output(output_dir, directory = TRUE), recursive = TRUE, showWarnings = FALSE)
find_preferred_file <- function(preferred, search_directory, pattern, recursive = TRUE) {
    if (ppi_file_exists(preferred)) {
        return(preferred)
    }
    candidates <- ppi_list_files(search_directory, pattern = pattern, recursive = recursive, full.names = TRUE)
    if (length(candidates) == 0L) {
        stop("No file matching '", pattern, "' was found under:\n", search_directory)
    }
    candidates <- candidates[order(nchar(basename(candidates)), basename(candidates))]
    warning("Preferred file not found. Using:\n", candidates[1])
    candidates[1]
}
association_file <- find_preferred_file(preferred = file.path(pipeline02_dir, "ALL_significant_associations_FDR0p1.csv"), 
    search_directory = pipeline02_dir, pattern = "^ALL_significant_associations_FDR0p1.*[.]csv$", recursive = TRUE)
gsea_file <- find_preferred_file(preferred = file.path(pipeline05_dir, "significant_gene_metabolite_GSEA_results.csv"), 
    search_directory = pipeline05_dir, pattern = "^significant_gene_metabolite_GSEA_results.*[.]csv$", 
    recursive = FALSE)
gene_files <- dplyr::mutate(tibble::tribble(~dataset_id, ~Organ, ~filename, "serum_lcms", "Liver", "gene_results__serum_lcms__Liver.csv.gz", 
    "serum_lcms", "Colon", "gene_results__serum_lcms__Colon.csv.gz", "cecum_lcms", "Liver", "gene_results__cecum_lcms__Liver.csv.gz", 
    "cecum_lcms", "Colon", "gene_results__cecum_lcms__Colon.csv.gz", "cecum_nmr", "Liver", "gene_results__cecum_nmr__Liver.csv.gz", 
    "cecum_nmr", "Colon", "gene_results__cecum_nmr__Colon.csv.gz"), path = file.path(gene_result_dir, 
    .data$filename))
missing_gene_files <- gene_files$path[!ppi_file_exists(gene_files$path)]
if (length(missing_gene_files) > 0L) {
    stop("The following DESeq2 gene-result files were not found:\n", paste(missing_gene_files, collapse = "\n"))
}
model_label <- "counts ~ Group + Metabolite"
family_fdr_cutoff <- 0.1
script_version <- "2026-07-21-v2-liver-filter-before-grid"
message("Running cross-layer figure script: ", script_version)
figure_width_mm <- 360
figure_height_mm <- 300
figure_dpi <- 400
theme_levels <- c("Aromatic amino-acid metabolism", "Lipid metabolism", "Bile acid and haem metabolism")
theme_palette <- c(`Aromatic amino-acid metabolism` = "#4C78A8", `Lipid metabolism` = "#F58518", `Bile acid and haem metabolism` = "#54A24B")
layer_labels <- c(serum_lcms = "serum LC-MS", cecum_lcms = "caecal LC-MS", cecum_nmr = "caecal NMR", 
    feces_lcms = "faecal LC-MS")
manuscript_theme <- ggplot2::theme_bw(base_size = 10) + ggplot2::theme(panel.grid = ggplot2::element_blank(), 
    panel.border = ggplot2::element_rect(colour = "#3A3A3A", fill = NA, linewidth = 0.4), axis.line = ggplot2::element_line(colour = "#3A3A3A", 
        linewidth = 0.35), axis.ticks = ggplot2::element_line(colour = "#3A3A3A", linewidth = 0.35), 
    axis.text = ggplot2::element_text(colour = "#2F2F2F"), axis.title = ggplot2::element_text(colour = "#2F2F2F"), 
    strip.background = ggplot2::element_rect(fill = "#D9D9D9", colour = "#3A3A3A", linewidth = 0.4), 
    strip.text = ggplot2::element_text(face = "bold", colour = "#2F2F2F"), legend.background = ggplot2::element_blank(), 
    legend.key = ggplot2::element_blank(), legend.title = ggplot2::element_text(face = "bold"), plot.title = ggplot2::element_blank(), 
    plot.subtitle = ggplot2::element_blank(), plot.margin = ggplot2::margin(5, 5, 5, 5))
association_results <- readr::read_csv(ppi_input(association_file), show_col_types = FALSE, progress = FALSE)
required_association_columns <- c("layer", "predictor_dataset", "predictor", "outcome", "adjusted_rank_beta", 
    "conf_low", "conf_high", "FDR")
missing_association_columns <- setdiff(required_association_columns, names(association_results))
if (length(missing_association_columns) > 0L) {
    stop("The association table is missing required columns: ", paste(missing_association_columns, collapse = ", "), 
        "\nAvailable columns:\n", paste(names(association_results), collapse = ", "))
}
genus_metabolite_all <- dplyr::transmute(dplyr::filter(association_results, .data$layer == "metabolite_genus"), 
    dataset_id = as.character(.data$predictor_dataset), metabolite = as.character(.data$predictor), genus = as.character(.data$outcome), 
    beta = as.numeric(.data$adjusted_rank_beta), conf_low = as.numeric(.data$conf_low), conf_high = as.numeric(.data$conf_high), 
    FDR = as.numeric(.data$FDR))
gsea_results <- readr::read_csv(ppi_input(gsea_file), show_col_types = FALSE, progress = FALSE)
required_gsea_columns <- c("dataset_id", "Organ", "metabolite", "model", "Library", "Pathway", "NES", 
    "GSEA_FDR_family")
missing_gsea_columns <- setdiff(required_gsea_columns, names(gsea_results))
if (length(missing_gsea_columns) > 0L) {
    stop("The DESeq2 GSEA table is missing required columns: ", paste(missing_gsea_columns, collapse = ", "), 
        "\nAvailable columns:\n", paste(names(gsea_results), collapse = ", "))
}
gsea_results <- dplyr::mutate(dplyr::filter(gsea_results, .data$model == model_label), dataset_id = as.character(.data$dataset_id), 
    Organ = as.character(.data$Organ), metabolite = as.character(.data$metabolite), Library = as.character(.data$Library), 
    Pathway = as.character(.data$Pathway), NES = as.numeric(.data$NES), GSEA_FDR_family = as.numeric(.data$GSEA_FDR_family))
selected_associations <- tibble::tribble(~theme, ~theme_order, ~plot_order, ~dataset_id, ~metabolite, 
    ~metabolite_short, ~genus, ~show_in_D, "Aromatic amino-acid metabolism", 1, 1, "cecum_lcms", "3-(3-Hydroxyphenyl)propionic acid", 
    "3-HPP", "Alistipes", TRUE, "Aromatic amino-acid metabolism", 1, 2, "cecum_nmr", "phenylalanine", 
    "Phenylalanine", "Alistipes", FALSE, "Aromatic amino-acid metabolism", 1, 3, "cecum_nmr", "tyrosine", 
    "Tyrosine", "Alistipes", TRUE, "Aromatic amino-acid metabolism", 1, 4, "cecum_nmr", "tryptophan", 
    "Tryptophan", "Alistipes", TRUE, "Aromatic amino-acid metabolism", 1, 5, "cecum_nmr", "phenylalanine", 
    "Phenylalanine", "Colidextribacter", FALSE, "Aromatic amino-acid metabolism", 1, 6, "cecum_nmr", 
    "propionate", "Propionate", "Colidextribacter", TRUE, "Lipid metabolism", 2, 1, "serum_lcms", "LPC 18:3", 
    "LPC 18:3", "Colidextribacter", TRUE, "Lipid metabolism", 2, 2, "serum_lcms", "LPI 20:4", "LPI 20:4", 
    "Colidextribacter", TRUE, "Lipid metabolism", 2, 3, "serum_lcms", "FA 17:0", "FA 17:0", "Colidextribacter", 
    TRUE, "Lipid metabolism", 2, 4, "serum_lcms", "FA 22:4", "FA 22:4", "Colidextribacter", TRUE, "Bile acid and haem metabolism", 
    3, 1, "cecum_lcms", "7-epiCA", "7-epiCA", "Agathobaculum", TRUE, "Bile acid and haem metabolism", 
    3, 2, "cecum_lcms", "CA-7S", "CA-7S", "Agathobaculum", TRUE, "Bile acid and haem metabolism", 3, 
    3, "cecum_lcms", "DCA-3-Sulfate", "DCA-3-sulfate", "Agathobaculum", TRUE, "Bile acid and haem metabolism", 
    3, 4, "serum_lcms", "DCA", "DCA", "Agathobaculum", TRUE)
selected_gm <- dplyr::left_join(selected_associations, genus_metabolite_all, by = c("dataset_id", "metabolite", 
    "genus"))
missing_selected_associations <- dplyr::filter(selected_gm, is.na(.data$beta))
if (nrow(missing_selected_associations) > 0L) {
    stop("The following curated genus-metabolite associations were not found:\n", paste(utils::capture.output(print(dplyr::select(missing_selected_associations, 
        .data$dataset_id, .data$metabolite, .data$genus))), collapse = "\n"))
}
selected_gm <- dplyr::mutate(selected_gm, theme = factor(.data$theme, levels = theme_levels), dataset_label = unname(layer_labels[.data$dataset_id]), 
    heatmap_row = paste0(.data$metabolite_short, " (", .data$dataset_label, ")"), detailed_label = paste0(.data$genus, 
        " - ", .data$metabolite_short, " (", .data$dataset_label, ")"))
profile_specs_base <- tibble::tribble(~profile_order, ~dataset_id, ~metabolite, ~profile_label, 1, "serum_lcms", 
    "3-(3-Hydroxyphenyl)propanoic acid", "3-HPP\nserum", 2, "cecum_lcms", "3-(3-Hydroxyphenyl)propionic acid", 
    "3-HPP\ncaecal LC-MS", 3, "cecum_nmr", "3-hydroxyphenylpropionate", "3-HPP\ncaecal NMR", 4, "serum_lcms", 
    "3-Phenylpropionic acid", "3-PPA\nserum", 5, "serum_lcms", "N-Cinnamoylglycine", "N-cinnamoyl-\nglycine", 
    6, "serum_lcms", "3-(4-Hydroxyphenyl)propionic acid", "4-HPP\nserum", 7, "serum_lcms", "FA 22:4", 
    "FA 22:4\nserum", 8, "cecum_lcms", "Stercobilin", "Stercobilin\ncaecum", 9, "cecum_lcms", "Mesobilirubinogen", 
    "Mesobilirubinogen\ncaecum", 10, "cecum_lcms", "Urobilinogen_Urobilin", "Urobilinogen/\nurobilin")
pathway_specs <- tibble::tribble(~pathway_order, ~Library, ~Pathway, ~pathway_label, 1, "HALLMARK", "HALLMARK_E2F_TARGETS", 
    "E2F targets", 2, "HALLMARK", "HALLMARK_G2M_CHECKPOINT", "G2M checkpoint", 3, "HALLMARK", "HALLMARK_MYC_TARGETS_V1", 
    "MYC targets", 4, "HALLMARK", "HALLMARK_MITOTIC_SPINDLE", "Mitotic spindle", 5, "HALLMARK", "HALLMARK_MTORC1_SIGNALING", 
    "mTORC1 signalling", 6, "HALLMARK", "HALLMARK_PROTEIN_SECRETION", "Protein secretion", 7, "HALLMARK", 
    "HALLMARK_CHOLESTEROL_HOMEOSTASIS", "Cholesterol homeostasis", 8, "REACTOME", "REACTOME_CHOLESTEROL_BIOSYNTHESIS", 
    "Cholesterol biosynthesis", 9, "HALLMARK", "HALLMARK_INTERFERON_ALPHA_RESPONSE", "Interferon-alpha response", 
    10, "HALLMARK", "HALLMARK_INTERFERON_GAMMA_RESPONSE", "Interferon-gamma response", 11, "HALLMARK", 
    "HALLMARK_TNFA_SIGNALING_VIA_NFKB", "TNF-NF-kB signalling", 12, "HALLMARK", "HALLMARK_INFLAMMATORY_RESPONSE", 
    "Inflammatory response", 13, "HALLMARK", "HALLMARK_IL6_JAK_STAT3_SIGNALING", "IL6-JAK-STAT3 signalling", 
    14, "HALLMARK", "HALLMARK_OXIDATIVE_PHOSPHORYLATION", "Oxidative phosphorylation")
build_selected_gsea_grid <- function(target_organ) {
    profile_specs <- dplyr::mutate(profile_specs_base, Organ = target_organ)
    selected_gsea_observed <- dplyr::select(dplyr::filter(dplyr::inner_join(dplyr::inner_join(gsea_results, 
        profile_specs, by = c("dataset_id", "Organ", "metabolite")), pathway_specs, by = c("Library", 
        "Pathway")), !is.na(.data$GSEA_FDR_family), .data$GSEA_FDR_family < family_fdr_cutoff), .data$Organ, 
        .data$profile_order, .data$profile_label, .data$pathway_order, .data$pathway_label, .data$NES, 
        .data$GSEA_FDR_family)
    duplicate_gsea_rows <- dplyr::filter(dplyr::count(selected_gsea_observed, .data$Organ, .data$profile_label, 
        .data$pathway_label), .data$n > 1L)
    if (nrow(duplicate_gsea_rows) > 0L) {
        stop("Duplicate selected GSEA rows were found for ", target_organ, ":\n", paste(utils::capture.output(print(duplicate_gsea_rows)), 
            collapse = "\n"))
    }
    selected_gsea_grid <- dplyr::mutate(dplyr::select(dplyr::mutate(dplyr::left_join(dplyr::left_join(dplyr::left_join(tidyr::expand_grid(Organ = target_organ, 
        profile_label = profile_specs$profile_label, pathway_label = pathway_specs$pathway_label), selected_gsea_observed, 
        by = c("Organ", "profile_label", "pathway_label")), dplyr::select(profile_specs, .data$Organ, 
        .data$profile_label, .data$profile_order), by = c("Organ", "profile_label"), suffix = c("", "_spec")), 
        dplyr::select(pathway_specs, .data$pathway_label, .data$pathway_order), by = "pathway_label", 
        suffix = c("", "_spec")), profile_order = dplyr::coalesce(.data$profile_order, .data$profile_order_spec), 
        pathway_order = dplyr::coalesce(.data$pathway_order, .data$pathway_order_spec)), -.data$profile_order_spec, 
        -.data$pathway_order_spec), Organ = factor(.data$Organ, levels = target_organ), profile_label = factor(.data$profile_label, 
        levels = profile_specs$profile_label), pathway_label = factor(.data$pathway_label, levels = rev(pathway_specs$pathway_label)), 
        significant = !is.na(.data$NES), NES_raw = .data$NES, NES_plot = .data$NES)
    selected_gsea_grid
}
selected_gsea_grid_colon <- build_selected_gsea_grid("Colon")
selected_gsea_grid_liver <- build_selected_gsea_grid("Liver")
panel_C_organ_label <- "Colon"
panel_E_organ_label <- "Liver"
make_gsea_panel <- function(gsea_grid, organ_label) {
    ggplot2::ggplot(gsea_grid, ggplot2::aes(x = .data$profile_label, y = .data$pathway_label, fill = .data$NES_plot)) + 
        ggplot2::geom_tile(colour = "white", linewidth = 0.45) + ggplot2::scale_fill_gradient2(low = "#2166AC", 
        mid = "white", high = "#B2182B", midpoint = 0, na.value = "#F2F2F2", name = "NES") + ggplot2::labs(x = NULL, 
        y = NULL) + manuscript_theme + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 40, 
        hjust = 1, vjust = 1, size = 7.4), axis.text.y = ggplot2::element_text(size = 7.5), legend.position = "bottom")
}
genus_levels <- c("Agathobaculum", "Alistipes", "Colidextribacter")
panel_B_rows <- dplyr::arrange(dplyr::summarise(dplyr::group_by(selected_gm, .data$theme, .data$theme_order, 
    .data$heatmap_row), plot_order = min(.data$plot_order), .groups = "drop"), .data$theme_order, .data$plot_order)
panel_B_values <- dplyr::distinct(dplyr::select(selected_gm, .data$theme, .data$heatmap_row, .data$genus, 
    .data$beta, .data$FDR))
duplicate_panel_B_values <- dplyr::filter(dplyr::count(panel_B_values, .data$theme, .data$heatmap_row, 
    .data$genus), .data$n > 1L)
if (nrow(duplicate_panel_B_values) > 0L) {
    stop("Duplicate genus-metabolite heatmap coordinates remain after ", "collapsing metabolite/platform rows:\n", 
        paste(utils::capture.output(print(duplicate_panel_B_values)), collapse = "\n"))
}
panel_B_row_levels <- dplyr::pull(panel_B_rows, .data$heatmap_row)
panel_B_data <- dplyr::mutate(dplyr::left_join(tidyr::expand_grid(panel_B_rows, genus = genus_levels), 
    panel_B_values, by = c("theme", "heatmap_row", "genus")), theme = factor(.data$theme, levels = theme_levels), 
    genus = factor(.data$genus, levels = genus_levels), heatmap_row = factor(.data$heatmap_row, levels = rev(panel_B_row_levels)))
if (anyDuplicated(dplyr::select(panel_B_data, .data$theme, .data$heatmap_row, .data$genus)) > 0L) {
    stop("Panel A still contains duplicated heatmap coordinates.")
}
panel_B <- ggplot2::ggplot(panel_B_data, ggplot2::aes(x = .data$genus, y = .data$heatmap_row, fill = .data$beta)) + 
    ggplot2::geom_tile(colour = "white", linewidth = 0.45) + ggplot2::facet_grid(rows = ggplot2::vars(.data$theme), 
    scales = "free_y", space = "free_y", switch = "y") + ggplot2::scale_fill_gradient2(low = "#2166AC", 
    mid = "white", high = "#B2182B", midpoint = 0, na.value = "#F2F2F2", name = "beta") + ggplot2::scale_x_discrete(drop = FALSE) + 
    ggplot2::labs(x = NULL, y = NULL) + manuscript_theme + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 35, 
    hjust = 1, size = 7.2), axis.text.y = ggplot2::element_text(size = 7), strip.placement = "outside", 
    strip.text.y.left = ggplot2::element_text(angle = 0, size = 7.2), legend.position = "bottom", plot.margin = ggplot2::margin(5, 
        5, 5, 18))
panel_C <- make_gsea_panel(gsea_grid = selected_gsea_grid_colon, organ_label = panel_C_organ_label)
panel_D_data <- dplyr::mutate(dplyr::arrange(dplyr::filter(selected_gm, .data$show_in_D), .data$theme_order, 
    .data$plot_order), term = paste0(.data$genus, " - ", .data$metabolite_short), term_key = paste0(.data$theme_order, 
    "__", .data$plot_order, "__", .data$term))
panel_D_levels <- rev(dplyr::pull(panel_D_data, .data$term_key))
panel_D_data <- dplyr::mutate(panel_D_data, term_key = factor(.data$term_key, levels = panel_D_levels))
panel_D <- ggplot2::ggplot(panel_D_data, ggplot2::aes(x = .data$beta, y = .data$term_key)) + ggplot2::geom_vline(xintercept = 0, 
    linetype = "dashed", linewidth = 0.4, colour = "#6F6F6F") + ggplot2::geom_segment(ggplot2::aes(x = .data$conf_low, 
    xend = .data$conf_high, yend = .data$term_key), linewidth = 0.55, colour = "#555555") + ggplot2::geom_point(ggplot2::aes(colour = .data$theme), 
    size = 2.4) + ggplot2::facet_grid(rows = ggplot2::vars(.data$theme), scales = "free_y", space = "free_y") + 
    ggplot2::scale_colour_manual(values = theme_palette, guide = "none") + ggplot2::scale_y_discrete(labels = function(x) {
    sub("^[0-9]+__[0-9]+__", "", x)
}) + ggplot2::labs(x = "Association coefficient (beta)", y = NULL) + manuscript_theme + ggplot2::theme(axis.text.x = ggplot2::element_text(size = 7.5), 
    axis.text.y = ggplot2::element_text(size = 6.9), axis.title.x = ggplot2::element_text(size = 8.5), 
    strip.text.y = ggplot2::element_text(angle = 0, size = 7))
hpp_names <- c(serum_lcms = "3-(3-Hydroxyphenyl)propanoic acid", cecum_lcms = "3-(3-Hydroxyphenyl)propionic acid", 
    cecum_nmr = "3-hydroxyphenylpropionate")
read_hpp_gene_statistics <- function(dataset_id, organ, path) {
    target_metabolite <- unname(hpp_names[[dataset_id]])
    dat <- dplyr::select(dplyr::filter(readr::read_csv(ppi_input(path), col_types = readr::cols_only(metabolite = readr::col_character(), 
        model = readr::col_character(), gene = readr::col_character(), statistic = readr::col_double()), 
        progress = FALSE), .data$model == model_label, .data$metabolite == target_metabolite, !is.na(.data$gene), 
        is.finite(.data$statistic)), .data$gene, .data$statistic)
    if (nrow(dat) == 0L) {
        stop("No 3-HPP rows were found for ", dataset_id, " / ", organ, " in:\n", path)
    }
    names(dat)[names(dat) == "statistic"] <- paste0("stat_", dataset_id)
    dat
}
gene_statistics <- dplyr::mutate(gene_files, hpp_data = Map(read_hpp_gene_statistics, .data$dataset_id, 
    .data$Organ, .data$path))
get_hpp_data <- function(dataset_id, organ) {
    idx <- which(gene_statistics$dataset_id == dataset_id & gene_statistics$Organ == organ)
    if (length(idx) != 1L) {
        stop("Could not uniquely identify ", dataset_id, " / ", organ)
    }
    gene_statistics$hpp_data[[idx]]
}
comparison_specs <- tibble::tribble(~Organ, ~comparison_dataset, ~comparison_label, "Liver", "cecum_lcms", 
    "Serum vs caecal LC-MS", "Liver", "cecum_nmr", "Serum LC-MS vs caecal NMR", "Colon", "cecum_lcms", 
    "Serum vs caecal LC-MS", "Colon", "cecum_nmr", "Serum LC-MS vs caecal NMR")
build_concordance_data <- function(organ, comparison_dataset, comparison_label) {
    serum_data <- get_hpp_data("serum_lcms", organ)
    comparison_data <- get_hpp_data(comparison_dataset, organ)
    merged <- dplyr::inner_join(serum_data, comparison_data, by = "gene")
    comparison_column <- paste0("stat_", comparison_dataset)
    rho <- stats::cor(merged$stat_serum_lcms, merged[[comparison_column]], method = "spearman", use = "complete.obs")
    tibble::tibble(gene = merged$gene, serum_statistic = merged$stat_serum_lcms, comparison_statistic = merged[[comparison_column]], 
        Organ = organ, comparison_dataset = comparison_dataset, comparison_label = comparison_label, 
        spearman_rho = rho, n_genes = nrow(merged))
}
panel_F_data <- dplyr::bind_rows(Map(build_concordance_data, comparison_specs$Organ, comparison_specs$comparison_dataset, 
    comparison_specs$comparison_label))
panel_F_labels <- dplyr::mutate(dplyr::distinct(panel_F_data, .data$Organ, .data$comparison_label, .data$spearman_rho, 
    .data$n_genes), label = paste0("rho = ", sprintf("%.3f", .data$spearman_rho)))
panel_F <- ggplot2::ggplot(panel_F_data, ggplot2::aes(x = .data$serum_statistic, y = .data$comparison_statistic)) + 
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", linewidth = 0.45, colour = "#6F6F6F") + 
    ggplot2::geom_point(size = 0.55, alpha = 0.18, colour = "#333333") + ggplot2::geom_text(data = panel_F_labels, 
    ggplot2::aes(x = -Inf, y = Inf, label = .data$label), inherit.aes = FALSE, hjust = -0.08, vjust = 1.3, 
    size = 3) + ggplot2::facet_grid(rows = ggplot2::vars(.data$Organ), cols = ggplot2::vars(.data$comparison_label), 
    scales = "free") + ggplot2::labs(x = "Serum 3-HPP Wald statistic", y = "Comparison-layer 3-HPP Wald statistic") + 
    manuscript_theme + ggplot2::theme(axis.text = ggplot2::element_text(size = 7.3), axis.title = ggplot2::element_text(size = 8.5), 
    strip.text = ggplot2::element_text(size = 8), panel.spacing = grid::unit(1.5, "mm"))
panel_E <- make_gsea_panel(gsea_grid = selected_gsea_grid_liver, organ_label = panel_E_organ_label)
row_1 <- cowplot::plot_grid(panel_D, panel_F, nrow = 1, rel_widths = c(1.1, 1.45), labels = c("A", "B"), 
    label_size = 14, label_fontface = "bold", label_x = 0, label_y = 1, hjust = 0, vjust = 1, align = "none")
row_2 <- cowplot::plot_grid(panel_C, nrow = 1, labels = "C", label_size = 14, label_fontface = "bold", 
    label_x = 0, label_y = 1, hjust = 0, vjust = 1, align = "none")
row_3 <- cowplot::plot_grid(panel_E, nrow = 1, labels = "D", label_size = 14, label_fontface = "bold", 
    label_x = 0, label_y = 1, hjust = 0, vjust = 1, align = "none")
final_figure <- cowplot::plot_grid(row_1, row_2, row_3, ncol = 1, rel_heights = c(1.1, 1.15, 1.15), align = "none")
png_file <- file.path(output_dir, "Figure_crosslayer_focus_concordance_colon_liver_v3_no_titles.png")
pdf_file <- file.path(output_dir, "Figure_crosslayer_focus_concordance_colon_liver_v3_no_titles.pdf")
if (requireNamespace("ragg", quietly = TRUE)) {
    ragg::agg_png(filename = ppi_output(png_file), width = figure_width_mm, height = figure_height_mm, 
        units = "mm", res = figure_dpi, background = "white")
    print(final_figure)
    grDevices::dev.off()
} else {
    grDevices::png(filename = ppi_output(png_file), width = figure_width_mm, height = figure_height_mm, 
        units = "mm", res = figure_dpi, bg = "white")
    print(final_figure)
    grDevices::dev.off()
}
grDevices::pdf(file = ppi_output(pdf_file), width = figure_width_mm/25.4, height = figure_height_mm/25.4, 
    family = "Helvetica", useDingbats = FALSE, paper = "special", onefile = TRUE)
print(final_figure)
grDevices::dev.off()
readr::write_csv(panel_D_data, ppi_output(file.path(output_dir, "panel_A_focused_effect_data.csv"), name = "panel_a_kieled_effect_data.csv"))
readr::write_csv(panel_F_labels, ppi_output(file.path(output_dir, "panel_B_3HPP_correlations.csv"), name = "panel_b_3_hpp_correlations.csv"))
readr::write_csv(selected_gsea_grid_colon, ppi_output(file.path(output_dir, "panel_C_colon_DESeq2_GSEA_data.csv"), 
    name = "panel_c_colon_de_seq2_gsea_data.csv"))
readr::write_csv(selected_gsea_grid_liver, ppi_output(file.path(output_dir, "panel_D_liver_DESeq2_GSEA_data.csv"), 
    name = "panel_d_liver_de_seq2_gsea_data.csv"))
source_audit <- tibble::tribble(~final_panel, ~analysis, ~organ, ~source_pipeline, ~source, "A", "Focused genus-metabolite effects", 
    NA_character_, "02_association_results", basename(association_file), "B", "3-HPP gene-statistic concordance", 
    "Colon and Liver", "05_gene_metabolite_DESeq2_Group_Metabolite_Hallmark_Reactome", "complete_gene_results/gene_results__*.csv.gz", 
    "C", "Metabolite-pathway heatmap", panel_C_organ_label, "05_gene_metabolite_DESeq2_Group_Metabolite_Hallmark_Reactome", 
    basename(gsea_file), "D", "Metabolite-pathway heatmap", panel_E_organ_label, "05_gene_metabolite_DESeq2_Group_Metabolite_Hallmark_Reactome", 
    basename(gsea_file))
readr::write_csv(source_audit, ppi_output(file.path(output_dir, "Figure_crosslayer_focus_concordance_colon_liver_v3_no_titles_source_audit.csv"), 
    name = "figure_crosslayer_kiel_concordance_colon_liver_no_titles_source_audit.csv"))
message("Figure created successfully:\n", png_file, "\n", pdf_file)
if (interactive()) {
    print(final_figure)
}
