# SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0
# Export of exact plotted matrices
# Run through run.R so input and output paths resolve within the package.
# Inputs: see data/input_manifest.csv and docs/RUN_ORDER.md.
# Writes generated files under results/; bundled inputs are read-only.

ppi_setwd(ppi_path("metabolomics"))
required_packages <- c("dplyr", "tidyr", "stringr", "tibble", "readr", "readxl", "ggplot2", "cowplot", 
    "vegan", "limma", "statmod", "forcats")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0L) {
    stop("Please install these packages first: ", paste(missing_packages, collapse = ", "))
}
config_file <- ppi_code("scripts/06_multiomics_integration/06_01_configure_multiomics.R")
ppi_source(config_file, local = FALSE)
required_config_objects <- c("project_dir", "metabolomics_root", "association_root", "input_paths")
missing_config_objects <- setdiff(required_config_objects, ls(envir = .GlobalEnv))
if (length(missing_config_objects) > 0L) {
    stop("The configuration file did not create: ", paste(missing_config_objects, collapse = ", "))
}
output_dir <- file.path(metabolomics_root, "Figure3_remade_mouse_dupCor_combined_significance")
dir.create(ppi_output(output_dir, directory = TRUE), recursive = TRUE, showWarnings = FALSE)
file_nmr_exvivo_raw <- input_paths$exvivo_nmr_raw
file_nmr_exvivo_map <- input_paths$exvivo_nmr_mapping
file_hilic_exvivo_rds <- input_paths$exvivo_hilic
file_cecum_nmr_workbook <- input_paths$cecum_nmr
file_cf_lcms_rds <- input_paths$cecum_feces_lcms
file_serum_rds <- input_paths$serum_lcms
mouse_lcms_heatmap_results <- "mouse"
resolve_exact_result <- function(candidates, label) {
    candidates <- unique(candidates)
    hits <- candidates[ppi_file_exists(candidates)]
    if (length(hits) != 1L) {
        stop("Could not uniquely identify ", label, ".\n", "Checked:\n", paste(candidates, collapse = "\n"))
    }
    hits[[1]]
}
file_nmr_exvivo_results <- resolve_exact_result(c(file.path(dirname(file_nmr_exvivo_map), "NMR_exvivo_annotated_group_results_v2.csv")), 
    "ex vivo NMR group results")
file_hilic_exvivo_results <- resolve_exact_result(c(file.path(dirname(file_hilic_exvivo_rds), "HILIC_exvivo_primary_rank1_rank2_primary_duplicateCorrelation.csv"), 
    file.path(dirname(file_hilic_exvivo_rds), "Limma_analysis", "HILIC_exvivo_primary_rank1_rank2_primary_duplicateCorrelation.csv")), 
    "ex vivo HILIC results")
serum_result_filename <- "Serum_limma_duplicateCorrelation_results.csv"
file_serum_results <- resolve_exact_result(c(file.path(dirname(file_serum_rds), "Limma_analysis", serum_result_filename), 
    file.path(dirname(file_serum_rds), serum_result_filename)), paste("serum", mouse_lcms_heatmap_results, 
    "results"))
cf_result_suffix <- "duplicateCorrelation_results.csv"
file_cf_cecum_results <- resolve_exact_result(c(file.path(dirname(file_cf_lcms_rds), "Limma_analysis", 
    paste0("cecum_one_per_metabolite_", cf_result_suffix))), paste("caecum LC-MS", mouse_lcms_heatmap_results, 
    "results"))
file_cf_feces_results <- resolve_exact_result(c(file.path(dirname(file_cf_lcms_rds), "Limma_analysis", 
    paste0("feces_one_per_metabolite_", cf_result_suffix))), paste("faeces LC-MS", mouse_lcms_heatmap_results, 
    "results"))
file_cf_combined_results <- resolve_exact_result(c(file.path(dirname(file_cf_lcms_rds), "Limma_analysis", 
    "Cecum_Feces_LCMS_combined_compartment_limma_FINAL_results.xlsx"), file.path(dirname(file_cf_lcms_rds), 
    "Cecum_Feces_LCMS_combined_compartment_limma_FINAL_results.xlsx")), "combined caecum/faeces LC-MS compartment-adjusted results")
heatmap_mode_exvivo <- "nominal_only"
heatmap_mode_intestinal <- "significant_only"
heatmap_mode_serum <- "significant_only"
intestinal_top_n_per_dataset <- 12L
serum_top_n <- 30L
cecum_nmr_heatmap_results <- "mouse"
figure_width_mm <- 360
figure_height_mm <- 340
figure_dpi <- 600
check_file <- function(path) {
    if (!ppi_file_exists(path)) 
        stop("File not found: ", path)
}
safe_read_csv <- function(path) {
    check_file(path)
    readr::read_csv(ppi_input(path), show_col_types = FALSE, progress = FALSE)
}
safe_read_xlsx <- function(path, sheet = 1) {
    check_file(path)
    readxl::read_xlsx(ppi_input(path), sheet = sheet, .name_repair = "minimal")
}
make_feature_pseudocount <- function(x) {
    positive <- x[is.finite(x) & x > 0]
    if (length(positive) == 0L) 
        return(1e-12)
    min(positive, na.rm = TRUE)/2
}
log2_transform_matrix <- function(X, negative_to_zero = FALSE) {
    X <- as.matrix(X)
    storage.mode(X) <- "numeric"
    if (negative_to_zero) 
        X[X < 0] <- 0
    pseudocounts <- apply(X, 1, make_feature_pseudocount)
    X_log <- X
    for (i in seq_len(nrow(X_log))) {
        X_log[i, ] <- log2(X_log[i, ] + pseudocounts[i])
    }
    X_log
}
zscore_rows <- function(X) {
    X <- as.matrix(X)
    t(scale(t(X)))
}
median_impute_cols <- function(df) {
    df <- as.data.frame(df, check.names = FALSE)
    for (j in seq_len(ncol(df))) {
        x <- df[[j]]
        if (is.numeric(x)) {
            med <- stats::median(x, na.rm = TRUE)
            if (!is.finite(med)) 
                med <- 0
            x[is.na(x)] <- med
            df[[j]] <- x
        }
    }
    df
}
scale_features_for_pca <- function(sample_feature_df) {
    X <- as.data.frame(sample_feature_df, check.names = FALSE)
    X <- median_impute_cols(X)
    keep <- vapply(X, function(x) is.numeric(x) && stats::sd(x, na.rm = TRUE) > 0, logical(1))
    X <- X[, keep, drop = FALSE]
    if (ncol(X) < 2L) 
        stop("Fewer than two variable features remain for PCA.")
    scale(X)
}
make_sig_stars <- function(fdr) {
    dplyr::case_when(is.na(fdr) ~ "", fdr < 0.01 ~ "***", fdr < 0.05 ~ "**", fdr < 0.1 ~ "*", TRUE ~ 
        "")
}
clean_metabolite_name_cecum <- function(x) {
    x <- stringr::str_squish(as.character(x))
    x <- dplyr::na_if(x, "")
    xl <- stringr::str_to_lower(x)
    dplyr::case_when(is.na(x) ~ NA_character_, xl == "3-(3-oh-phenyl)propionate" ~ "3-hydroxyphenylpropionate", 
        xl == "5-aminopentanoic acid" ~ "5-aminopentanoate", xl == "phosphocoline" ~ "phosphocholine", 
        TRUE ~ xl)
}
capitalize_display_label <- function(x) {
    vapply(x, function(label) {
        if (is.na(label) || !nzchar(label)) 
            return(NA_character_)
        parts <- stringr::str_split(label, "\\s*/\\s*", simplify = FALSE)[[1]]
        parts <- vapply(parts, function(part) {
            location <- stringr::str_locate(part, "[[:alpha:]]")[1, ]
            if (any(is.na(location))) 
                return(part)
            paste0(stringr::str_sub(part, 1, location[[1]] - 1L), stringr::str_to_upper(stringr::str_sub(part, 
                location[[1]], location[[2]])), stringr::str_sub(part, location[[2]] + 1L))
        }, FUN.VALUE = character(1))
        paste(parts, collapse = " / ")
    }, FUN.VALUE = character(1))
}
standardize_metabolite_label <- function(x) {
    x <- stringr::str_squish(as.character(x))
    x <- dplyr::na_if(x, "")
    key <- stringr::str_to_lower(x)
    standardized <- dplyr::case_when(is.na(x) ~ NA_character_, key %in% c("3-hydroxyphenylpropionate", 
        "3-(3-hydroxyphenyl)propionate", "3-(3-hydroxyphenyl)propionic acid", "3-(3-hydroxyphenyl)propanoic acid", 
        "3-(3-oh-phenyl)propionate", "3-hpp") ~ "3-HPP", key %in% c("3-(4-hydroxyphenyl)propionate", 
        "3-(4-hydroxyphenyl)propionic acid", "3-(4-hydroxyphenyl)propanoic acid", "4-hydroxyphenylpropionate", 
        "4-hydroxyphenylpropionic acid", "hppa", "4-hppa") ~ "4-HPPA", key %in% c("isocaproate", "isocaproic acid", 
        "4-methylvalerate", "4-methylvaleric acid", "4-methylpentanoate", "4-methylpentanoic acid") ~ 
        "Isocaproate", stringr::str_detect(key, "^3-\\(4-hydroxyphenyl\\)[ -]?propionate\\s*/\\s*3-phenyllactate") ~ 
        "4-HPPA / 3-Phenyllactate (unresolved)", TRUE ~ x)
    standardized <- stringr::str_replace_all(standardized, stringr::regex("4-methylvaleric acid|4-methylvalerate|4-methylpentanoic acid|4-methylpentanoate|isocaproic acid", 
        ignore_case = TRUE), "Isocaproate")
    capitalize_display_label(standardized)
}
base_theme <- ggplot2::theme_bw(base_size = 10) + ggplot2::theme(panel.grid = ggplot2::element_line(colour = "#E6E6E6", 
    linewidth = 0.35), panel.border = ggplot2::element_rect(colour = "#4A4A4A", fill = NA, linewidth = 0.55), 
    axis.ticks = ggplot2::element_line(colour = "#4A4A4A", linewidth = 0.35), axis.line = ggplot2::element_line(colour = "#4A4A4A", 
        linewidth = 0.35), strip.background = ggplot2::element_rect(fill = "#DDDDDD", colour = "#4A4A4A", 
        linewidth = 0.45), strip.text = ggplot2::element_text(face = "plain", size = 8), legend.background = ggplot2::element_blank(), 
    legend.key = ggplot2::element_blank(), plot.margin = ggplot2::margin(4, 4, 4, 4))
make_pca_panel <- function(sample_feature_df, metadata, colour_col = "Group", shape_col = NULL, panel_label_text = NULL, 
    permanova_rhs = "Group", annotation_position = c("top_right", "bottom_right", "top_left", "bottom_left")) {
    annotation_position <- match.arg(annotation_position)
    X_scaled <- scale_features_for_pca(sample_feature_df)
    pca <- stats::prcomp(X_scaled, center = FALSE, scale. = FALSE)
    var_exp <- 100 * (pca$sdev^2/sum(pca$sdev^2))
    scores <- as.data.frame(pca$x[, 1:2, drop = FALSE])
    scores <- dplyr::bind_cols(metadata, scores)
    dist_mat <- stats::dist(X_scaled)
    permanova_formula <- stats::as.formula(paste("dist_mat ~", permanova_rhs))
    ad <- vegan::adonis2(permanova_formula, data = metadata, permutations = 999)
    r2 <- unname(ad$R2[1])
    pval <- unname(ad$`Pr(>F)`[1])
    p <- ggplot2::ggplot(scores, ggplot2::aes(x = .data$PC1, y = .data$PC2, colour = .data[[colour_col]]))
    if (!is.null(shape_col) && shape_col %in% names(scores)) {
        p <- p + ggplot2::geom_point(ggplot2::aes(shape = .data[[shape_col]]), size = 2.5)
    }
    else {
        p <- p + ggplot2::geom_point(size = 2.5)
    }
    p <- p + ggplot2::scale_colour_manual(values = c(Ctrl = "#E67E22", PPI = "#1F4E8C", ctrl = "#E67E22", 
        ppi = "#1F4E8C"), name = "Group")
    if (!is.null(shape_col) && shape_col %in% names(scores)) {
        shape_values <- NULL
        vals <- unique(as.character(stats::na.omit(scores[[shape_col]])))
        if (all(vals %in% c("Day 1", "Day 2", "Day 3", "Day 4"))) {
            shape_values <- c(`Day 1` = 16, `Day 2` = 15, `Day 3` = 17, `Day 4` = 3)
        }
        else if (all(vals %in% c("Caecum", "Faeces"))) {
            shape_values <- c(Caecum = 16, Faeces = 17)
        }
        if (!is.null(shape_values)) {
            p <- p + ggplot2::scale_shape_manual(values = shape_values, name = shape_col)
        }
    }
    ann <- paste0(if (!is.null(panel_label_text)) 
        paste0(panel_label_text, "\n")
    else "", "PERMANOVA\nR² = ", sprintf("%.3f", r2), "\np = ", formatC(pval, digits = 3, format = "f"))
    annotation_parameters <- switch(annotation_position, top_right = list(x = Inf, y = Inf, hjust = 1.05, 
        vjust = 1.1), bottom_right = list(x = Inf, y = -Inf, hjust = 1.05, vjust = -0.1), top_left = list(x = -Inf, 
        y = Inf, hjust = -0.05, vjust = 1.1), bottom_left = list(x = -Inf, y = -Inf, hjust = -0.05, vjust = -0.1))
    p + ggplot2::annotate("text", x = annotation_parameters$x, y = annotation_parameters$y, label = ann, 
        hjust = annotation_parameters$hjust, vjust = annotation_parameters$vjust, size = 3) + ggplot2::labs(x = paste0("PC1 (", 
        sprintf("%.1f", var_exp[1]), "%)"), y = paste0("PC2 (", sprintf("%.1f", var_exp[2]), "%)")) + 
        base_theme + ggplot2::theme(axis.text = ggplot2::element_text(size = 7.2), axis.title = ggplot2::element_text(size = 8.2), 
        legend.title = ggplot2::element_text(size = 8), legend.text = ggplot2::element_text(size = 7.4))
}
build_exvivo_nmr_matrix <- function(raw_csv, mapping_csv) {
    raw_df <- safe_read_csv(raw_csv)
    map_df <- safe_read_csv(mapping_csv)
    sample_cols <- names(raw_df)[stringr::str_detect(names(raw_df), "^[1-4]_[CP][1-3]$")]
    if (length(sample_cols) == 0L) {
        stop("Could not identify ex vivo NMR day 1-4 sample columns.")
    }
    feature_id_col <- names(raw_df)[stringr::str_detect(tolower(names(raw_df)), "classy")][1]
    if (is.na(feature_id_col)) 
        stop("Could not identify CLASSY feature-ID column in ex vivo raw NMR file.")
    raw_use <- dplyr::transmute(raw_df, feature_id = as.character(.data[[feature_id_col]]), dplyr::across(dplyr::all_of(sample_cols), 
        as.numeric))
    map_use <- dplyr::filter(dplyr::mutate(map_df, feature_id = as.character(.data$feature_id)), .data$include_group_level)
    merged <- dplyr::inner_join(raw_use, map_use, by = "feature_id")
    if (nrow(merged) == 0L) 
        stop("No ex vivo NMR features matched the cleaned feature mapping.")
    X <- as.matrix(merged[, sample_cols, drop = FALSE])
    rownames(X) <- merged$feature_id
    X_log <- log2_transform_matrix(X, negative_to_zero = TRUE)
    X_z <- zscore_rows(X_log)
    groups <- split(merged$feature_id, merged$analysis_group)
    score_mat <- vapply(groups, function(ids) {
        ids <- intersect(ids, rownames(X_z))
        if (length(ids) == 1L) {
            as.numeric(X_z[ids, ])
        }
        else {
            apply(X_z[ids, , drop = FALSE], 2, stats::median, na.rm = TRUE)
        }
    }, numeric(length(sample_cols)))
    score_mat <- t(score_mat)
    colnames(score_mat) <- sample_cols
    score_mat <- t(score_mat)
    meta <- dplyr::mutate(tibble::tibble(sample_col = sample_cols), day = stringr::str_extract(.data$sample_col, 
        "^[1-4]"), donor = stringr::str_extract(.data$sample_col, "[CP][1-3]$"), Group = dplyr::if_else(stringr::str_starts(.data$donor, 
        "C"), "Ctrl", "PPI"), Day = factor(paste("Day", .data$day), levels = c("Day 1", "Day 2", "Day 3", 
        "Day 4")))
    list(matrix = as.data.frame(score_mat, check.names = FALSE), metadata = meta)
}
load_cecum_nmr_from_results_workbook <- function(workbook_path) {
    check_file(workbook_path)
    mouse_values <- readxl::read_xlsx(ppi_input(workbook_path), sheet = "Mouse_level_values", .name_repair = "minimal")
    sample_metadata <- readxl::read_xlsx(ppi_input(workbook_path), sheet = "Sample_metadata", .name_repair = "minimal")
    required_value_columns <- c("sample_id", "analysis_group", "score")
    missing_value_columns <- setdiff(required_value_columns, names(mouse_values))
    if (length(missing_value_columns) > 0L) {
        stop("Caecum NMR Mouse_level_values is missing: ", paste(missing_value_columns, collapse = ", "))
    }
    duplicate_values <- dplyr::filter(dplyr::count(mouse_values, .data$sample_id, .data$analysis_group), 
        .data$n > 1L)
    if (nrow(duplicate_values) > 0L) {
        stop("Duplicate sample/metabolite rows were found in the caecum NMR Mouse_level_values sheet.")
    }
    wide_values <- tidyr::pivot_wider(dplyr::transmute(mouse_values, sample_id = as.character(.data$sample_id), 
        analysis_group = as.character(.data$analysis_group), score = as.numeric(.data$score)), names_from = .data$analysis_group, 
        values_from = .data$score)
    required_metadata_columns <- c("sample_id", "inoculum_donor", "group")
    missing_metadata_columns <- setdiff(required_metadata_columns, names(sample_metadata))
    if (length(missing_metadata_columns) > 0L) {
        stop("Caecum NMR Sample_metadata is missing: ", paste(missing_metadata_columns, collapse = ", "))
    }
    metadata <- dplyr::transmute(sample_metadata, sample_id = as.character(.data$sample_id), donor = as.character(.data$inoculum_donor), 
        Group = dplyr::if_else(stringr::str_to_lower(as.character(.data$group)) == "ctrl", "Ctrl", "PPI"))
    metadata <- metadata[match(wide_values$sample_id, metadata$sample_id), , drop = FALSE]
    if (anyNA(metadata$sample_id)) {
        stop("Some caecum NMR samples could not be matched to Sample_metadata.")
    }
    X <- as.data.frame(dplyr::select(wide_values, -.data$sample_id), check.names = FALSE)
    rownames(X) <- wide_values$sample_id
    list(matrix = X, metadata = metadata)
}
load_serum_matrix <- function(rds_path) {
    check_file(rds_path)
    obj <- readRDS(ppi_input(rds_path))
    if (!all(c("serum_clean_bio", "metadata_clean_bio") %in% names(obj))) {
        stop("Serum RDS does not contain serum_clean_bio and metadata_clean_bio.")
    }
    X <- as.data.frame(obj$serum_clean_bio, check.names = FALSE)
    meta <- as.data.frame(obj$metadata_clean_bio, check.names = FALSE)
    sample_col <- if ("Sample_Name" %in% names(meta)) 
        "Sample_Name"
    else names(meta)[1]
    if (!identical(rownames(X), as.character(meta[[sample_col]]))) {
        m <- match(rownames(X), as.character(meta[[sample_col]]))
        meta <- meta[m, , drop = FALSE]
    }
    group_col <- if ("Separation" %in% names(meta)) 
        "Separation"
    else if ("group" %in% names(meta)) 
        "group"
    else stop("Could not find serum group column.")
    meta$Group <- dplyr::if_else(stringr::str_to_lower(as.character(meta[[group_col]])) == "ctrl", "Ctrl", 
        "PPI")
    X_matrix <- as.matrix(X)
    storage.mode(X_matrix) <- "numeric"
    pseudocounts <- apply(X_matrix, 2, make_feature_pseudocount)
    X_log <- log2(sweep(X_matrix, 2, pseudocounts, "+"))
    X_samples <- as.data.frame(X_log, check.names = FALSE)
    list(matrix = X_samples, metadata = meta)
}
load_cfeces_lcms_for_pca <- function(rds_path, representation = "one_per_metabolite") {
    check_file(rds_path)
    obj <- readRDS(ppi_input(rds_path))
    needed <- c("rt_resolved", "one_per_metabolite", "metadata")
    if (!all(needed %in% names(obj))) {
        stop("Cecum/Faeces LC-MS RDS does not contain expected top-level objects: ", paste(needed, collapse = ", "))
    }
    Xc <- as.data.frame(obj[[representation]]$cecum, check.names = FALSE)
    Xf <- as.data.frame(obj[[representation]]$feces, check.names = FALSE)
    mc <- as.data.frame(obj$metadata$cecum, check.names = FALSE)
    mf <- as.data.frame(obj$metadata$feces, check.names = FALSE)
    common_cols <- intersect(colnames(Xc), colnames(Xf))
    Xc <- Xc[, common_cols, drop = FALSE]
    Xf <- Xf[, common_cols, drop = FALSE]
    group_col_c <- if ("Separation" %in% names(mc)) 
        "Separation"
    else "group"
    group_col_f <- if ("Separation" %in% names(mf)) 
        "Separation"
    else "group"
    mc$Group <- dplyr::if_else(stringr::str_to_lower(as.character(mc[[group_col_c]])) == "ctrl", "Ctrl", 
        "PPI")
    mf$Group <- dplyr::if_else(stringr::str_to_lower(as.character(mf[[group_col_f]])) == "ctrl", "Ctrl", 
        "PPI")
    mc$Organ <- "Caecum"
    mf$Organ <- "Faeces"
    Xall <- dplyr::bind_rows(Xc, Xf)
    mall <- dplyr::bind_rows(mc, mf)
    Xall_m <- as.matrix(Xall)
    storage.mode(Xall_m) <- "numeric"
    pseudo <- apply(Xall_m, 2, function(x) {
        pos <- x[is.finite(x) & x > 0]
        if (length(pos) == 0) 
            1e-12
        else min(pos)/2
    })
    Xall_log <- log2(sweep(Xall_m, 2, pseudo, "+"))
    list(matrix = as.data.frame(Xall_log, check.names = FALSE), metadata = mall)
}
load_hilic_exvivo_for_pca <- function(rds_path) {
    check_file(rds_path)
    obj <- readRDS(ppi_input(rds_path))
    required_objects <- c("primary_matrix_all_samples", "metadata_all")
    missing_objects <- setdiff(required_objects, names(obj))
    if (length(missing_objects) > 0L) {
        stop("HILIC ex vivo RDS is missing: ", paste(missing_objects, collapse = ", "))
    }
    X <- as.matrix(obj$primary_matrix_all_samples)
    metadata <- as.data.frame(obj$metadata_all, check.names = FALSE)
    storage.mode(X) <- "numeric"
    required_metadata <- c("Sample_Name", "is_analysis_sample", "is_medium", "is_day0", "Day", "Donor", 
        "Separation")
    missing_metadata <- setdiff(required_metadata, names(metadata))
    if (length(missing_metadata) > 0L) {
        stop("HILIC metadata is missing: ", paste(missing_metadata, collapse = ", "))
    }
    if (!identical(rownames(X), as.character(metadata$Sample_Name))) {
        sample_match <- match(rownames(X), as.character(metadata$Sample_Name))
        if (anyNA(sample_match)) {
            stop("HILIC matrix rows cannot be matched to metadata$Sample_Name.")
        }
        metadata <- metadata[sample_match, , drop = FALSE]
    }
    keep <- dplyr::coalesce(as.logical(metadata$is_analysis_sample), FALSE) & !dplyr::coalesce(as.logical(metadata$is_medium), 
        FALSE) & !dplyr::coalesce(as.logical(metadata$is_day0), FALSE)
    X <- X[keep, , drop = FALSE]
    metadata <- metadata[keep, , drop = FALSE]
    day <- as.integer(stringr::str_extract(as.character(metadata$Day), "[0-9]+"))
    metadata_out <- tibble::tibble(sample_id = rownames(X), donor = as.character(metadata$Donor), Group = dplyr::if_else(stringr::str_to_lower(as.character(metadata$Separation)) == 
        "ctrl", "Ctrl", "PPI"), Day = factor(paste("Day", day), levels = c("Day 1", "Day 2", "Day 3", 
        "Day 4")))
    pseudocounts <- apply(X, 2, make_feature_pseudocount)
    X_log <- log2(sweep(X, 2, pseudocounts, "+"))
    list(matrix = as.data.frame(X_log, check.names = FALSE), metadata = metadata_out)
}
run_lcms_duplicate_correlation <- function(X, metadata) {
    X <- as.matrix(X)
    storage.mode(X) <- "numeric"
    sample_col <- if ("Sample_Name" %in% names(metadata)) 
        "Sample_Name"
    else names(metadata)[1]
    group_col <- if ("Separation" %in% names(metadata)) 
        "Separation"
    else if ("group" %in% names(metadata)) 
        "group"
    else stop("No group column found.")
    donor_col <- if ("Donor" %in% names(metadata)) 
        "Donor"
    else stop("No donor column found.")
    if (!identical(rownames(X), as.character(metadata[[sample_col]]))) {
        mm <- match(rownames(X), as.character(metadata[[sample_col]]))
        metadata <- metadata[mm, , drop = FALSE]
    }
    group <- factor(stringr::str_to_lower(as.character(metadata[[group_col]])), levels = c("ctrl", "ppi"))
    donor <- factor(as.character(metadata[[donor_col]]))
    keep <- apply(X, 2, function(x) {
        x_finite <- x[is.finite(x)]
        variable <- length(x_finite) >= 2 && stats::sd(x_finite, na.rm = TRUE) > 0
        n_pos <- tapply(is.finite(x) & x > 0, group, sum)
        enough <- all(n_pos >= 2)
        variable && enough
    })
    X <- X[, keep, drop = FALSE]
    pseudo <- apply(X, 2, function(x) {
        pos <- x[is.finite(x) & x > 0]
        if (length(pos) == 0) 
            1e-12
        else min(pos)/2
    })
    X_log <- log2(sweep(X, 2, pseudo, "+"))
    design <- stats::model.matrix(~group)
    dupcor <- limma::duplicateCorrelation(t(X_log), design = design, block = donor)
    fit <- limma::lmFit(t(X_log), design = design, block = donor, correlation = dupcor$consensus)
    fit <- limma::eBayes(fit, robust = TRUE)
    coef_name <- "groupppi"
    coef_idx <- match(coef_name, colnames(fit$coefficients))
    est <- fit$coefficients[, coef_idx]
    se <- fit$stdev.unscaled[, coef_idx] * sqrt(fit$s2.post)
    df_total <- fit$df.total
    crit <- stats::qt(0.975, df = df_total)
    dplyr::arrange(tibble::tibble(analyte = rownames(fit$coefficients), log2_FC = as.numeric(est), fold_change = 2^log2_FC, 
        moderated_SE = as.numeric(se), conf_low = as.numeric(est - crit * se), conf_high = as.numeric(est + 
            crit * se), p_value = as.numeric(fit$p.value[, coef_idx]), FDR = stats::p.adjust(fit$p.value[, 
            coef_idx], method = "BH")), .data$p_value)
}
prepare_heatmap_long <- function(df, metabolite_col, effect_col, fdr_col, dataset_label, mode = c("top_ranked", 
    "significant_only", "nominal_only"), top_n = 10L, p_col = NULL, significance_source = "dataset-specific mouse duplicateCorrelation") {
    mode <- match.arg(mode)
    p_values <- if (is.null(p_col)) {
        rep(NA_real_, nrow(df))
    }
    else {
        if (!p_col %in% names(df)) {
            stop("Nominal P-value column '", p_col, "' was not found for dataset ", dataset_label, ".")
        }
        as.numeric(df[[p_col]])
    }
    out <- dplyr::filter(dplyr::transmute(df, Metabolite_raw = as.character(.data[[metabolite_col]]), 
        Metabolite = standardize_metabolite_label(.data[[metabolite_col]]), effect = as.numeric(.data[[effect_col]]), 
        P_value = p_values, FDR = as.numeric(.data[[fdr_col]]), Dataset = dataset_label, significance_source = significance_source), 
        !is.na(.data$Metabolite), !is.na(.data$effect))
    duplicate_labels <- dplyr::filter(dplyr::count(out, .data$Dataset, .data$Metabolite), .data$n > 1L)
    if (nrow(duplicate_labels) > 0L) {
        warning("Synonym harmonisation collapsed duplicate metabolite labels in ", dataset_label, ": ", 
            paste(duplicate_labels$Metabolite, collapse = ", "))
        out <- dplyr::select(dplyr::ungroup(dplyr::slice_head(dplyr::group_by(dplyr::arrange(dplyr::mutate(out, 
            .rank_p = dplyr::coalesce(.data$P_value, Inf), .rank_fdr = dplyr::coalesce(.data$FDR, Inf)), 
            .data$Dataset, .data$Metabolite, .data$.rank_p, .data$.rank_fdr, dplyr::desc(abs(.data$effect))), 
            .data$Dataset, .data$Metabolite), n = 1L)), -.data$.rank_p, -.data$.rank_fdr)
    }
    if (mode == "significant_only") {
        out <- dplyr::filter(out, !is.na(.data$FDR), .data$FDR < 0.1)
    }
    else if (mode == "nominal_only") {
        if (is.null(p_col)) {
            stop("mode = 'nominal_only' requires p_col for dataset ", dataset_label, ".")
        }
        out <- dplyr::arrange(dplyr::filter(out, !is.na(.data$P_value), .data$P_value < 0.05), dplyr::desc(.data$effect), 
            .data$P_value, .data$Metabolite)
    }
    else {
        out <- dplyr::slice_head(dplyr::arrange(out, dplyr::desc(abs(.data$effect)), .data$FDR, .data$Metabolite), 
            n = top_n)
    }
    out
}
build_heatmap_panel <- function(long_df, dataset_order, top_n_total = NULL, y_title = "Metabolite", tile_width = 0.78, 
    metabolite_font_size = 8) {
    long_df <- dplyr::ungroup(dplyr::mutate(dplyr::group_by(long_df, .data$Dataset), scaled_effect = {
        effect_sd <- stats::sd(.data$effect, na.rm = TRUE)
        if (dplyr::n() > 1L && is.finite(effect_sd) && effect_sd > 0) {
            as.numeric(scale(.data$effect))
        }
        else {
            rep(0, dplyr::n())
        }
    }, sig_label = make_sig_stars(.data$FDR)))
    if (!is.null(top_n_total)) {
        keep_metabs <- dplyr::pull(dplyr::slice_head(dplyr::arrange(dplyr::summarise(dplyr::group_by(long_df, 
            .data$Metabolite), max_abs_effect = max(abs(.data$effect), na.rm = TRUE), .groups = "drop"), 
            dplyr::desc(.data$max_abs_effect), .data$Metabolite), n = top_n_total), .data$Metabolite)
        long_df <- dplyr::filter(long_df, .data$Metabolite %in% keep_metabs)
    }
    metabolite_order <- dplyr::arrange(dplyr::summarise(dplyr::group_by(long_df, .data$Metabolite), order_effect = mean(.data$scaled_effect, 
        na.rm = TRUE), .groups = "drop"), dplyr::desc(.data$order_effect), .data$Metabolite)
    metabolite_levels <- rev(metabolite_order$Metabolite)
    grid <- dplyr::mutate(dplyr::left_join(tidyr::expand_grid(Metabolite = metabolite_levels, Dataset = dataset_order), 
        long_df, by = c("Metabolite", "Dataset")), Dataset = factor(.data$Dataset, levels = dataset_order), 
        Metabolite = factor(.data$Metabolite, levels = metabolite_levels))
    ggplot2::ggplot(grid, ggplot2::aes(x = .data$Dataset, y = .data$Metabolite, fill = .data$scaled_effect)) + 
        ggplot2::geom_tile(width = tile_width, height = 0.88, colour = NA, linewidth = 0) + ggplot2::geom_text(ggplot2::aes(label = .data$sig_label), 
        size = 2.8, na.rm = TRUE) + ggplot2::scale_fill_gradient2(low = "#3B4CC0", mid = "#F3F3F3", high = "#D7191C", 
        midpoint = 0, na.value = "#D9D9D9", name = "Scaled\neffect") + ggplot2::scale_x_discrete(expand = ggplot2::expansion(add = 0.55)) + 
        ggplot2::labs(x = NULL, y = y_title) + base_theme + ggplot2::theme(panel.grid = ggplot2::element_blank(), 
        panel.border = ggplot2::element_blank(), axis.line = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(), 
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1, size = 7.3), axis.text.y = ggplot2::element_text(size = metabolite_font_size), 
        axis.title.y = ggplot2::element_text(size = 8.2), legend.position = "none")
}
exvivo_nmr <- build_exvivo_nmr_matrix(file_nmr_exvivo_raw, file_nmr_exvivo_map)
panel_A <- make_pca_panel(sample_feature_df = exvivo_nmr$matrix, metadata = exvivo_nmr$metadata, colour_col = "Group", 
    shape_col = "Day", panel_label_text = "NMR Ex vivo", permanova_rhs = "Group", annotation_position = "bottom_right")
hilic_exvivo <- load_hilic_exvivo_for_pca(file_hilic_exvivo_rds)
panel_B <- make_pca_panel(sample_feature_df = hilic_exvivo$matrix, metadata = hilic_exvivo$metadata, 
    colour_col = "Group", shape_col = if ("Day" %in% names(hilic_exvivo$metadata)) "Day" else NULL, panel_label_text = "HILIC Ex vivo", 
    permanova_rhs = "Group", annotation_position = "top_right")
cecum_nmr <- load_cecum_nmr_from_results_workbook(file_cecum_nmr_workbook)
panel_C <- make_pca_panel(sample_feature_df = cecum_nmr$matrix, metadata = cecum_nmr$metadata, colour_col = "Group", 
    shape_col = NULL, panel_label_text = "NMR Caecum", permanova_rhs = "Group", annotation_position = "bottom_right")
cf_lcms <- load_cfeces_lcms_for_pca(file_cf_lcms_rds, representation = "one_per_metabolite")
panel_D <- make_pca_panel(sample_feature_df = cf_lcms$matrix, metadata = cf_lcms$metadata, colour_col = "Group", 
    shape_col = "Organ", panel_label_text = "LCMS Caecum/Faeces", permanova_rhs = "Group", annotation_position = "top_left")
serum <- load_serum_matrix(file_serum_rds)
panel_E <- make_pca_panel(sample_feature_df = serum$matrix, metadata = serum$metadata, colour_col = "Group", 
    shape_col = NULL, panel_label_text = "LCMS Serum", permanova_rhs = "Group", annotation_position = "bottom_right")
res_nmr_exvivo <- safe_read_csv(file_nmr_exvivo_results)
res_hilic_exvivo <- safe_read_csv(file_hilic_exvivo_results)
long_F <- dplyr::bind_rows(prepare_heatmap_long(res_hilic_exvivo, metabolite_col = "metabolite", effect_col = "log2FC", 
    fdr_col = "FDR", dataset_label = "HILIC (Ex vivo)", mode = heatmap_mode_exvivo, p_col = "p_value", 
    significance_source = "Donor-blocked limma + duplicateCorrelation; displayed when nominal P < 0.05"), 
    prepare_heatmap_long(res_nmr_exvivo, metabolite_col = "group_label", effect_col = "score_difference", 
        fdr_col = "FDR", dataset_label = "NMR (Ex vivo)", mode = heatmap_mode_exvivo, p_col = "p_value", 
        significance_source = "Donor-blocked limma + duplicateCorrelation; displayed when nominal P < 0.05"))
cf_cecum_res <- safe_read_csv(file_cf_cecum_results)
cf_feces_res <- safe_read_csv(file_cf_feces_results)
cf_combined_res <- dplyr::transmute(dplyr::filter(readxl::read_xlsx(ppi_input(file_cf_combined_results), 
    sheet = 1, .name_repair = "minimal"), .data$representation == "one_per_metabolite", .data$model == 
    "all_samples_donor_blocked_sensitivity", .data$contrast == "PPI_vs_control_adjusted_for_compartment"), 
    analyte_standardized = standardize_metabolite_label(.data$analyte), combined_log2_FC = as.numeric(.data$log2_FC), 
    combined_p_value = as.numeric(.data$p_value), combined_FDR = as.numeric(.data$FDR))
combined_duplicates <- dplyr::filter(dplyr::count(cf_combined_res, .data$analyte_standardized), .data$n > 
    1L)
if (nrow(combined_duplicates) > 0L) {
    stop("Duplicate standardized analytes were found in the combined caecum/faeces results:\n", paste(utils::capture.output(print(combined_duplicates)), 
        collapse = "\n"))
}
attach_combined_significance <- function(organ_results, organ_label) {
    out <- dplyr::left_join(dplyr::mutate(organ_results, analyte_standardized = standardize_metabolite_label(.data$analyte), 
        organ_specific_FDR = as.numeric(.data$FDR)), cf_combined_res, by = "analyte_standardized")
    missing_combined <- dplyr::filter(out, is.na(.data$combined_FDR))
    if (nrow(missing_combined) > 0L) {
        warning(nrow(missing_combined), " ", organ_label, " analytes could not be matched to the combined-compartment results.")
    }
    out
}
cf_cecum_res <- attach_combined_significance(cf_cecum_res, "caecum")
cf_feces_res <- attach_combined_significance(cf_feces_res, "faeces")
res_cecum_nmr <- readxl::read_xlsx(ppi_input(file_cecum_nmr_workbook), sheet = "Mouse_dupCor_sensitivity", 
    .name_repair = "minimal")
long_G <- dplyr::bind_rows(prepare_heatmap_long(cf_cecum_res, metabolite_col = "analyte_standardized", 
    effect_col = "log2_FC", fdr_col = "combined_FDR", dataset_label = "LCMS (Caecum)", mode = heatmap_mode_intestinal, 
    top_n = intestinal_top_n_per_dataset, p_col = "combined_p_value", significance_source = "Combined caecum/faeces model adjusted for compartment"), 
    prepare_heatmap_long(cf_feces_res, metabolite_col = "analyte_standardized", effect_col = "log2_FC", 
        fdr_col = "combined_FDR", dataset_label = "LCMS (Faeces)", mode = heatmap_mode_intestinal, top_n = intestinal_top_n_per_dataset, 
        p_col = "combined_p_value", significance_source = "Combined caecum/faeces model adjusted for compartment"), 
    prepare_heatmap_long(res_cecum_nmr, metabolite_col = "analysis_group", effect_col = "score_difference", 
        fdr_col = "FDR", dataset_label = "NMR (Caecum)", mode = heatmap_mode_intestinal, top_n = intestinal_top_n_per_dataset, 
        p_col = "p_value", significance_source = "Mouse-level limma + duplicateCorrelation"))
res_serum <- safe_read_csv(file_serum_results)
long_H <- prepare_heatmap_long(res_serum, metabolite_col = "metabolite", effect_col = "log2_FC", fdr_col = "FDR", 
    dataset_label = "LCMS (Serum)", mode = heatmap_mode_serum, top_n = serum_top_n, p_col = "p_value", 
    significance_source = "Mouse-level limma + duplicateCorrelation")
if (heatmap_mode_serum == "significant_only") {
    long_H <- dplyr::slice_head(dplyr::arrange(long_H, dplyr::desc(abs(.data$effect)), .data$FDR, .data$Metabolite), 
        n = serum_top_n)
}
panel_F <- build_heatmap_panel(long_df = long_F, dataset_order = c("HILIC (Ex vivo)", "NMR (Ex vivo)"), 
    top_n_total = NULL, y_title = "Metabolite", tile_width = 0.72, metabolite_font_size = 8.3)
panel_G <- build_heatmap_panel(long_df = long_G, dataset_order = c("LCMS (Caecum)", "LCMS (Faeces)", 
    "NMR (Caecum)"), top_n_total = NULL, y_title = NULL, tile_width = 0.72, metabolite_font_size = 8.1)
panel_H <- build_heatmap_panel(long_df = long_H, dataset_order = c("LCMS (Serum)"), top_n_total = serum_top_n, 
    y_title = NULL, tile_width = 0.42, metabolite_font_size = 8.5)
legend_plot <- build_heatmap_panel(long_df = long_H, dataset_order = c("LCMS (Serum)"), top_n_total = min(nrow(long_H), 
    serum_top_n), y_title = NULL, tile_width = 0.42, metabolite_font_size = 8.5) + ggplot2::theme(legend.position = "bottom")
heatmap_legend <- cowplot::get_legend(legend_plot)
pca_top_left <- cowplot::plot_grid(panel_A, panel_C, nrow = 1, rel_widths = c(1, 1), labels = c("A", 
    "C"), label_size = 14, label_fontface = "bold", label_x = 0, label_y = 1, hjust = 0, vjust = 1)
pca_bottom_left <- cowplot::plot_grid(panel_B, panel_D, nrow = 1, rel_widths = c(1, 1), labels = c("B", 
    "D"), label_size = 14, label_fontface = "bold", label_x = 0, label_y = 1, hjust = 0, vjust = 1)
pca_left_block <- cowplot::plot_grid(pca_top_left, pca_bottom_left, ncol = 1, rel_heights = c(1, 1))
pca_block <- cowplot::plot_grid(pca_left_block, panel_E, nrow = 1, rel_widths = c(2, 1), labels = c("", 
    "E"), label_size = 14, label_fontface = "bold", label_x = 0, label_y = 1, hjust = 0, vjust = 1)
panel_H_compact <- cowplot::plot_grid(NULL, panel_H, NULL, nrow = 1, rel_widths = c(0.22, 0.56, 0.22))
heatmap_row <- cowplot::plot_grid(panel_F, panel_G, panel_H_compact, nrow = 1, rel_widths = c(1, 1, 1), 
    labels = c("F", "G", "H"), label_size = 14, label_fontface = "bold", label_x = 0, label_y = 1, hjust = 0, 
    vjust = 1)
final_figure <- cowplot::plot_grid(pca_block, heatmap_row, heatmap_legend, ncol = 1, rel_heights = c(1.82, 
    1.42, 0.12))
png_file <- file.path(output_dir, "Figure3_multimodal_metabolomics_mouse_dupCor_v8_heatmaps_no_frame_wider_H_600dpi.png")
pdf_file <- file.path(output_dir, "Figure3_multimodal_metabolomics_mouse_dupCor_v8_heatmaps_no_frame_wider_H_600dpi.pdf")
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
    family = "Helvetica", useDingbats = FALSE, paper = "special")
print(final_figure)
grDevices::dev.off()
readr::write_csv(exvivo_nmr$metadata, ppi_output(file.path(output_dir, "panel_A_metadata.csv"), name = "panel_a_metadata.csv"))
readr::write_csv(dplyr::bind_cols(exvivo_nmr$metadata, exvivo_nmr$matrix), ppi_output(file.path(output_dir, 
    "panel_A_plot_matrix.csv"), name = "panel_a_plot_matrix.csv"))
readr::write_csv(hilic_exvivo$metadata, ppi_output(file.path(output_dir, "panel_B_metadata.csv"), name = "panel_b_metadata.csv"))
readr::write_csv(dplyr::bind_cols(hilic_exvivo$metadata, hilic_exvivo$matrix), ppi_output(file.path(output_dir, 
    "panel_B_plot_matrix.csv"), name = "panel_b_plot_matrix.csv"))
readr::write_csv(cecum_nmr$metadata, ppi_output(file.path(output_dir, "panel_C_metadata.csv"), name = "panel_c_metadata.csv"))
readr::write_csv(dplyr::bind_cols(cecum_nmr$metadata, cecum_nmr$matrix), ppi_output(file.path(output_dir, 
    "panel_C_plot_matrix.csv"), name = "panel_c_plot_matrix.csv"))
readr::write_csv(cf_lcms$metadata, ppi_output(file.path(output_dir, "panel_D_metadata.csv"), name = "panel_d_metadata.csv"))
readr::write_csv(dplyr::bind_cols(cf_lcms$metadata, cf_lcms$matrix), ppi_output(file.path(output_dir, 
    "panel_D_plot_matrix.csv"), name = "panel_d_plot_matrix.csv"))
readr::write_csv(serum$metadata, ppi_output(file.path(output_dir, "panel_E_metadata.csv"), name = "panel_e_metadata.csv"))
readr::write_csv(dplyr::bind_cols(serum$metadata, serum$matrix), ppi_output(file.path(output_dir, "panel_E_plot_matrix.csv"), 
    name = "panel_e_plot_matrix.csv"))
readr::write_csv(long_F, ppi_output(file.path(output_dir, "panel_F_heatmap_long.csv"), name = "panel_f_heatmap_long.csv"))
readr::write_csv(long_G, ppi_output(file.path(output_dir, "panel_G_heatmap_long.csv"), name = "panel_g_heatmap_long.csv"))
readr::write_csv(long_H, ppi_output(file.path(output_dir, "panel_H_heatmap_long.csv"), name = "panel_h_heatmap_long.csv"))
readr::write_csv(cf_combined_res, ppi_output(file.path(output_dir, "panel_G_combined_compartment_significance.csv"), 
    name = "panel_g_combined_compartment_significance.csv"))
readr::write_csv(cf_cecum_res, ppi_output(file.path(output_dir, "panel_G_caecum_mouse_effects_with_combined_FDR.csv"), 
    name = "panel_g_caecum_mouse_effects_with_combined_fdr.csv"))
readr::write_csv(cf_feces_res, ppi_output(file.path(output_dir, "panel_G_faeces_mouse_effects_with_combined_FDR.csv"), 
    name = "panel_g_faeces_mouse_effects_with_combined_fdr.csv"))
input_path_audit <- dplyr::mutate(tibble::tribble(~input, ~path, "configuration", config_file, "ex vivo NMR raw", 
    file_nmr_exvivo_raw, "ex vivo NMR mapping", file_nmr_exvivo_map, "ex vivo NMR results", file_nmr_exvivo_results, 
    "ex vivo HILIC RDS", file_hilic_exvivo_rds, "ex vivo HILIC results", file_hilic_exvivo_results, "caecum NMR workbook", 
    file_cecum_nmr_workbook, "caecum/faeces LC-MS RDS", file_cf_lcms_rds, "caecum LC-MS results", file_cf_cecum_results, 
    "faeces LC-MS results", file_cf_feces_results, "combined caecum/faeces LC-MS results", file_cf_combined_results, 
    "serum LC-MS RDS", file_serum_rds, "serum LC-MS results", file_serum_results), exists = ppi_file_exists(.data$path))
readr::write_csv(input_path_audit, ppi_output(file.path(output_dir, "Figure3_input_path_audit.csv"), 
    name = "figure3_input_path_audit.csv"))
message("Figure exported to:\n", png_file, "\n", pdf_file, "\n\nAll primary input paths were taken from: ", 
    config_file)
