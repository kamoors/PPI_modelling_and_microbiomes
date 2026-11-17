# SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0
# PCA preparation

source("scripts/helpers/setup.R")
invisible(run_analysis(
  section = "04_mouse_metabolomics",
  directory = "metabolomics",
  script = "04_10_prepare_metabolomics_ordination.R",
  restricted = FALSE,
  code = {
suppressPackageStartupMessages({
    library(dplyr)
    library(tidyr)
    library(readr)
    library(readxl)
    library(stringr)
    library(tibble)
    library(purrr)
    library(ggplot2)
    library(ggrepel)
    library(patchwork)
    library(scales)
})
set.seed(20260714)
message("Using PCA build script v3 bugfix")
project_dir <- project_path("")
metabolomics_root <- file.path(project_dir, "metabolomics")
input_paths <- list(serum_lcms_rds = file.path(metabolomics_root, "Serum_metabolomics", "20260121_data",
    "Serum_cleaned_after_collaborator_review", "Serum_LCMS_limma_input_objects.rds"), cecum_feces_lcms_rds = file.path(metabolomics_root,
    "20250714_MassSpecData", "Cecum_Feces_cleaned_after_collaborator_review", "Cecum_Feces_LCMS_limma_input_objects.rds"),
    cecum_nmr_results_xlsx = file.path(metabolomics_root, "NMR_data", "NMR_cecum_annotated_limma", "NMR_cecum_annotated_limma_results.xlsx"),
    exvivo_nmr_raw_csv = file.path(metabolomics_root, "NMR_data", "NMR_NEW_raw_data.csv"), exvivo_nmr_mapping_csv = file.path(metabolomics_root,
        "NMR_data", "NMR_exvivo_annotated_limma_dupCor", "NMR_exvivo_cleaned_feature_mapping_v2.csv"),
    exvivo_hilic_rds = file.path(metabolomics_root, "20250714_MassSpecData", "Ex_vivo_cleaned_one_per_metabolite",
        "HILIC_exvivo_cleaned_analysis_objects.rds"))
out_dir <- file.path(metabolomics_root, "Metabolomics_PCA")
dir.create(project_output(out_dir, directory = TRUE), recursive = TRUE, showWarnings = FALSE)
run_dataset <- list(serum_lcms = TRUE, intestinal_lcms = TRUE, cecum_nmr = TRUE, exvivo_nmr = TRUE, exvivo_hilic = TRUE)
enabled_input_keys <- c(if (run_dataset$serum_lcms) "serum_lcms_rds", if (run_dataset$intestinal_lcms) "cecum_feces_lcms_rds",
    if (run_dataset$cecum_nmr) "cecum_nmr_results_xlsx", if (run_dataset$exvivo_nmr) c("exvivo_nmr_raw_csv",
        "exvivo_nmr_mapping_csv"), if (run_dataset$exvivo_hilic) "exvivo_hilic_rds")
missing_input_keys <- enabled_input_keys[!vapply(input_paths[enabled_input_keys], project_file_exists, logical(1))]
if (length(missing_input_keys) > 0L) {
    missing_lines <- paste0("  - ", missing_input_keys, ": ", unlist(input_paths[missing_input_keys],
        use.names = FALSE))
    stop("The following enabled PCA input files were not found:\n", paste(missing_lines, collapse = "\n"),
        "\n\nCorrect the paths above or set the corresponding run_dataset flag to FALSE.")
}
make_separate_intestinal_organ_pcas <- TRUE
minimum_detected_fraction <- list(serum_lcms = 0.5, intestinal_lcms = 0.5, exvivo_hilic = 0.5)
ppi_colours <- c(ctrl = "#FF6E00", ppi = "#08306B")
organ_shapes <- c(cecum = 16, feces = 17)
output_width <- 7.2
output_height <- 5.8
output_dpi <- 450
label_donors_on_static_pca <- TRUE
`%||%` <- function(x, y) {
    if (is.null(x) || length(x) == 0L)
        y
    else x
}
first_existing_column <- function(data, candidates, required = TRUE) {
    hit <- candidates[candidates %in% names(data)]
    if (length(hit) > 0L) {
        return(hit[[1]])
    }
    if (required) {
        stop("None of these columns were found: ", paste(candidates, collapse = ", "), ". Available columns: ",
            paste(names(data), collapse = ", "))
    }
    NULL
}
normalize_group <- function(x) {
    x <- str_to_lower(str_squish(as.character(x)))
    dplyr::case_when(x %in% c("ctrl", "control", "ctlr", "c", "0", "no") ~ "ctrl", x %in% c("ppi", "p",
        "1", "yes") ~ "ppi", str_detect(x, "^ctrl|^control|^ctlr|^c[0-9]+$") ~ "ctrl", str_detect(x,
        "^ppi|^p[0-9]+$") ~ "ppi", TRUE ~ NA_character_)
}
normalize_donor <- function(x, sample_id = NULL) {
    x <- as.character(x)
    if (!is.null(sample_id)) {
        x[is.na(x) | x == ""] <- as.character(sample_id)[is.na(x) | x == ""]
    }
    raw <- str_to_lower(x)
    donor <- str_extract(raw, "(?:ctrl|ctlr|control|ppi|[cp])_?[0-9]+")
    donor <- str_replace_all(donor, "control|ctlr", "ctrl")
    donor <- str_replace_all(donor, "_", "")
    donor <- case_when(str_detect(donor, "^c[0-9]+$") ~ str_replace(donor, "^c", "ctrl"), str_detect(donor,
        "^p[0-9]+$") ~ str_replace(donor, "^p", "ppi"), TRUE ~ donor)
    donor
}
core_donor_code <- function(donor) {
    donor <- normalize_donor(donor)
    case_when(str_detect(donor, "^ctrl") ~ paste0("C", str_extract(donor, "[0-9]+")), str_detect(donor,
        "^ppi") ~ paste0("P", str_extract(donor, "[0-9]+")), TRUE ~ NA_character_)
}
as_numeric_matrix <- function(x) {
    x <- as.data.frame(x, check.names = FALSE)
    out <- vapply(x, function(z) suppressWarnings(as.numeric(z)), numeric(nrow(x)))
    out <- as.matrix(out)
    rownames(out) <- rownames(x)
    colnames(out) <- colnames(x)
    out
}
align_matrix_and_metadata <- function(X, metadata, sample_column = NULL) {
    X <- as_numeric_matrix(X)
    metadata <- as.data.frame(metadata, check.names = FALSE)
    if (is.null(sample_column)) {
        sample_column <- first_existing_column(metadata, c("Sample_ID", "Sample_Name", "sample_id", "sample",
            "ID"), required = FALSE)
    }
    if (!is.null(sample_column)) {
        sample_ids <- as.character(metadata[[sample_column]])
        if (!is.null(rownames(X)) && all(rownames(X) %in% sample_ids)) {
            metadata <- metadata[match(rownames(X), sample_ids), , drop = FALSE]
        }
        else if (nrow(X) == nrow(metadata)) {
            rownames(X) <- sample_ids
        }
        else {
            stop("Could not align feature matrix and metadata by sample ID.")
        }
    }
    else if (nrow(X) != nrow(metadata)) {
        stop("Matrix and metadata row counts differ and no sample column was found.")
    }
    if (is.null(rownames(X))) {
        rownames(X) <- paste0("sample_", seq_len(nrow(X)))
    }
    metadata$Sample_ID <- rownames(X)
    rownames(metadata) <- metadata$Sample_ID
    list(X = X, metadata = metadata)
}
calculate_half_minimum <- function(x) {
    positive <- x[is.finite(x) & x > 0]
    if (length(positive) == 0L)
        return(NA_real_)
    min(positive)/2
}
median_impute_columns <- function(X) {
    for (j in seq_len(ncol(X))) {
        bad <- !is.finite(X[, j])
        if (any(bad)) {
            med <- median(X[, j][!bad], na.rm = TRUE)
            if (!is.finite(med))
                med <- 0
            X[bad, j] <- med
        }
    }
    X
}
prepare_lcms_for_pca <- function(X, minimum_fraction = 0.5, label = "LCMS") {
    X <- as_numeric_matrix(X)
    if (is.null(colnames(X))) {
        colnames(X) <- paste0("feature_", seq_len(ncol(X)))
    }
    if (anyDuplicated(colnames(X))) {
        stop(label, ": duplicated metabolite column names were found. ", "The PCA input must contain one uniquely named column per metabolite.")
    }
    X[!is.finite(X)] <- NA_real_
    X[X < 0] <- 0
    detected_fraction <- as.numeric(colMeans(X > 0, na.rm = TRUE))
    variance_raw <- as.numeric(apply(X, 2, var, na.rm = TRUE))
    keep_initial <- detected_fraction >= minimum_fraction & is.finite(variance_raw) & variance_raw >
        0
    if (!any(keep_initial)) {
        stop(label, ": no features passed PCA prevalence/variance filtering.")
    }
    X <- X[, keep_initial, drop = FALSE]
    detected_fraction <- detected_fraction[keep_initial]
    pseudocounts <- vapply(seq_len(ncol(X)), function(j) calculate_half_minimum(X[, j]), numeric(1))
    keep_pc <- is.finite(pseudocounts) & pseudocounts > 0
    if (!any(keep_pc)) {
        stop(label, ": no retained features had a usable positive-value pseudocount.")
    }
    X <- X[, keep_pc, drop = FALSE]
    detected_fraction <- detected_fraction[keep_pc]
    pseudocounts <- pseudocounts[keep_pc]
    audit_names <- colnames(X)
    for (j in seq_len(ncol(X))) {
        X[, j] <- log2(X[, j] + pseudocounts[j])
    }
    X <- median_impute_columns(X)
    variance_log <- as.numeric(apply(X, 2, var, na.rm = TRUE))
    keep_log <- is.finite(variance_log) & variance_log > 0
    audit <- tibble(metabolite = audit_names, detected_fraction = detected_fraction, pseudocount = pseudocounts,
        retained_after_log_variance = keep_log)
    X <- X[, keep_log, drop = FALSE]
    if (ncol(X) < 2L) {
        stop(label, ": fewer than two metabolites remained after log-variance filtering.")
    }
    list(X = X, audit = audit)
}
prepare_score_matrix_for_pca <- function(X, label = "score matrix") {
    X <- as_numeric_matrix(X)
    X <- median_impute_columns(X)
    variance_x <- apply(X, 2, var, na.rm = TRUE)
    keep <- is.finite(variance_x) & variance_x > 0
    if (!any(keep))
        stop(label, ": no non-constant columns available for PCA.")
    list(X = X[, keep, drop = FALSE], audit = tibble(metabolite = colnames(X), retained_after_variance = keep))
}
run_pca <- function(dataset_id, X, metadata, input_type, audit = NULL) {
    if (nrow(X) != nrow(metadata))
        stop(dataset_id, ": matrix/metadata mismatch.")
    if (ncol(X) < 2L)
        stop(dataset_id, ": fewer than two metabolites available.")
    pca <- prcomp(X, center = TRUE, scale. = TRUE)
    explained <- (pca$sdev^2)/sum(pca$sdev^2)
    scores <- as_tibble(pca$x, rownames = "Sample_ID") %>% left_join(metadata %>% mutate(Sample_ID = as.character(Sample_ID)),
        by = "Sample_ID")
    loadings <- as_tibble(pca$rotation, rownames = "Metabolite")
    variance <- tibble(PC = paste0("PC", seq_along(explained)), variance_explained = explained, cumulative_variance = cumsum(explained))
    list(dataset_id = dataset_id, input_type = input_type, pca = pca, X_processed = X, metadata = metadata,
        scores = scores, loadings = loadings, variance = variance, preprocessing_audit = audit)
}
make_exvivo_delta <- function(X_transformed, metadata) {
    metadata <- metadata %>% mutate(Day_num = as.integer(Day_num))
    donors <- unique(metadata$Donor)
    delta_rows <- list()
    delta_meta <- list()
    for (donor in donors) {
        idx <- which(metadata$Donor == donor)
        baseline_idx <- idx[metadata$Day_num[idx] == 0]
        follow_idx <- idx[metadata$Day_num[idx] %in% 1:4]
        if (length(baseline_idx) != 1L || length(follow_idx) == 0L) {
            warning("Skipping donor ", donor, ": incomplete day-0/follow-up structure.")
            next
        }
        baseline <- X_transformed[baseline_idx, , drop = FALSE]
        delta <- sweep(X_transformed[follow_idx, , drop = FALSE], 2, as.numeric(baseline), FUN = "-")
        delta_rows[[donor]] <- delta
        delta_meta[[donor]] <- metadata[follow_idx, , drop = FALSE]
    }
    X_delta <- do.call(rbind, delta_rows)
    meta_delta <- bind_rows(delta_meta)
    rownames(X_delta) <- meta_delta$Sample_ID
    list(X = X_delta, metadata = meta_delta)
}
variance_label <- function(pca_object, pc) {
    value <- pca_object$variance$variance_explained[match(pc, pca_object$variance$PC)]
    paste0(pc, " (", percent(value, accuracy = 0.1), ")")
}
base_pca_theme <- function() {
    theme_classic(base_size = 11) + theme(legend.title = element_text(face = "bold"), axis.title = element_text(face = "bold"),
        plot.title = element_text(face = "bold", size = 12), plot.subtitle = element_text(size = 10),
        strip.background = element_blank(), strip.text = element_text(face = "bold"))
}
plot_static_pca <- function(obj, title, shape_variable = NULL, label_variable = NULL, facet_variable = NULL) {
    dat <- obj$scores
    if (!is.null(shape_variable) && shape_variable %in% names(dat)) {
        p <- ggplot(dat, aes(x = PC1, y = PC2, colour = Group, shape = .data[[shape_variable]]))
    }
    else {
        p <- ggplot(dat, aes(x = PC1, y = PC2, colour = Group))
    }
    p <- p + geom_hline(yintercept = 0, linewidth = 0.25, colour = "grey85") + geom_vline(xintercept = 0,
        linewidth = 0.25, colour = "grey85") + geom_point(size = 3.1, alpha = 0.9) + scale_colour_manual(values = ppi_colours,
        drop = FALSE) + labs(title = title, x = variance_label(obj, "PC1"), y = variance_label(obj, "PC2"),
        colour = "Inoculum") + base_pca_theme()
    if (!is.null(shape_variable) && shape_variable == "Organ") {
        p <- p + scale_shape_manual(values = organ_shapes, drop = FALSE)
    }
    if (!is.null(label_variable) && label_variable %in% names(dat)) {
        p <- p + ggrepel::geom_text_repel(aes(label = .data[[label_variable]]), size = 2.7, show.legend = FALSE,
            max.overlaps = Inf, box.padding = 0.25, point.padding = 0.15)
    }
    if (!is.null(facet_variable) && facet_variable %in% names(dat)) {
        p <- p + facet_wrap(stats::as.formula(paste("~", facet_variable)))
    }
    p
}
plot_exvivo_trajectory <- function(obj, title) {
    dat <- obj$scores %>% arrange(Donor, Day_num)
    ggplot(dat, aes(PC1, PC2, colour = Group, group = Donor)) + geom_hline(yintercept = 0, linewidth = 0.25,
        colour = "grey85") + geom_vline(xintercept = 0, linewidth = 0.25, colour = "grey85") + geom_path(linewidth = 0.8,
        alpha = 0.75, arrow = grid::arrow(length = grid::unit(0.11, "inches"), type = "closed")) + geom_point(aes(shape = DayF),
        size = 3, alpha = 0.95) + ggrepel::geom_text_repel(data = dat %>% group_by(Donor) %>% filter(Day_num ==
        max(Day_num)) %>% ungroup(), aes(label = Donor), size = 2.8, show.legend = FALSE, max.overlaps = Inf) +
        scale_colour_manual(values = ppi_colours, drop = FALSE) + labs(title = title, subtitle = "Lines connect repeated measurements from the same donor culture",
        x = variance_label(obj, "PC1"), y = variance_label(obj, "PC2"), colour = "Inoculum", shape = "Day") +
        base_pca_theme()
}
plot_paired_compartments <- function(obj, title) {
    dat <- obj$scores
    ggplot(dat, aes(PC1, PC2, colour = Group)) + geom_hline(yintercept = 0, linewidth = 0.25, colour = "grey85") +
        geom_vline(xintercept = 0, linewidth = 0.25, colour = "grey85") + geom_line(aes(group = Pair_ID),
        linewidth = 0.45, alpha = 0.35, colour = "grey45") + geom_point(aes(shape = Organ), size = 3,
        alpha = 0.92) + scale_colour_manual(values = ppi_colours, drop = FALSE) + scale_shape_manual(values = organ_shapes,
        drop = FALSE) + labs(title = title, subtitle = "Lines connect caecum and faeces from the same recipient mouse",
        x = variance_label(obj, "PC1"), y = variance_label(obj, "PC2"), colour = "Inoculum", shape = "Compartment") +
        base_pca_theme()
}
plot_scree <- function(obj, title, n_pcs = 10L) {
    n_keep <- min(as.integer(n_pcs), nrow(obj$variance))
    dat <- obj$variance %>% slice_head(n = n_keep)
    ggplot(dat, aes(x = factor(PC, levels = PC), y = variance_explained)) + geom_col(width = 0.75, fill = "grey40") +
        geom_line(aes(group = 1), linewidth = 0.5) + geom_point(size = 1.8) + scale_y_continuous(labels = percent_format(accuracy = 1)) +
        labs(title = title, x = NULL, y = "Variance explained") + base_pca_theme()
}
plot_loading_map <- function(obj, title, n_labels = 18L) {
    n_keep <- min(as.integer(n_labels), nrow(obj$loadings))
    dat <- obj$loadings %>% mutate(importance = sqrt(PC1^2 + PC2^2)) %>% slice_max(order_by = importance,
        n = n_keep, with_ties = FALSE)
    ggplot(dat, aes(PC1, PC2)) + geom_hline(yintercept = 0, linewidth = 0.25, colour = "grey80") + geom_vline(xintercept = 0,
        linewidth = 0.25, colour = "grey80") + geom_segment(aes(x = 0, y = 0, xend = PC1, yend = PC2),
        arrow = grid::arrow(length = grid::unit(0.08, "inches")), linewidth = 0.45, colour = "grey35") +
        ggrepel::geom_text_repel(aes(label = Metabolite), size = 2.7, max.overlaps = Inf) + labs(title = title,
        x = "PC1 loading", y = "PC2 loading") + base_pca_theme()
}
save_pca_object_outputs <- function(obj, plot_main, title_stub) {
    dataset_dir <- file.path(out_dir, obj$dataset_id)
    dir.create(project_output(dataset_dir, directory = TRUE), recursive = TRUE, showWarnings = FALSE)
    write_csv(obj$scores, project_output(file.path(dataset_dir, "PCA_scores.csv"), name = "pca_scores.csv"))
    write_csv(obj$loadings, project_output(file.path(dataset_dir, "PCA_loadings.csv"), name = "pca_loadings.csv"))
    write_csv(obj$variance, project_output(file.path(dataset_dir, "PCA_variance_explained.csv"), name = "pca_variance_explained.csv"))
    write_csv(as_tibble(obj$X_processed, rownames = "Sample_ID"), project_output(file.path(dataset_dir, "PCA_processed_matrix.csv"),
        name = "pca_processed_matrix.csv"))
    if (!is.null(obj$preprocessing_audit)) {
        write_csv(obj$preprocessing_audit, project_output(file.path(dataset_dir, "PCA_preprocessing_audit.csv"),
            name = "pca_preprocessing_audit.csv"))
    }
    scree <- plot_scree(obj, paste0(title_stub, ": scree plot"))
    loadings <- plot_loading_map(obj, paste0(title_stub, ": leading PC1/PC2 loadings"))
    ggsave(project_output(file.path(dataset_dir, "PCA_PC1_PC2.png"), name = "pca_pc1_pc2.png"), plot_main,
        width = output_width, height = output_height, dpi = output_dpi, bg = "white")
    ggsave(project_output(file.path(dataset_dir, "PCA_PC1_PC2.pdf"), name = "pca_pc1_pc2.pdf"), plot_main,
        width = output_width, height = output_height, device = grDevices::cairo_pdf)
    ggsave(project_output(file.path(dataset_dir, "PCA_scree.png"), name = "pca_scree.png"), scree, width = 6.5,
        height = 4.5, dpi = output_dpi, bg = "white")
    ggsave(project_output(file.path(dataset_dir, "PCA_loadings_PC1_PC2.png"), name = "pca_loadings_pc1_pc2.png"),
        loadings, width = 7.5, height = 6, dpi = output_dpi, bg = "white")
}
load_serum_lcms <- function(path) {
    obj <- readRDS(project_input(path))
    X <- obj$serum_clean_bio %||% obj$limma_feature_matrix
    meta <- obj$metadata_clean_bio %||% obj$limma_metadata
    if (is.null(X) || is.null(meta))
        stop("Serum RDS does not contain expected objects.")
    aligned <- align_matrix_and_metadata(X, meta)
    X <- aligned$X
    meta <- aligned$metadata
    group_col <- first_existing_column(meta, c("Separation", "group", "grp", "PPI"))
    donor_col <- first_existing_column(meta, c("Donor", "inoculum_donor"), required = FALSE)
    mouse_col <- first_existing_column(meta, c("Mouse_ID", "mouse", "sample_id"), required = FALSE)
    meta <- meta %>% mutate(Group = factor(normalize_group(.data[[group_col]]), levels = c("ctrl", "ppi")),
        Donor = normalize_donor(if (!is.null(donor_col))
            .data[[donor_col]]
        else NA_character_, Sample_ID), Mouse_ID = if (!is.null(mouse_col))
            as.character(.data[[mouse_col]])
        else Sample_ID, Pair_ID = paste(Donor, Mouse_ID, sep = "__"), Dataset = "Serum LC-MS")
    prep <- prepare_lcms_for_pca(X, minimum_fraction = minimum_detected_fraction$serum_lcms, label = "Serum LC-MS")
    run_pca("serum_lcms", prep$X, meta, "LC-MS log2", prep$audit)
}
load_intestinal_lcms <- function(path) {
    obj <- readRDS(project_input(path))
    if (is.null(obj$one_per_metabolite))
        stop("Intestinal RDS lacks one_per_metabolite.")
    X_c <- as_numeric_matrix(obj$one_per_metabolite$cecum)
    X_f <- as_numeric_matrix(obj$one_per_metabolite$feces)
    meta_c <- as.data.frame(obj$metadata$cecum)
    meta_f <- as.data.frame(obj$metadata$feces)
    common_metabolites <- intersect(colnames(X_c), colnames(X_f))
    if (length(common_metabolites) < 2L)
        stop("No shared caecum/faeces metabolites.")
    X_c <- X_c[, common_metabolites, drop = FALSE]
    X_f <- X_f[, common_metabolites, drop = FALSE]
    a_c <- align_matrix_and_metadata(X_c, meta_c)
    a_f <- align_matrix_and_metadata(X_f, meta_f)
    a_c$metadata$Organ <- "cecum"
    a_f$metadata$Organ <- "feces"
    X <- rbind(a_c$X, a_f$X)
    meta <- bind_rows(a_c$metadata, a_f$metadata)
    rownames(X) <- meta$Sample_ID
    group_col <- first_existing_column(meta, c("Separation", "group", "grp"))
    donor_col <- first_existing_column(meta, c("Donor", "inoculum_donor"), required = FALSE)
    mouse_col <- first_existing_column(meta, c("Mouse_ID", "mouse"), required = FALSE)
    pair_col <- first_existing_column(meta, c("Donor_Mouse_ID", "Pair_ID"), required = FALSE)
    meta <- meta %>% mutate(Group = factor(normalize_group(.data[[group_col]]), levels = c("ctrl", "ppi")),
        Donor = normalize_donor(if (!is.null(donor_col))
            .data[[donor_col]]
        else NA_character_, Sample_ID), Mouse_ID = if (!is.null(mouse_col))
            as.character(.data[[mouse_col]])
        else str_extract(Sample_ID, "mouse[_-]?[0-9]+|[0-9]+"), Pair_ID = if (!is.null(pair_col))
            as.character(.data[[pair_col]])
        else paste(Donor, Mouse_ID, sep = "__"), Organ = factor(str_to_lower(Organ), levels = c("cecum",
            "feces")), Dataset = "Caecum/faeces LC-MS")
    prep <- prepare_lcms_for_pca(X, minimum_fraction = minimum_detected_fraction$intestinal_lcms, label = "Caecum/faeces LC-MS")
    combined <- run_pca("intestinal_lcms_combined", prep$X, meta, "LC-MS log2", prep$audit)
    separate <- list()
    if (make_separate_intestinal_organ_pcas) {
        for (organ in c("cecum", "feces")) {
            idx <- meta$Organ == organ
            separate[[organ]] <- run_pca(paste0("intestinal_lcms_", organ), prep$X[idx, , drop = FALSE],
                meta[idx, , drop = FALSE], "LC-MS log2", prep$audit)
        }
    }
    list(combined = combined, separate = separate)
}
load_cecum_nmr <- function(path) {
    long <- read_excel(project_input(path), sheet = "Mouse_level_values") %>% mutate(sample_id = as.character(sample_id),
        analysis_group = as.character(analysis_group), score = as.numeric(score))
    meta <- long %>% distinct(sample_id, inoculum_donor, grp, group, isolator, age_weeks, weight_mg,
        mouse) %>% transmute(Sample_ID = sample_id, Group = factor(normalize_group(group), levels = c("ctrl",
        "ppi")), Donor = normalize_donor(inoculum_donor), Mouse_ID = as.character(mouse), Pair_ID = paste(Donor,
        Mouse_ID, sep = "__"), Isolator = isolator, Age_weeks = age_weeks, Weight_mg = weight_mg, Dataset = "Caecum NMR")
    wide <- long %>% select(sample_id, analysis_group, score) %>% pivot_wider(names_from = analysis_group,
        values_from = score)
    X <- wide %>% select(-sample_id) %>% as.data.frame(check.names = FALSE) %>% as.matrix()
    rownames(X) <- wide$sample_id
    meta <- meta[match(rownames(X), meta$Sample_ID), , drop = FALSE]
    prep <- prepare_score_matrix_for_pca(X, "Caecum NMR")
    run_pca("cecum_nmr", prep$X, meta, "standardized NMR group scores", prep$audit)
}
load_exvivo_nmr <- function(raw_path, mapping_path) {
    raw <- read_csv(project_input(raw_path), show_col_types = FALSE)
    mapping <- read_csv(project_input(mapping_path), show_col_types = FALSE)
    mapping_id_col <- first_existing_column(mapping, c("feature_id", "CLASSYID"))
    if (mapping_id_col != "CLASSYID") {
        mapping <- mapping %>% rename(CLASSYID = all_of(mapping_id_col))
    }
    mapping <- mapping %>% mutate(include_group_level = as.logical(include_group_level), CLASSYID = as.character(CLASSYID),
        analysis_group = as.character(analysis_group)) %>% filter(include_group_level, !is.na(analysis_group),
        analysis_group != "")
    sample_cols <- names(raw)[str_detect(names(raw), "^[0-4]_[CP][1-3]$")]
    if (length(sample_cols) != 30L) {
        warning("Expected 30 ex vivo NMR biological samples (6 donors x 5 days), found ", length(sample_cols))
    }
    feature_table <- raw %>% mutate(CLASSYID = as.character(CLASSYID)) %>% inner_join(mapping %>% select(CLASSYID,
        analysis_group), by = "CLASSYID")
    X_raw <- feature_table %>% select(all_of(sample_cols)) %>% as.data.frame(check.names = FALSE) %>%
        as.matrix()
    storage.mode(X_raw) <- "numeric"
    rownames(X_raw) <- feature_table$CLASSYID
    X_raw[!is.finite(X_raw)] <- NA_real_
    X_raw <- pmax(X_raw, 0)
    pseudocount <- apply(X_raw, 1, calculate_half_minimum)
    keep <- is.finite(pseudocount) & pseudocount > 0
    X_raw <- X_raw[keep, , drop = FALSE]
    pseudocount <- pseudocount[keep]
    feature_table <- feature_table[keep, , drop = FALSE]
    X_log <- X_raw
    for (i in seq_len(nrow(X_log))) {
        X_log[i, ] <- log2(X_raw[i, ] + pseudocount[[i]])
    }
    X_z <- t(scale(t(X_log)))
    X_z[!is.finite(X_z)] <- 0
    groups <- split(feature_table$CLASSYID, feature_table$analysis_group)
    group_matrix <- vapply(groups, function(ids) {
        ids <- intersect(ids, rownames(X_z))
        if (length(ids) == 1L) {
            as.numeric(X_z[ids, ])
        }
        else {
            apply(X_z[ids, , drop = FALSE], 2, median, na.rm = TRUE)
        }
    }, numeric(ncol(X_z)))
    X <- as.matrix(group_matrix)
    rownames(X) <- colnames(X_z)
    meta <- tibble(Sample_ID = rownames(X)) %>% extract(Sample_ID, into = c("Day_num", "Core_Donor"),
        regex = "^([0-4])_([CP][1-3])$", remove = FALSE, convert = TRUE) %>% mutate(Group = factor(if_else(str_starts(Core_Donor,
        "C"), "ctrl", "ppi"), levels = c("ctrl", "ppi")), Donor = normalize_donor(Core_Donor), DayF = factor(Day_num,
        levels = 0:4), Match_ID = paste0(Day_num, "_", Core_Donor), Dataset = "Ex vivo NMR")
    prep <- prepare_score_matrix_for_pca(X, "Ex vivo NMR")
    raw_pca <- run_pca("exvivo_nmr_raw", prep$X, meta, "standardized NMR group scores", prep$audit)
    delta <- make_exvivo_delta(prep$X, meta)
    delta_prep <- prepare_score_matrix_for_pca(delta$X, "Ex vivo NMR baseline-centred")
    delta_pca <- run_pca("exvivo_nmr_delta", delta_prep$X, delta$metadata, "day-0-centred NMR group scores",
        delta_prep$audit)
    list(raw = raw_pca, delta = delta_pca)
}
load_exvivo_hilic <- function(path) {
    obj <- readRDS(project_input(path))
    X <- obj$primary_matrix_all_samples
    meta <- as.data.frame(obj$metadata_all)
    if (is.null(X) || is.null(meta))
        stop("HILIC RDS lacks primary matrix or metadata.")
    aligned <- align_matrix_and_metadata(X, meta)
    X <- aligned$X
    meta <- aligned$metadata
    medium_col <- first_existing_column(meta, c("is_medium"), required = FALSE)
    keep <- if (!is.null(medium_col))
        !as.logical(meta[[medium_col]])
    else !str_detect(str_to_lower(meta$Sample_ID), "medium")
    X <- X[keep, , drop = FALSE]
    meta <- meta[keep, , drop = FALSE]
    group_col <- first_existing_column(meta, c("Separation", "group", "grp"))
    donor_col <- first_existing_column(meta, c("Donor", "donor"))
    day_col <- first_existing_column(meta, c("Day", "day"))
    meta <- meta %>% mutate(Group = factor(normalize_group(.data[[group_col]]), levels = c("ctrl", "ppi")),
        Donor = normalize_donor(.data[[donor_col]], Sample_ID), Day_num = as.integer(str_extract(as.character(.data[[day_col]]),
            "[0-9]+")), DayF = factor(Day_num, levels = 0:4), Core_Donor = core_donor_code(Donor), Match_ID = paste0(Day_num,
            "_", Core_Donor), Dataset = "Ex vivo HILIC LC-MS")
    prep <- prepare_lcms_for_pca(X, minimum_fraction = minimum_detected_fraction$exvivo_hilic, label = "Ex vivo HILIC")
    raw_pca <- run_pca("exvivo_hilic_raw", prep$X, meta, "HILIC LC-MS log2", prep$audit)
    delta <- make_exvivo_delta(prep$X, meta)
    delta_prep <- prepare_score_matrix_for_pca(delta$X, "Ex vivo HILIC baseline-centred")
    delta_pca <- run_pca("exvivo_hilic_delta", delta_prep$X, delta$metadata, "day-0-centred HILIC log2 abundance",
        delta_prep$audit)
    list(raw = raw_pca, delta = delta_pca)
}
pca_objects <- list()
if (run_dataset$serum_lcms) {
    message("Running serum LC-MS PCA...")
    obj <- load_serum_lcms(input_paths$serum_lcms_rds)
    pca_objects[[obj$dataset_id]] <- obj
    p <- plot_static_pca(obj, "Mouse serum LC-MS", label_variable = if (label_donors_on_static_pca)
        "Donor"
    else NULL)
    save_pca_object_outputs(obj, p, "Mouse serum LC-MS")
}
if (run_dataset$intestinal_lcms) {
    message("Running caecum/faeces LC-MS PCA...")
    intestinal <- load_intestinal_lcms(input_paths$cecum_feces_lcms_rds)
    obj <- intestinal$combined
    pca_objects[[obj$dataset_id]] <- obj
    p <- plot_paired_compartments(obj, "Mouse intestinal LC-MS")
    save_pca_object_outputs(obj, p, "Mouse intestinal LC-MS")
    for (organ in names(intestinal$separate)) {
        organ_obj <- intestinal$separate[[organ]]
        pca_objects[[organ_obj$dataset_id]] <- organ_obj
        p_organ <- plot_static_pca(organ_obj, paste0("Mouse ", organ, " LC-MS"), label_variable = if (label_donors_on_static_pca)
            "Donor"
        else NULL)
        save_pca_object_outputs(organ_obj, p_organ, paste0("Mouse ", organ, " LC-MS"))
    }
}
if (run_dataset$cecum_nmr) {
    message("Running caecum NMR PCA...")
    obj <- load_cecum_nmr(input_paths$cecum_nmr_results_xlsx)
    pca_objects[[obj$dataset_id]] <- obj
    p <- plot_static_pca(obj, "Mouse caecum NMR", label_variable = if (label_donors_on_static_pca)
        "Donor"
    else NULL)
    save_pca_object_outputs(obj, p, "Mouse caecum NMR")
}
if (run_dataset$exvivo_nmr) {
    message("Running ex vivo NMR PCAs...")
    ex_nmr <- load_exvivo_nmr(input_paths$exvivo_nmr_raw_csv, input_paths$exvivo_nmr_mapping_csv)
    for (kind in names(ex_nmr)) {
        obj <- ex_nmr[[kind]]
        pca_objects[[obj$dataset_id]] <- obj
        p <- plot_exvivo_trajectory(obj, if (kind == "raw") {
            "Ex vivo NMR metabolic trajectories"
        }
        else {
            "Ex vivo NMR trajectories relative to day 0"
        })
        save_pca_object_outputs(obj, p, if (kind == "raw")
            "Ex vivo NMR"
        else "Ex vivo NMR, day-0 centred")
    }
}
if (run_dataset$exvivo_hilic) {
    message("Running ex vivo HILIC PCAs...")
    ex_hilic <- load_exvivo_hilic(input_paths$exvivo_hilic_rds)
    for (kind in names(ex_hilic)) {
        obj <- ex_hilic[[kind]]
        pca_objects[[obj$dataset_id]] <- obj
        p <- plot_exvivo_trajectory(obj, if (kind == "raw") {
            "Ex vivo HILIC LC-MS metabolic trajectories"
        }
        else {
            "Ex vivo HILIC LC-MS trajectories relative to day 0"
        })
        save_pca_object_outputs(obj, p, if (kind == "raw")
            "Ex vivo HILIC LC-MS"
        else "Ex vivo HILIC LC-MS, day-0 centred")
    }
}
overview <- imap_dfr(pca_objects, function(obj, id) {
    tibble(dataset_id = id, n_samples = nrow(obj$X_processed), n_metabolites = ncol(obj$X_processed),
        PC1_variance = obj$variance$variance_explained[1], PC2_variance = obj$variance$variance_explained[2],
        PC1_PC2_cumulative = sum(obj$variance$variance_explained[1:2]), input_type = obj$input_type)
})
write_csv(overview, project_output(file.path(out_dir, "PCA_dataset_overview.csv"), name = "pca_dataset_overview.csv"))
saveRDS(list(pca_objects = pca_objects, overview = overview, settings = list(minimum_detected_fraction = minimum_detected_fraction,
    ppi_colours = ppi_colours, note = "Use PPI_metabolomics_PCA_02_statistical_analysis.R for donor-aware inference.")),
    project_output(file.path(out_dir, "PPI_metabolomics_PCA_bundle.rds"), name = "ppi_metabolomics_pca_bundle.rds"))
message("PCA plotting complete. Outputs written to: ", out_dir)
}))
