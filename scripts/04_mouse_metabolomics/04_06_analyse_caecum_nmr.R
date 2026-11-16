# SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0
# Caecum NMR
# Run through run.R so input and output paths resolve within the package.
# Inputs: see data/input_manifest.csv and docs/RUN_ORDER.md.
# Writes generated files under results/; bundled inputs are read-only.

required_packages <- c("dplyr", "tidyr", "stringr", "tibble", "readxl", "readr", "limma", "statmod", 
    "writexl", "ggplot2", "purrr")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]
if (length(missing_packages) > 0) {
    stop("Install these packages first: ", paste(missing_packages, collapse = ", "))
}
cecum_xlsx <- ppi_path("metabolomics/NMR_data/cecum_Karlis.xlsx")
out_dir <- ppi_path("metabolomics/NMR_data/NMR_cecum_annotated_limma")
dir.create(ppi_output(out_dir, directory = TRUE), recursive = TRUE, showWarnings = FALSE)
group_levels <- c("ctrl", "ppi")
run_adjusted_mouse_sensitivity <- T
clean_metabolite_name <- function(x) {
    x <- stringr::str_squish(as.character(x))
    x <- dplyr::na_if(x, "")
    xl <- stringr::str_to_lower(x)
    dplyr::case_when(is.na(x) ~ NA_character_, xl == "3-(3-oh-phenyl)propionate" ~ "3-hydroxyphenylpropionate", 
        xl == "5-aminopentanoic acid" ~ "5-aminopentanoate", xl == "phosphocoline" ~ "phosphocholine", 
        TRUE ~ xl)
}
make_feature_pseudocount <- function(x) {
    positive <- x[is.finite(x) & x > 0]
    if (length(positive) == 0) {
        return(NA_real_)
    }
    min(positive, na.rm = TRUE)/2
}
make_group_score <- function(feature_ids, X_z) {
    feature_ids <- intersect(feature_ids, rownames(X_z))
    if (length(feature_ids) == 0) {
        return(rep(NA_real_, ncol(X_z)))
    }
    if (length(feature_ids) == 1) {
        return(as.numeric(X_z[feature_ids, , drop = TRUE]))
    }
    apply(X_z[feature_ids, , drop = FALSE], 2, stats::median, na.rm = TRUE)
}
calculate_group_concordance <- function(analysis_group, feature_ids, X_log2) {
    feature_ids <- intersect(feature_ids, rownames(X_log2))
    if (length(feature_ids) == 1) {
        return(tibble::tibble(analysis_group = analysis_group, n_CLASSY_features = 1L, member_CLASSY_IDs = feature_ids, 
            median_pairwise_spearman = NA_real_, minimum_pairwise_spearman = NA_real_, maximum_pairwise_spearman = NA_real_, 
            low_concordance_flag = FALSE))
    }
    R <- stats::cor(t(X_log2[feature_ids, , drop = FALSE]), method = "spearman", use = "pairwise.complete.obs")
    pairwise_rho <- R[upper.tri(R)]
    tibble::tibble(analysis_group = analysis_group, n_CLASSY_features = length(feature_ids), member_CLASSY_IDs = paste(feature_ids, 
        collapse = "; "), median_pairwise_spearman = stats::median(pairwise_rho, na.rm = TRUE), minimum_pairwise_spearman = min(pairwise_rho, 
        na.rm = TRUE), maximum_pairwise_spearman = max(pairwise_rho, na.rm = TRUE), low_concordance_flag = stats::median(pairwise_rho, 
        na.rm = TRUE) < 0.5 | min(pairwise_rho, na.rm = TRUE) < 0)
}
add_limma_confidence_intervals <- function(fit, coefficient_name, result_table, id_column) {
    coefficient_index <- match(coefficient_name, colnames(fit$coefficients))
    moderated_se <- fit$stdev.unscaled[, coefficient_index] * sqrt(fit$s2.post)
    df_total <- fit$df.total
    critical_value <- stats::qt(0.975, df = df_total)
    ci_table <- tibble::tibble(analysis_id = rownames(fit$coefficients), moderated_SE = moderated_se, 
        df_total = df_total, conf_low = fit$coefficients[, coefficient_index] - critical_value * moderated_se, 
        conf_high = fit$coefficients[, coefficient_index] + critical_value * moderated_se)
    names(ci_table)[1] <- id_column
    result_table %>% dplyr::left_join(ci_table, by = id_column)
}
cecum_raw <- readxl::read_xlsx(ppi_input(cecum_xlsx), sheet = 1, .name_repair = "minimal")
cecum_meta_raw <- readxl::read_xlsx(ppi_input(cecum_xlsx), sheet = 2, .name_repair = "minimal")
names(cecum_raw) <- stringr::str_squish(names(cecum_raw))
names(cecum_meta_raw) <- stringr::str_squish(names(cecum_meta_raw))
cecum_feature_info_all <- tibble::tibble(feature_id = as.character(cecum_raw[[2]]), original_metabolite = stringr::str_squish(as.character(cecum_raw[[1]])), 
    canonical_metabolite = clean_metabolite_name(cecum_raw[[1]]), ppm = suppressWarnings(as.numeric(cecum_raw[[3]])))
if (anyDuplicated(cecum_feature_info_all$feature_id)) {
    stop("CLASSY IDs in caecum sheet 1 are not unique.")
}
sample_column_labels <- names(cecum_raw)[4:ncol(cecum_raw)]
names(cecum_meta_raw)[1] <- "sample_index"
if (length(sample_column_labels) != nrow(cecum_meta_raw)) {
    stop("Sheet-1 sample columns do not match sheet-2 metadata rows.\n", "Sheet-1 sample columns: ", 
        length(sample_column_labels), "\nSheet-2 metadata rows: ", nrow(cecum_meta_raw))
}
sample_info_all <- cecum_meta_raw %>% dplyr::mutate(sample_index = as.character(sample_index), sample_type = dplyr::case_when(stringr::str_to_upper(as.character(Inoculum)) == 
    "QC" ~ "QC", TRUE ~ "biological"), sample_id = dplyr::case_when(sample_type == "QC" ~ "QC", TRUE ~ 
    paste0("Mouse_", stringr::str_pad(sample_index, width = 2, pad = "0"))), original_matrix_label = sample_column_labels, 
    inoculum_donor = as.character(Inoculum), group = dplyr::case_when(stringr::str_detect(stringr::str_to_lower(inoculum_donor), 
        "^ctrl") ~ "ctrl", stringr::str_detect(stringr::str_to_lower(inoculum_donor), "^ppi") ~ "ppi", 
        sample_type == "QC" ~ "QC", TRUE ~ NA_character_), grp = factor(group, levels = c(group_levels, 
        "QC")), inoculum_donor = factor(inoculum_donor), mouse = as.character(mouse), weight_mg = suppressWarnings(as.numeric(`Weight [mg]`)), 
    isolator = factor(as.character(Isolator)), cage_position = factor(as.character(Käfigplatz)), age_weeks = suppressWarnings(as.numeric(`Alter (w)`)))
sample_info <- sample_info_all %>% dplyr::filter(sample_type == "biological") %>% dplyr::mutate(grp = factor(group, 
    levels = group_levels), inoculum_donor = droplevels(inoculum_donor), isolator = droplevels(isolator))
if (dplyr::n_distinct(sample_info$inoculum_donor) != 6) {
    warning("Expected six independent inoculum donors but found ", dplyr::n_distinct(sample_info$inoculum_donor), 
        ".")
}
feature_info <- cecum_feature_info_all %>% dplyr::filter(!is.na(canonical_metabolite), !canonical_metabolite %in% 
    c("baseline", "tsp")) %>% dplyr::mutate(analysis_group = canonical_metabolite, group_type = dplyr::case_when(analysis_group %in% 
    c("3-hydroxyphenylpropionate", "glutamine") ~ "combined_metabolite", analysis_group == "arginine+lysine" ~ 
    "ambiguous_combined_annotation", analysis_group %in% c("bile1", "bile2") ~ "unresolved_annotation", 
    TRUE ~ "metabolite"), interpretation_note = dplyr::case_when(analysis_group == "3-hydroxyphenylpropionate" ~ 
    paste("Multiple highly concordant annotated resonances", "combined into one standardized score."), 
    analysis_group == "glutamine" ~ paste("Three highly concordant annotated resonances", "combined into one standardized score."), 
    analysis_group == "arginine+lysine" ~ paste("Combined/ambiguous collaborator annotation retained", 
        "under the supplied compound label."), analysis_group %in% c("bile1", "bile2") ~ paste("Unresolved bile-associated annotation retained", 
        "under the collaborator-provided label."), TRUE ~ ""))
message("Annotated CLASSY features retained: ", nrow(feature_info))
message("Cleaned metabolite groups: ", dplyr::n_distinct(feature_info$analysis_group))
X_all <- as.matrix(cecum_raw[, 4:ncol(cecum_raw), drop = FALSE])
storage.mode(X_all) <- "numeric"
rownames(X_all) <- cecum_feature_info_all$feature_id
colnames(X_all) <- sample_info_all$sample_id
X_raw <- X_all[feature_info$feature_id, sample_info$sample_id, drop = FALSE]
negative_summary <- tibble::tibble(feature_id = rownames(X_raw), n_negative = rowSums(X_raw < 0, na.rm = TRUE), 
    proportion_negative = rowMeans(X_raw < 0, na.rm = TRUE)) %>% dplyr::left_join(feature_info, by = "feature_id")
X_nonnegative <- pmax(X_raw, 0)
feature_pseudocount <- apply(X_nonnegative, 1, make_feature_pseudocount)
valid_feature <- is.finite(feature_pseudocount)
if (any(!valid_feature)) {
    warning(sum(!valid_feature), " annotated features had no positive values and were removed.")
}
X_nonnegative <- X_nonnegative[valid_feature, , drop = FALSE]
feature_pseudocount <- feature_pseudocount[valid_feature]
feature_info <- feature_info %>% dplyr::filter(feature_id %in% rownames(X_nonnegative))
X_log2 <- X_nonnegative
for (i in seq_len(nrow(X_log2))) {
    X_log2[i, ] <- log2(X_nonnegative[i, ] + feature_pseudocount[i])
}
X_z <- t(scale(t(X_log2)))
group_members <- split(feature_info$feature_id, feature_info$analysis_group)
Y_group <- vapply(group_members, make_group_score, X_z = X_z, FUN.VALUE = numeric(ncol(X_z)))
Y_group <- t(Y_group)
colnames(Y_group) <- colnames(X_z)
group_variance <- apply(Y_group, 1, stats::var, na.rm = TRUE)
Y_group <- Y_group[is.finite(group_variance) & group_variance > 0, , drop = FALSE]
group_concordance <- dplyr::bind_rows(lapply(names(group_members), function(analysis_group) {
    calculate_group_concordance(analysis_group = analysis_group, feature_ids = group_members[[analysis_group]], 
        X_log2 = X_log2)
}))
group_summary <- feature_info %>% dplyr::filter(analysis_group %in% rownames(Y_group)) %>% dplyr::group_by(analysis_group, 
    group_type, interpretation_note) %>% dplyr::summarise(n_CLASSY_features = dplyr::n_distinct(feature_id), 
    member_CLASSY_IDs = paste(feature_id, collapse = "; "), original_names = paste(sort(unique(original_metabolite)), 
        collapse = "; "), ppm_values = paste(sort(unique(round(ppm, 5))), collapse = "; "), ppm_min = min(ppm, 
        na.rm = TRUE), ppm_max = max(ppm, na.rm = TRUE), .groups = "drop") %>% dplyr::left_join(group_concordance, 
    by = c("analysis_group", "n_CLASSY_features", "member_CLASSY_IDs"))
donor_levels <- levels(sample_info$inoculum_donor)
Y_donor <- vapply(donor_levels, function(donor_id) {
    rowMeans(Y_group[, sample_info$inoculum_donor == donor_id, drop = FALSE], na.rm = TRUE)
}, FUN.VALUE = numeric(nrow(Y_group)))
rownames(Y_donor) <- rownames(Y_group)
colnames(Y_donor) <- donor_levels
donor_info <- sample_info %>% dplyr::distinct(inoculum_donor, grp, group) %>% dplyr::arrange(match(inoculum_donor, 
    donor_levels))
design_donor <- stats::model.matrix(~grp, data = donor_info)
colnames(design_donor) <- make.names(colnames(design_donor))
group_coefficient <- "grpppi"
fit_donor <- limma::lmFit(Y_donor, design_donor)
fit_donor <- limma::eBayes(fit_donor, robust = TRUE)
donor_results <- limma::topTable(fit_donor, coef = group_coefficient, number = Inf, sort.by = "P") %>% 
    tibble::rownames_to_column("analysis_group") %>% dplyr::rename(score_difference = logFC, average_score = AveExpr, 
    t_statistic = t, p_value = P.Value, FDR = adj.P.Val, log_odds = B)
donor_results <- add_limma_confidence_intervals(fit = fit_donor, coefficient_name = group_coefficient, 
    result_table = donor_results, id_column = "analysis_group") %>% dplyr::left_join(group_summary, by = "analysis_group") %>% 
    dplyr::mutate(direction = dplyr::case_when(score_difference > 0 ~ "higher in PPI", score_difference < 
        0 ~ "lower in PPI", TRUE ~ "no difference"), analysis = "Primary: donor-averaged limma", independent_unit = "microbiome inoculum donor", 
        effect_scale = "difference in donor-averaged standardized metabolite score") %>% dplyr::arrange(p_value)
design_mouse <- stats::model.matrix(~grp, data = sample_info)
colnames(design_mouse) <- make.names(colnames(design_mouse))
dupcor_mouse <- limma::duplicateCorrelation(Y_group, design = design_mouse, block = sample_info$inoculum_donor)
mouse_consensus_correlation <- dupcor_mouse$consensus.correlation
if (!is.finite(mouse_consensus_correlation)) {
    warning("Mouse-level duplicateCorrelation was not finite; ", "using correlation = 0.")
    mouse_consensus_correlation <- 0
}
fit_mouse <- limma::lmFit(Y_group, design = design_mouse, block = sample_info$inoculum_donor, correlation = mouse_consensus_correlation)
fit_mouse <- limma::eBayes(fit_mouse, robust = TRUE)
mouse_dupcor_results <- limma::topTable(fit_mouse, coef = group_coefficient, number = Inf, sort.by = "P") %>% 
    tibble::rownames_to_column("analysis_group") %>% dplyr::rename(score_difference = logFC, average_score = AveExpr, 
    t_statistic = t, p_value = P.Value, FDR = adj.P.Val, log_odds = B)
mouse_dupcor_results <- add_limma_confidence_intervals(fit = fit_mouse, coefficient_name = group_coefficient, 
    result_table = mouse_dupcor_results, id_column = "analysis_group") %>% dplyr::left_join(group_summary, 
    by = "analysis_group") %>% dplyr::mutate(direction = dplyr::case_when(score_difference > 0 ~ "higher in PPI", 
    score_difference < 0 ~ "lower in PPI", TRUE ~ "no difference"), analysis = "Sensitivity: mouse-level limma + duplicateCorrelation", 
    independent_unit = "microbiome inoculum donor", effect_scale = "difference in standardized metabolite score") %>% 
    dplyr::arrange(p_value)
adjusted_mouse_results <- tibble::tibble()
adjusted_consensus_correlation <- NA_real_
if (run_adjusted_mouse_sensitivity) {
    design_mouse_adjusted <- stats::model.matrix(~grp + isolator + scale(age_weeks), data = sample_info)
    colnames(design_mouse_adjusted) <- make.names(colnames(design_mouse_adjusted))
    dupcor_mouse_adjusted <- limma::duplicateCorrelation(Y_group, design = design_mouse_adjusted, block = sample_info$inoculum_donor)
    adjusted_consensus_correlation <- dupcor_mouse_adjusted$consensus.correlation
    if (!is.finite(adjusted_consensus_correlation)) {
        adjusted_consensus_correlation <- 0
    }
    fit_mouse_adjusted <- limma::lmFit(Y_group, design = design_mouse_adjusted, block = sample_info$inoculum_donor, 
        correlation = adjusted_consensus_correlation)
    fit_mouse_adjusted <- limma::eBayes(fit_mouse_adjusted, robust = TRUE)
    adjusted_mouse_results <- limma::topTable(fit_mouse_adjusted, coef = group_coefficient, number = Inf, 
        sort.by = "P") %>% tibble::rownames_to_column("analysis_group") %>% dplyr::rename(score_difference = logFC, 
        average_score = AveExpr, t_statistic = t, p_value = P.Value, FDR = adj.P.Val, log_odds = B)
    adjusted_mouse_results <- add_limma_confidence_intervals(fit = fit_mouse_adjusted, coefficient_name = group_coefficient, 
        result_table = adjusted_mouse_results, id_column = "analysis_group") %>% dplyr::left_join(group_summary, 
        by = "analysis_group") %>% dplyr::mutate(direction = dplyr::case_when(score_difference > 0 ~ 
        "higher in PPI", score_difference < 0 ~ "lower in PPI", TRUE ~ "no difference"), analysis = paste("Sensitivity: mouse-level duplicateCorrelation", 
        "adjusted for recipient age and isolator"), independent_unit = "microbiome inoculum donor", effect_scale = "adjusted difference in standardized metabolite score") %>% 
        dplyr::arrange(p_value)
}
dupcor_feature <- limma::duplicateCorrelation(X_log2, design = design_mouse, block = sample_info$inoculum_donor)
feature_consensus_correlation <- dupcor_feature$consensus.correlation
if (!is.finite(feature_consensus_correlation)) {
    feature_consensus_correlation <- 0
}
fit_feature <- limma::lmFit(X_log2, design = design_mouse, block = sample_info$inoculum_donor, correlation = feature_consensus_correlation)
fit_feature <- limma::eBayes(fit_feature, robust = TRUE)
feature_results <- limma::topTable(fit_feature, coef = group_coefficient, number = Inf, sort.by = "P") %>% 
    tibble::rownames_to_column("feature_id") %>% dplyr::rename(log2_difference = logFC, p_value = P.Value, 
    FDR = adj.P.Val) %>% dplyr::left_join(feature_info, by = "feature_id") %>% dplyr::mutate(direction = dplyr::case_when(log2_difference > 
    0 ~ "higher in PPI", log2_difference < 0 ~ "lower in PPI", TRUE ~ "no difference"), analysis = "Annotated CLASSY feature-level duplicateCorrelation") %>% 
    dplyr::arrange(p_value)
group_long <- as.data.frame(t(Y_group)) %>% tibble::rownames_to_column("sample_id") %>% dplyr::left_join(sample_info %>% 
    dplyr::select(sample_id, inoculum_donor, grp, group, isolator, age_weeks, weight_mg, mouse), by = "sample_id") %>% 
    tidyr::pivot_longer(cols = -c(sample_id, inoculum_donor, grp, group, isolator, age_weeks, weight_mg, 
        mouse), names_to = "analysis_group", values_to = "score")
donor_average_long <- group_long %>% dplyr::group_by(analysis_group, inoculum_donor, grp, group) %>% 
    dplyr::summarise(donor_mean_score = mean(score, na.rm = TRUE), n_recipient_mice = dplyr::n(), .groups = "drop")
plot_cecum_metabolite_mice <- function(metabolite, data = group_long) {
    plot_df <- data %>% dplyr::filter(analysis_group == metabolite)
    if (nrow(plot_df) == 0) {
        stop("Metabolite not found: ", metabolite)
    }
    ggplot2::ggplot(plot_df, ggplot2::aes(x = grp, y = score, colour = inoculum_donor)) + ggplot2::geom_hline(yintercept = 0, 
        linetype = "dashed", linewidth = 0.4) + ggplot2::geom_point(size = 3, position = ggplot2::position_jitter(width = 0.07, 
        height = 0)) + ggplot2::stat_summary(ggplot2::aes(group = grp), fun = mean, geom = "point", shape = 18, 
        size = 5, colour = "black") + ggplot2::scale_x_discrete(labels = c(ctrl = "Control", ppi = "PPI")) + 
        ggplot2::labs(title = stringr::str_replace_all(metabolite, "_", " "), subtitle = "Recipient mice coloured by microbiome inoculum donor", 
            x = NULL, y = "Standardized metabolite score", colour = "Inoculum donor") + ggplot2::theme_bw() + 
        ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"), legend.position = "top")
}
plot_cecum_metabolite_donor_average <- function(metabolite, data = donor_average_long, results = donor_results) {
    plot_df <- data %>% dplyr::filter(analysis_group == metabolite)
    if (nrow(plot_df) == 0) {
        stop("Metabolite not found: ", metabolite)
    }
    result_row <- results %>% dplyr::filter(analysis_group == metabolite)
    subtitle_text <- NULL
    if (nrow(result_row) == 1) {
        subtitle_text <- paste0("Donor-averaged limma: p = ", format.pval(result_row$p_value, digits = 3, 
            eps = 0.001), "; FDR = ", format.pval(result_row$FDR, digits = 3, eps = 0.001))
    }
    ggplot2::ggplot(plot_df, ggplot2::aes(x = grp, y = donor_mean_score, colour = grp, label = inoculum_donor)) + 
        ggplot2::geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.4) + ggplot2::geom_point(size = 4, 
        position = ggplot2::position_jitter(width = 0.03, height = 0)) + ggplot2::geom_text(colour = "black", 
        nudge_x = 0.05, hjust = 0, size = 3) + ggplot2::stat_summary(fun = mean, geom = "point", shape = 18, 
        size = 5, colour = "black") + ggplot2::scale_x_discrete(labels = c(ctrl = "Control", ppi = "PPI")) + 
        ggplot2::guides(colour = "none") + ggplot2::labs(title = stringr::str_replace_all(metabolite, 
        "_", " "), subtitle = subtitle_text, x = NULL, y = "Mean standardized score across recipient mice") + 
        ggplot2::theme_bw() + ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"))
}
analysis_summary <- tibble::tibble(item = c("Biological recipient mice", "Independent inoculum donors", 
    "Control inoculum donors", "PPI inoculum donors", "Annotated CLASSY features retained", "Cleaned metabolite groups", 
    "Primary analysis", "Mouse-level sensitivity", "Mouse-level consensus correlation", "Adjusted consensus correlation", 
    "Feature-level consensus correlation", "Negative-value handling", "Duplicate-feature combination"), 
    value = c(nrow(sample_info), dplyr::n_distinct(sample_info$inoculum_donor), dplyr::n_distinct(sample_info$inoculum_donor[sample_info$group == 
        "ctrl"]), dplyr::n_distinct(sample_info$inoculum_donor[sample_info$group == "ppi"]), nrow(feature_info), 
        nrow(Y_group), "donor-averaged limma on six inoculum-donor means", "limma + duplicateCorrelation on 16 mice, block = inoculum donor", 
        mouse_consensus_correlation, adjusted_consensus_correlation, feature_consensus_correlation, "negative peak areas truncated to zero before log2 transformation", 
        "median of z-scored CLASSY features sharing a cleaned metabolite name"))
settings <- tibble::tibble(setting = c("cecum_xlsx", "matrix_sheet", "metadata_sheet", "annotated_only", 
    "excluded_annotations", "manual_name_cleaning", "primary_multiple_testing", "primary_independent_unit", 
    "endpoint_weight_adjustment"), value = c(cecum_xlsx, "Tabelle1", "Tabelle2", "TRUE", "baseline; TSP", 
    paste("3-(3-OH-phenyl)propionate -> 3-hydroxyphenylpropionate;", "5-Aminopentanoic acid -> 5-aminopentanoate;", 
        "phosphocoline -> phosphocholine"), "BH FDR across cleaned annotated metabolite groups", "microbiome inoculum donor", 
    "not included in primary model"))
output_xlsx <- file.path(out_dir, "NMR_cecum_annotated_limma_results.xlsx")
writexl::write_xlsx(list(Donor_averaged_primary = donor_results, Mouse_dupCor_sensitivity = mouse_dupcor_results, 
    Mouse_adjusted_sensitivity = adjusted_mouse_results, CLASSY_feature_sensitivity = feature_results, 
    Cleaned_feature_mapping = feature_info, Group_summary = group_summary, Group_concordance = group_concordance, 
    Negative_value_summary = negative_summary, Donor_average_values = donor_average_long, Mouse_level_values = group_long, 
    Sample_metadata = sample_info, Analysis_summary = analysis_summary, Settings = settings), ppi_output(output_xlsx))
readr::write_csv(donor_results, ppi_output(file.path(out_dir, "NMR_cecum_donor_averaged_primary_results.csv"), 
    name = "nmr_caecum_donor_averaged_primary_results.csv"))
readr::write_csv(mouse_dupcor_results, ppi_output(file.path(out_dir, "NMR_cecum_mouse_duplicateCorrelation_results.csv"), 
    name = "nmr_caecum_mouse_duplicate_correlation_results.csv"))
readr::write_csv(adjusted_mouse_results, ppi_output(file.path(out_dir, "NMR_cecum_mouse_adjusted_duplicateCorrelation_results.csv"), 
    name = "nmr_caecum_mouse_adjusted_duplicate_correlation_results.csv"))
readr::write_csv(feature_info, ppi_output(file.path(out_dir, "NMR_cecum_cleaned_feature_mapping.csv"), 
    name = "nmr_caecum_cleaned_feature_mapping.csv"))
message("\nCaecum analysis complete.\n", "Main workbook:\n", output_xlsx, "\n\nPrimary independent unit: inoculum donor\n", 
    "Number of independent donors: ", ncol(Y_donor), "\nMouse-level consensus correlation: ", round(mouse_consensus_correlation, 
        4), "\n")
print(donor_results %>% dplyr::select(analysis_group, n_CLASSY_features, score_difference, conf_low, 
    conf_high, p_value, FDR, direction, low_concordance_flag) %>% head(20))
