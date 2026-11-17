# SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0
# Ex vivo HILIC limma: combined-mode primary and mode-separated sensitivity.

source("scripts/helpers/setup.R")
invisible(run_analysis(
  section = "04_mouse_metabolomics",
  directory = "metabolomics/20250714_MassSpecData",
  script = "04_09_analyse_ex_vivo_hilic.R",
  restricted = FALSE,
  code = {
required_packages <- c("limma", "statmod", "dplyr", "tidyr", "purrr", "stringr", "tibble", "readr", "writexl")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]
if (length(missing_packages) > 0) {
    stop("Install these packages first: ", paste(missing_packages, collapse = ", "), "\nFor limma, use BiocManager::install('limma').")
}
input_rds <- project_path("metabolomics/20250714_MassSpecData/Ex_vivo_cleaned_one_per_metabolite/HILIC_exvivo_cleaned_analysis_objects.rds")
out_dir <- project_path("metabolomics/20250714_MassSpecData/Ex_vivo_cleaned_one_per_metabolite/Limma_analysis")
dir.create(project_output(out_dir, directory = TRUE), recursive = TRUE, showWarnings = FALSE)
analysis_set <- "primary_rank1_rank2"
keep_days <- c("d1", "d2", "d3", "d4")
group_levels <- c("ctrl", "ppi")
min_positive_per_group <- 4L
min_donors_detected_per_group <- 2L
strict_min_positive_per_group <- 6L
robust_ebayes <- TRUE
safe_min_positive <- function(x) {
    x <- x[is.finite(x) & x > 0]
    if (length(x) == 0)
        NA_real_
    else min(x)
}
extract_limma_result <- function(fit, coefficient, analysis_label) {
    if (!coefficient %in% colnames(fit$coefficients)) {
        stop("Coefficient ", coefficient, " not present. Available: ", paste(colnames(fit$coefficients),
            collapse = ", "))
    }
    coef_index <- match(coefficient, colnames(fit$coefficients))
    out <- limma::topTable(fit, coef = coefficient, number = Inf, sort.by = "P") %>% tibble::rownames_to_column("metabolite") %>%
        dplyr::rename(log2FC = logFC, average_log2_abundance = AveExpr, t_statistic = t, p_value = P.Value,
            FDR = adj.P.Val, log_odds = B)
    moderated_se <- fit$stdev.unscaled[, coef_index] * sqrt(fit$s2.post)
    df_total <- fit$df.total
    critical <- stats::qt(0.975, df = df_total)
    ci <- tibble::tibble(metabolite = rownames(fit$coefficients), moderated_SE = moderated_se, df_total = df_total,
        CI_low_log2 = fit$coefficients[, coef_index] - critical * moderated_se, CI_high_log2 = fit$coefficients[,
            coef_index] + critical * moderated_se)
    out %>% dplyr::left_join(ci, by = "metabolite") %>% dplyr::mutate(fold_change = 2^log2FC, CI_low_fold_change = 2^CI_low_log2,
        CI_high_fold_change = 2^CI_high_log2, analysis = analysis_label)
}
fit_duplicate_correlation <- function(Y, design, block, coefficient, label) {
    corfit <- limma::duplicateCorrelation(object = Y, design = design, block = block)
    consensus <- corfit$consensus.correlation
    if (!is.finite(consensus)) {
        warning(label, ": duplicateCorrelation was not finite; using 0.")
        consensus <- 0
    }
    fit <- limma::lmFit(object = Y, design = design, block = block, correlation = consensus)
    fit <- limma::eBayes(fit, robust = robust_ebayes)
    list(fit = fit, results = extract_limma_result(fit, coefficient, label), consensus_correlation = consensus)
}
fit_simple_limma <- function(Y, design, coefficient, label) {
    fit <- limma::lmFit(Y, design)
    fit <- limma::eBayes(fit, robust = robust_ebayes)
    list(fit = fit, results = extract_limma_result(fit, coefficient, label))
}
fit_primary_by_mode <- function(Y, design, block, coefficient, metabolite_mode_map, label_prefix) {
    required_map_columns <- c("metabolite", "mode")
    missing_map_columns <- setdiff(required_map_columns, names(metabolite_mode_map))
    if (length(missing_map_columns) > 0L) {
        stop("metabolite_mode_map is missing: ", paste(missing_map_columns, collapse = ", "))
    }
    fitted_modes <- lapply(c("Positive", "Negative"), function(mode_value) {
        metabolites_mode <- metabolite_mode_map %>% dplyr::filter(.data$mode == mode_value) %>% dplyr::pull(metabolite) %>%
            intersect(rownames(Y))
        if (length(metabolites_mode) < 2L) {
            warning("Only ", length(metabolites_mode), " metabolites available in ", mode_value, " mode; this mode was not fitted.")
            return(NULL)
        }
        fit_mode <- fit_duplicate_correlation(Y = Y[metabolites_mode, , drop = FALSE], design = design,
            block = block, coefficient = coefficient, label = paste0(label_prefix, " [", mode_value,
                " mode]"))
        results_mode <- fit_mode$results %>% dplyr::rename(FDR_within_mode = FDR) %>% dplyr::mutate(mode_fit = mode_value,
            n_metabolites_in_mode_fit = length(metabolites_mode), consensus_donor_correlation_mode = fit_mode$consensus_correlation)
        list(mode = mode_value, fit = fit_mode$fit, results = results_mode, consensus_correlation = fit_mode$consensus_correlation,
            n_metabolites = length(metabolites_mode))
    })
    fitted_modes <- Filter(Negate(is.null), fitted_modes)
    if (length(fitted_modes) == 0L) {
        stop("Neither ionization mode could be fitted.")
    }
    combined_results <- dplyr::bind_rows(lapply(fitted_modes, function(x) x$results)) %>% dplyr::mutate(FDR_global = stats::p.adjust(.data$p_value,
        method = "BH"), FDR = .data$FDR_global) %>% dplyr::arrange(.data$FDR_global, .data$p_value)
    list(results = combined_results, positive_results = combined_results %>% dplyr::filter(.data$mode_fit ==
        "Positive"), negative_results = combined_results %>% dplyr::filter(.data$mode_fit == "Negative"),
        fits = fitted_modes)
}
if (!project_file_exists(input_rds)) {
    stop("Input RDS not found:\n", input_rds)
}
obj <- readRDS(project_input(input_rds))
if (analysis_set == "primary_rank1_rank2") {
    features_all <- as.matrix(obj$primary_matrix_all_samples)
    selected_map <- as.data.frame(obj$primary_selected_map, check.names = FALSE)
} else if (analysis_set == "exploratory_all_ranks") {
    features_all <- as.matrix(obj$exploratory_matrix_all_samples)
    selected_map <- as.data.frame(obj$exploratory_selected_map, check.names = FALSE)
} else {
    stop("analysis_set must be 'primary_rank1_rank2' or 'exploratory_all_ranks'.")
}
metadata_all <- as.data.frame(obj$metadata_all, check.names = FALSE)
storage.mode(features_all) <- "numeric"
if (!identical(rownames(features_all), metadata_all$Sample_Name)) {
    metadata_order <- match(rownames(features_all), metadata_all$Sample_Name)
    if (anyNA(metadata_order)) {
        stop("Feature-matrix sample names cannot be matched to metadata.")
    }
    metadata_all <- metadata_all[metadata_order, , drop = FALSE]
}
metadata_all$Separation <- factor(as.character(metadata_all$Separation), levels = c(group_levels, "medium"))
metadata_all$Donor <- factor(as.character(metadata_all$Donor))
metadata_all$Day <- as.character(metadata_all$Day)
negative_count <- sum(features_all < 0, na.rm = TRUE)
if (negative_count > 0) {
    warning(negative_count, " negative intensities were truncated to zero before transformation.")
    features_all[features_all < 0] <- 0
}
analysis_index <- which(metadata_all$Separation %in% group_levels & metadata_all$Day %in% keep_days)
biological_all_index <- which(metadata_all$Separation %in% group_levels & metadata_all$Day %in% c("d0",
    keep_days))
metadata <- metadata_all[analysis_index, , drop = FALSE]
features_raw <- features_all[analysis_index, , drop = FALSE]
metadata$Separation <- factor(as.character(metadata$Separation), levels = group_levels)
metadata$DayF <- factor(metadata$Day, levels = keep_days)
metadata$Donor <- factor(as.character(metadata$Donor))
if (nrow(metadata) != 24L) {
    warning("Expected 24 d1-d4 biological samples, but found ", nrow(metadata), ".")
}
prevalence <- purrr::map_dfr(colnames(features_raw), function(metabolite) {
    values <- features_raw[, metabolite]
    one_group <- function(group_value) {
        idx <- which(metadata$Separation == group_value)
        donors <- unique(metadata$Donor[idx])
        donor_detected <- vapply(donors, function(donor_value) {
            donor_idx <- idx[metadata$Donor[idx] == donor_value]
            any(values[donor_idx] > 0, na.rm = TRUE)
        }, FUN.VALUE = logical(1))
        c(n_positive = sum(values[idx] > 0, na.rm = TRUE), n_donors = sum(donor_detected))
    }
    ctrl <- one_group("ctrl")
    ppi <- one_group("ppi")
    tibble::tibble(metabolite = metabolite, n_positive_ctrl = unname(ctrl["n_positive"]), n_positive_ppi = unname(ppi["n_positive"]),
        n_donors_detected_ctrl = unname(ctrl["n_donors"]), n_donors_detected_ppi = unname(ppi["n_donors"]),
        passes_primary_prevalence = n_positive_ctrl >= min_positive_per_group & n_positive_ppi >= min_positive_per_group &
            n_donors_detected_ctrl >= min_donors_detected_per_group & n_donors_detected_ppi >= min_donors_detected_per_group,
        passes_strict_prevalence = n_positive_ctrl >= strict_min_positive_per_group & n_positive_ppi >=
            strict_min_positive_per_group & n_donors_detected_ctrl >= min_donors_detected_per_group &
            n_donors_detected_ppi >= min_donors_detected_per_group)
})
kept_metabolites <- prevalence$metabolite[prevalence$passes_primary_prevalence]
if (length(kept_metabolites) < 2L) {
    stop("Fewer than two metabolites passed prevalence filtering.")
}
features_all_kept <- features_all[, kept_metabolites, drop = FALSE]
features_raw <- features_raw[, kept_metabolites, drop = FALSE]
pseudocounts <- vapply(kept_metabolites, function(metabolite) {
    minimum_positive <- safe_min_positive(features_all[biological_all_index, metabolite])
    if (!is.finite(minimum_positive))
        NA_real_
    else 0.5 * minimum_positive
}, FUN.VALUE = numeric(1))
if (any(!is.finite(pseudocounts))) {
    bad <- names(pseudocounts)[!is.finite(pseudocounts)]
    stop("No positive biological value available for: ", paste(bad, collapse = ", "))
}
features_log2_all <- sweep(features_all_kept, 2, pseudocounts, FUN = "+")
features_log2_all <- log2(features_log2_all)
features_log2 <- features_log2_all[analysis_index, , drop = FALSE]
Y <- t(features_log2)
pseudocount_table <- tibble::tibble(metabolite = names(pseudocounts), pseudocount = as.numeric(pseudocounts))
primary_design <- stats::model.matrix(~Separation + DayF, data = metadata)
colnames(primary_design) <- make.names(colnames(primary_design))
primary_coefficient <- "Separationppi"
primary_fit <- fit_duplicate_correlation(Y = Y, design = primary_design, block = metadata$Donor, coefficient = primary_coefficient,
    label = "Unaveraged d1-d4: group + categorical day; block donor")
primary_results <- primary_fit$results %>% dplyr::mutate(consensus_donor_correlation = primary_fit$consensus_correlation)
metabolite_mode_map <- selected_map %>% dplyr::transmute(metabolite = .data$canonical_name, mode = stringr::str_to_title(as.character(.data$mode))) %>%
    dplyr::distinct()
mode_map_duplicates <- metabolite_mode_map %>% dplyr::count(.data$metabolite, name = "n_mode_rows") %>%
    dplyr::filter(.data$n_mode_rows > 1L)
if (nrow(mode_map_duplicates) > 0L) {
    stop("The selected one-per-metabolite map contains more than one selected mode ", "for: ", paste(mode_map_duplicates$metabolite,
        collapse = ", "))
}
missing_modes <- setdiff(rownames(Y), metabolite_mode_map$metabolite)
if (length(missing_modes) > 0L) {
    stop("Selected ionization mode is missing for: ", paste(missing_modes, collapse = ", "))
}
mode_separated_primary <- fit_primary_by_mode(Y = Y, design = primary_design, block = metadata$Donor,
    coefficient = primary_coefficient, metabolite_mode_map = metabolite_mode_map, label_prefix = "Mode-separated unaveraged d1-d4: group + categorical day; block donor")
mode_separated_primary_results <- mode_separated_primary$results
mode_separated_positive_results <- mode_separated_primary$positive_results
mode_separated_negative_results <- mode_separated_primary$negative_results
features_log2_df <- tibble::as_tibble(features_log2, .name_repair = "minimal") %>% dplyr::mutate(Donor = metadata$Donor,
    Separation = metadata$Separation)
donor_metadata <- metadata %>% dplyr::select(Donor, Separation) %>% dplyr::distinct() %>% dplyr::arrange(Separation,
    Donor)
donor_profiles <- features_log2_df %>% dplyr::group_by(Donor, Separation) %>% dplyr::summarise(dplyr::across(dplyr::all_of(kept_metabolites),
    mean, na.rm = TRUE), .groups = "drop")
donor_profiles$Separation <- factor(as.character(donor_profiles$Separation), levels = group_levels)
donor_design <- stats::model.matrix(~Separation, data = donor_profiles)
colnames(donor_design) <- make.names(colnames(donor_design))
Y_donor <- t(as.matrix(donor_profiles[, kept_metabolites, drop = FALSE]))
colnames(Y_donor) <- as.character(donor_profiles$Donor)
donor_fit <- fit_simple_limma(Y = Y_donor, design = donor_design, coefficient = "Separationppi", label = "Sensitivity: donor-averaged d1-d4 profiles")
donor_results <- donor_fit$results
baseline_rows <- which(metadata_all$Separation %in% group_levels & metadata_all$Day == "d0")
baseline_metadata <- metadata_all[baseline_rows, , drop = FALSE]
baseline_matrix <- features_log2_all[baseline_rows, , drop = FALSE]
if (anyDuplicated(baseline_metadata$Donor)) {
    stop("More than one day-0 sample exists for at least one donor.")
}
baseline_match <- match(as.character(metadata$Donor), as.character(baseline_metadata$Donor))
if (anyNA(baseline_match)) {
    stop("At least one d1-d4 sample has no matching day-0 donor sample.")
}
baseline_adjusted <- features_log2 - baseline_matrix[baseline_match, , drop = FALSE]
Y_baseline <- t(baseline_adjusted)
baseline_fit <- fit_duplicate_correlation(Y = Y_baseline, design = primary_design, block = metadata$Donor,
    coefficient = primary_coefficient, label = "Sensitivity: change from donor-specific day 0")
baseline_results <- baseline_fit$results %>% dplyr::mutate(consensus_donor_correlation = baseline_fit$consensus_correlation)
metadata$Cell <- factor(paste(metadata$Separation, metadata$Day, sep = "_"), levels = as.vector(outer(group_levels,
    keep_days, paste, sep = "_")))
cell_design <- stats::model.matrix(~0 + Cell, data = metadata)
colnames(cell_design) <- sub("^Cell", "", colnames(cell_design))
cell_corfit <- limma::duplicateCorrelation(Y, design = cell_design, block = metadata$Donor)
cell_consensus <- cell_corfit$consensus.correlation
if (!is.finite(cell_consensus)) cell_consensus <- 0
cell_fit <- limma::lmFit(Y, design = cell_design, block = metadata$Donor, correlation = cell_consensus)
contrast_matrix <- limma::makeContrasts(PPI_effect_d1 = ppi_d1 - ctrl_d1, PPI_effect_d2 = ppi_d2 - ctrl_d2,
    PPI_effect_d3 = ppi_d3 - ctrl_d3, PPI_effect_d4 = ppi_d4 - ctrl_d4, levels = cell_design)
cell_fit2 <- limma::contrasts.fit(cell_fit, contrast_matrix)
cell_fit2 <- limma::eBayes(cell_fit2, robust = robust_ebayes)
day_specific_results <- purrr::map_dfr(colnames(contrast_matrix), function(coef_name) {
    extract_limma_result(cell_fit2, coefficient = coef_name, analysis_label = paste("Day-specific", coef_name)) %>%
        dplyr::mutate(day = stringr::str_remove(coef_name, "^PPI_effect_"))
})
interaction_design <- stats::model.matrix(~Separation * DayF, data = metadata)
colnames(interaction_design) <- make.names(colnames(interaction_design))
interaction_corfit <- limma::duplicateCorrelation(Y, design = interaction_design, block = metadata$Donor)
interaction_consensus <- interaction_corfit$consensus.correlation
if (!is.finite(interaction_consensus)) interaction_consensus <- 0
interaction_fit <- limma::lmFit(Y, design = interaction_design, block = metadata$Donor, correlation = interaction_consensus)
interaction_fit <- limma::eBayes(interaction_fit, robust = robust_ebayes)
interaction_columns <- grep("^Separationppi.DayF", colnames(interaction_fit$coefficients), value = TRUE)
if (length(interaction_columns) > 0) {
    interaction_results <- limma::topTable(interaction_fit, coef = interaction_columns, number = Inf,
        sort.by = "F") %>% tibble::rownames_to_column("metabolite") %>% dplyr::rename(average_log2_abundance = AveExpr,
        F_statistic = F, p_value = P.Value, FDR = adj.P.Val) %>% dplyr::mutate(analysis = "Group-by-day interaction",
        consensus_donor_correlation = interaction_consensus)
} else {
    interaction_results <- tibble::tibble(metabolite = character(), analysis = character())
}
selected_annotation_export <- selected_map %>% dplyr::select(canonical_name, feature_key, mode, rank,
    selected_annotation_raw, annotation_flag, duplicate_concordance_class, median_pairwise_spearman,
    passes_representative_prevalence) %>% dplyr::rename(metabolite = canonical_name)
add_annotation <- function(x) {
    x %>% dplyr::left_join(selected_annotation_export, by = "metabolite") %>% dplyr::left_join(prevalence %>%
        dplyr::select(metabolite, passes_primary_prevalence, passes_strict_prevalence), by = "metabolite")
}
primary_results <- add_annotation(primary_results)
mode_separated_primary_results <- add_annotation(mode_separated_primary_results)
mode_separated_positive_results <- add_annotation(mode_separated_positive_results)
mode_separated_negative_results <- add_annotation(mode_separated_negative_results)
donor_results <- add_annotation(donor_results)
baseline_results <- add_annotation(baseline_results)
day_specific_results <- add_annotation(day_specific_results)
interaction_results <- add_annotation(interaction_results)
mode_separated_vs_combined <- primary_results %>% dplyr::select(metabolite, selected_mode = mode, combined_mode_log2FC = log2FC,
    combined_mode_p = p_value, combined_mode_FDR = FDR, combined_mode_consensus_correlation = consensus_donor_correlation) %>%
    dplyr::left_join(mode_separated_primary_results %>% dplyr::select(metabolite, mode_fit, mode_separated_log2FC = log2FC,
        mode_separated_p = p_value, FDR_within_mode, FDR_global, mode_separated_consensus_correlation = consensus_donor_correlation_mode),
        by = "metabolite") %>% dplyr::mutate(log2FC_difference_mode_separated_minus_combined = .data$mode_separated_log2FC -
    .data$combined_mode_log2FC, same_effect_direction = sign(.data$mode_separated_log2FC) == sign(.data$combined_mode_log2FC))
method_comparison <- primary_results %>% dplyr::select(metabolite, primary_log2FC = log2FC, primary_p = p_value,
    primary_FDR = FDR) %>% dplyr::left_join(donor_results %>% dplyr::select(metabolite, donor_average_log2FC = log2FC,
    donor_average_p = p_value, donor_average_FDR = FDR), by = "metabolite") %>% dplyr::left_join(baseline_results %>%
    dplyr::select(metabolite, baseline_adjusted_log2FC = log2FC, baseline_adjusted_p = p_value, baseline_adjusted_FDR = FDR),
    by = "metabolite") %>% dplyr::mutate(same_direction_primary_donor = sign(primary_log2FC) == sign(donor_average_log2FC),
    same_direction_primary_baseline = sign(primary_log2FC) == sign(baseline_adjusted_log2FC))
analysis_summary <- tibble::tribble(~item, ~value, "Analysis set", analysis_set, "Biological samples analysed",
    as.character(nrow(metadata)), "Independent donors", as.character(dplyr::n_distinct(metadata$Donor)),
    "Days analysed", paste(keep_days, collapse = ", "), "Primary model", "log2 abundance ~ group + categorical day; block = donor",
    "Primary combined-mode consensus donor correlation", as.character(primary_fit$consensus_correlation),
    "Mode-separated positive consensus donor correlation", as.character(unique(mode_separated_positive_results$consensus_donor_correlation_mode)),
    "Mode-separated negative consensus donor correlation", as.character(unique(mode_separated_negative_results$consensus_donor_correlation_mode)),
    "Mode-separated FDR procedure", "Positive and negative modes fitted separately; BH recalculated globally after recombination",
    "Baseline-adjusted consensus donor correlation", as.character(baseline_fit$consensus_correlation),
    "Metabolites before prevalence filtering", as.character(ncol(features_all)), "Metabolites after prevalence filtering",
    as.character(length(kept_metabolites)), "Minimum positive samples per group", as.character(min_positive_per_group),
    "Minimum detected donors per group", as.character(min_donors_detected_per_group), "Medium sample",
    "Not modelled and not subtracted", "Positive log2FC", "Higher in PPI-associated cultures")
metadata_export <- metadata %>% dplyr::select(Sample_Name, Donor, Day, Separation, Sample_Name_Positive,
    Sample_Name_Negative)
writexl::write_xlsx(list(Unaveraged_d1_d4_dupCor = primary_results, Primary_dupCor = primary_results,
    Mode_sep_primary_globalFDR = mode_separated_primary_results, Mode_sep_positive = mode_separated_positive_results,
    Mode_sep_negative = mode_separated_negative_results, Mode_sep_vs_combined = mode_separated_vs_combined,
    Donor_averaged = donor_results, Baseline_adjusted = baseline_results, Day_specific = day_specific_results,
    Group_day_interaction = interaction_results, Method_comparison = method_comparison, Prevalence = prevalence,
    Pseudocounts = pseudocount_table, Selected_feature_map = selected_annotation_export, Sample_metadata = metadata_export,
    Analysis_summary = analysis_summary), path = project_output(file.path(out_dir, paste0("HILIC_exvivo_limma_",
    analysis_set, "_results.xlsx"))))
readr::write_csv(primary_results, project_output(file.path(out_dir, paste0("HILIC_exvivo_", analysis_set,
    "_primary_duplicateCorrelation.csv"))))
readr::write_csv(primary_results, project_output(file.path(out_dir, paste0("HILIC_exvivo_", analysis_set,
    "_unaveraged_d1_d4_duplicateCorrelation.csv"))))
readr::write_csv(mode_separated_primary_results, project_output(file.path(out_dir, paste0("HILIC_exvivo_",
    analysis_set, "_mode_separated_primary_globalFDR.csv"))))
readr::write_csv(mode_separated_positive_results, project_output(file.path(out_dir, paste0("HILIC_exvivo_",
    analysis_set, "_mode_separated_positive.csv"))))
readr::write_csv(mode_separated_negative_results, project_output(file.path(out_dir, paste0("HILIC_exvivo_",
    analysis_set, "_mode_separated_negative.csv"))))
readr::write_csv(mode_separated_vs_combined, project_output(file.path(out_dir, paste0("HILIC_exvivo_", analysis_set,
    "_mode_separated_vs_combined.csv"))))
readr::write_csv(donor_results, project_output(file.path(out_dir, paste0("HILIC_exvivo_", analysis_set, "_donor_averaged.csv"))))
readr::write_csv(baseline_results, project_output(file.path(out_dir, paste0("HILIC_exvivo_", analysis_set,
    "_baseline_adjusted.csv"))))
readr::write_csv(method_comparison, project_output(file.path(out_dir, paste0("HILIC_exvivo_", analysis_set,
    "_method_comparison.csv"))))
saveRDS(list(primary_results = primary_results, mode_separated_primary_results = mode_separated_primary_results,
    mode_separated_positive_results = mode_separated_positive_results, mode_separated_negative_results = mode_separated_negative_results,
    mode_separated_vs_combined = mode_separated_vs_combined, donor_results = donor_results, baseline_results = baseline_results,
    day_specific_results = day_specific_results, interaction_results = interaction_results, method_comparison = method_comparison,
    prevalence = prevalence, pseudocounts = pseudocount_table, metadata = metadata, features_log2 = features_log2,
    selected_feature_map = selected_annotation_export, settings = analysis_summary), project_output(file.path(out_dir,
    paste0("HILIC_exvivo_limma_", analysis_set, "_objects.rds"))))
message("Completed HILIC limma analysis for ", analysis_set, ". Metabolites tested: ", length(kept_metabolites),
    ".")
}))
