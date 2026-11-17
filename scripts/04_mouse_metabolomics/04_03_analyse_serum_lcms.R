# SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0
# Serum limma

source("scripts/helpers/setup.R")
invisible(run_analysis(
  section = "04_mouse_metabolomics",
  directory = "metabolomics",
  script = "04_03_analyse_serum_lcms.R",
  restricted = FALSE,
  code = {
required_packages <- c("limma", "statmod", "dplyr", "tibble", "stringr", "readr", "writexl")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]
if (length(missing_packages) > 0) {
    stop("Install these packages first: ", paste(missing_packages, collapse = ", "), "\nFor limma, use BiocManager::install('limma').")
}
input_rds <- project_path("metabolomics/Serum_metabolomics/20260121_data/Serum_cleaned_after_collaborator_review/Serum_LCMS_limma_input_objects.rds")
out_dir <- project_path("metabolomics/Serum_metabolomics/20260121_data/Serum_cleaned_after_collaborator_review/Limma_analysis")
group_col <- "Separation"
sample_col <- "Sample_Name"
group_levels <- c("ctrl", "ppi")
donor_col <- NULL
dir.create(project_output(out_dir, directory = TRUE), recursive = TRUE, showWarnings = FALSE)
if (!project_file_exists(input_rds)) {
    stop("Input RDS not found:\n", input_rds)
}
clean_obj <- readRDS(project_input(input_rds))
required_objects <- c("serum_clean_bio", "metadata_clean_bio")
missing_objects <- setdiff(required_objects, names(clean_obj))
if (length(missing_objects) > 0) {
    stop("The RDS is missing: ", paste(missing_objects, collapse = ", "))
}
features_raw <- as.matrix(clean_obj$serum_clean_bio)
metadata <- as.data.frame(clean_obj$metadata_clean_bio, check.names = FALSE)
storage.mode(features_raw) <- "numeric"
if (!sample_col %in% names(metadata)) {
    stop("Metadata does not contain sample column: ", sample_col)
}
if (!group_col %in% names(metadata)) {
    stop("Metadata does not contain group column: ", group_col)
}
if (is.null(rownames(features_raw))) {
    stop("The feature matrix has no sample row names.")
}
if (anyDuplicated(rownames(features_raw))) {
    stop("The feature matrix contains duplicated sample row names.")
}
if (anyDuplicated(metadata[[sample_col]])) {
    stop("Metadata contains duplicated sample names.")
}
if (!identical(rownames(features_raw), as.character(metadata[[sample_col]]))) {
    sample_match <- match(rownames(features_raw), as.character(metadata[[sample_col]]))
    if (anyNA(sample_match)) {
        stop("Feature-matrix row names cannot be matched to metadata sample names.")
    }
    metadata <- metadata[sample_match, , drop = FALSE]
}
stopifnot(identical(rownames(features_raw), as.character(metadata[[sample_col]])))
if (is.null(colnames(features_raw))) {
    stop("The feature matrix has no metabolite column names.")
}
if (anyNA(colnames(features_raw)) || any(colnames(features_raw) == "")) {
    stop("At least one metabolite name is missing or empty.")
}
if (anyDuplicated(colnames(features_raw))) {
    duplicate_names <- unique(colnames(features_raw)[duplicated(colnames(features_raw))])
    stop("Duplicate metabolite names remain in serum_clean_bio: ", paste(head(duplicate_names, 20), collapse = ", "),
        if (length(duplicate_names) > 20)
            " ..."
        else "", "\nRerun the corrected serum-cleaning script first.")
}
metadata[[group_col]] <- factor(as.character(metadata[[group_col]]), levels = group_levels)
if (anyNA(metadata[[group_col]])) {
    bad_groups <- unique(as.character(clean_obj$metadata_clean_bio[[group_col]][is.na(metadata[[group_col]])]))
    stop("Unexpected or missing group labels. Expected: ", paste(group_levels, collapse = ", "), ". Observed problematic values: ",
        paste(bad_groups, collapse = ", "))
}
if (is.null(donor_col) || !donor_col %in% names(metadata)) {
    metadata$Donor <- stringr::str_extract(as.character(metadata[[sample_col]]), "(?:ctrl|ppi)[0-9]+")
    donor_col_use <- "Donor"
} else {
    donor_col_use <- donor_col
}
metadata[[donor_col_use]] <- factor(as.character(metadata[[donor_col_use]]))
if (anyNA(metadata[[donor_col_use]])) {
    stop("At least one donor ID is missing. Check donor_col or Sample_Name.")
}
donor_group_check <- metadata %>% dplyr::distinct(Donor_check = .data[[donor_col_use]], Group_check = .data[[group_col]]) %>%
    dplyr::count(Donor_check, name = "n_groups") %>% dplyr::filter(n_groups != 1L)
if (nrow(donor_group_check) > 0) {
    stop("At least one donor maps to more than one treatment group.")
}
donor_counts <- metadata %>% dplyr::distinct(Donor = .data[[donor_col_use]], Group = .data[[group_col]]) %>%
    dplyr::count(Group, name = "n_donors")
print(donor_counts)
if (nrow(donor_counts) != 2 || any(donor_counts$n_donors != 3)) {
    warning("Expected three control and three PPI donors, but obtained:\n", paste(paste0(donor_counts$Group,
        "=", donor_counts$n_donors), collapse = ", "))
}
group_factor <- metadata[[group_col]]
usable_feature <- apply(features_raw, 2, function(x) {
    variable <- is.finite(stats::sd(x, na.rm = TRUE)) && stats::sd(x, na.rm = TRUE) > 0
    n_per_group <- tapply(is.finite(x), group_factor, sum)
    enough_data <- all(n_per_group >= 2)
    isTRUE(variable) && isTRUE(enough_data)
})
excluded_unusable <- tibble::tibble(metabolite = colnames(features_raw), retained = usable_feature) %>%
    dplyr::filter(!retained)
features_use <- features_raw[, usable_feature, drop = FALSE]
if (ncol(features_use) == 0) {
    stop("No usable metabolites remain after filtering.")
}
pseudocounts <- apply(features_use, 2, function(x) {
    positive <- x[is.finite(x) & x > 0]
    if (length(positive) == 0) {
        1e-12
    }
    else {
        min(positive)/2
    }
})
features_log <- log2(sweep(features_use, 2, pseudocounts, "+"))
rownames(features_log) <- rownames(features_use)
colnames(features_log) <- colnames(features_use)
extract_limma_results <- function(fit_ebayes, coef_name, method_name) {
    coef_index <- match(coef_name, colnames(fit_ebayes$coefficients))
    if (is.na(coef_index)) {
        stop("Coefficient not found: ", coef_name, ". Available coefficients: ", paste(colnames(fit_ebayes$coefficients),
            collapse = ", "))
    }
    estimate <- fit_ebayes$coefficients[, coef_index]
    moderated_se <- fit_ebayes$stdev.unscaled[, coef_index] * sqrt(fit_ebayes$s2.post)
    df_total <- fit_ebayes$df.total
    critical_t <- stats::qt(0.975, df = df_total)
    conf_low <- estimate - critical_t * moderated_se
    conf_high <- estimate + critical_t * moderated_se
    tibble::tibble(metabolite = rownames(fit_ebayes$coefficients), method = method_name, log2_FC = as.numeric(estimate),
        fold_change = 2^log2_FC, moderated_SE = as.numeric(moderated_se), df_total = as.numeric(df_total),
        conf_low = as.numeric(conf_low), conf_high = as.numeric(conf_high), fc_ci_low = 2^conf_low, fc_ci_high = 2^conf_high,
        AveExpr = as.numeric(fit_ebayes$Amean), moderated_t = as.numeric(fit_ebayes$t[, coef_index]),
        p_value = as.numeric(fit_ebayes$p.value[, coef_index]), FDR = stats::p.adjust(p_value, method = "BH"),
        B = as.numeric(fit_ebayes$lods[, coef_index])) %>% dplyr::arrange(p_value)
}
donor_ids <- levels(droplevels(metadata[[donor_col_use]]))
donor_by_metabolite <- vapply(donor_ids, function(current_donor) {
    donor_rows <- metadata[[donor_col_use]] == current_donor
    values <- colMeans(features_log[donor_rows, , drop = FALSE], na.rm = TRUE)
    values[is.nan(values)] <- NA_real_
    values
}, FUN.VALUE = numeric(ncol(features_log)))
Y_donor <- donor_by_metabolite
rownames(Y_donor) <- colnames(features_log)
colnames(Y_donor) <- donor_ids
donor_metadata <- metadata %>% dplyr::transmute(Donor = as.character(.data[[donor_col_use]]), Separation = .data[[group_col]]) %>%
    dplyr::distinct() %>% dplyr::slice(match(donor_ids, Donor))
if (!identical(donor_metadata$Donor, colnames(Y_donor))) {
    stop("Donor-level feature matrix and donor metadata are misaligned.")
}
donor_metadata$Separation <- factor(donor_metadata$Separation, levels = group_levels)
design_donor <- stats::model.matrix(~Separation, data = donor_metadata)
coef_name <- paste0("Separation", group_levels[2])
fit_donor <- limma::lmFit(Y_donor, design_donor)
fit_donor_ebayes <- limma::eBayes(fit_donor, robust = TRUE)
results_donor <- extract_limma_results(fit_ebayes = fit_donor_ebayes, coef_name = coef_name, method_name = "Donor-aggregated limma")
Y_mouse <- t(features_log)
design_mouse <- stats::model.matrix(~Separation, data = data.frame(Separation = metadata[[group_col]]))
block_mouse <- metadata[[donor_col_use]]
corfit <- limma::duplicateCorrelation(object = Y_mouse, design = design_mouse, block = block_mouse)
consensus_correlation <- corfit$consensus.correlation
if (is.null(consensus_correlation)) {
    consensus_correlation <- corfit$consensus
}
if (length(consensus_correlation) != 1 || !is.finite(consensus_correlation)) {
    stop("duplicateCorrelation did not return a finite consensus correlation.")
}
fit_mouse <- limma::lmFit(object = Y_mouse, design = design_mouse, block = block_mouse, correlation = consensus_correlation)
fit_mouse_ebayes <- limma::eBayes(fit_mouse, robust = TRUE)
results_duplicate_correlation <- extract_limma_results(fit_ebayes = fit_mouse_ebayes, coef_name = coef_name,
    method_name = "Mouse-level duplicateCorrelation limma")
comparison <- results_donor %>% dplyr::select(metabolite, donor_log2_FC = log2_FC, donor_fold_change = fold_change,
    donor_conf_low = conf_low, donor_conf_high = conf_high, donor_p_value = p_value, donor_FDR = FDR) %>%
    dplyr::full_join(results_duplicate_correlation %>% dplyr::select(metabolite, duplicateCorrelation_log2_FC = log2_FC,
        duplicateCorrelation_fold_change = fold_change, duplicateCorrelation_conf_low = conf_low, duplicateCorrelation_conf_high = conf_high,
        duplicateCorrelation_p_value = p_value, duplicateCorrelation_FDR = FDR), by = "metabolite") %>%
    dplyr::mutate(same_direction = sign(donor_log2_FC) == sign(duplicateCorrelation_log2_FC), absolute_log2FC_difference = abs(donor_log2_FC -
        duplicateCorrelation_log2_FC), donor_FDR_below_0_10 = donor_FDR < 0.1, duplicateCorrelation_FDR_below_0_10 = duplicateCorrelation_FDR <
        0.1) %>% dplyr::arrange(donor_p_value)
diagnostics <- tibble::tibble(item = c("Number of recipient mice", "Number of independent donors", "Number of metabolites before usability filtering",
    "Number of metabolites analyzed", "Number of unusable metabolites excluded", "duplicateCorrelation consensus correlation",
    "Donor-level limma FDR < 0.05", "Donor-level limma FDR < 0.10", "duplicateCorrelation limma FDR < 0.05",
    "duplicateCorrelation limma FDR < 0.10", "Effect-direction agreement"), value = c(nrow(features_raw),
    nrow(donor_metadata), ncol(features_raw), ncol(features_log), nrow(excluded_unusable), consensus_correlation,
    sum(results_donor$FDR < 0.05, na.rm = TRUE), sum(results_donor$FDR < 0.1, na.rm = TRUE), sum(results_duplicate_correlation$FDR <
        0.05, na.rm = TRUE), sum(results_duplicate_correlation$FDR < 0.1, na.rm = TRUE), mean(comparison$same_direction,
        na.rm = TRUE)))
donor_profile_export <- data.frame(metabolite = rownames(Y_donor), Y_donor, check.names = FALSE)
pseudocount_export <- tibble::tibble(metabolite = names(pseudocounts), pseudocount = as.numeric(pseudocounts))
output_xlsx <- file.path(out_dir, "Serum_limma_donor_and_duplicateCorrelation_results.xlsx")
writexl::write_xlsx(list(Diagnostics = diagnostics, Donor_level_limma = results_donor, Mouse_duplicateCorrelation = results_duplicate_correlation,
    Method_comparison = comparison, Donor_metadata = donor_metadata, Donor_log2_profiles = donor_profile_export,
    Pseudocounts = pseudocount_export, Excluded_unusable = excluded_unusable), path = project_output(output_xlsx))
readr::write_csv(results_donor, project_output(file.path(out_dir, "Serum_limma_donor_level_results.csv"),
    name = "serum_limma_donor_level_results.csv"))
readr::write_csv(results_duplicate_correlation, project_output(file.path(out_dir, "Serum_limma_duplicateCorrelation_results.csv"),
    name = "serum_limma_duplicate_correlation_results.csv"))
readr::write_csv(comparison, project_output(file.path(out_dir, "Serum_limma_method_comparison.csv"), name = "serum_limma_method_comparison.csv"))
saveRDS(list(input_rds = input_rds, metadata_mouse = metadata, features_raw = features_raw, features_log = features_log,
    pseudocounts = pseudocounts, donor_metadata = donor_metadata, Y_donor = Y_donor, Y_mouse = Y_mouse,
    design_donor = design_donor, design_mouse = design_mouse, duplicateCorrelation = corfit, consensus_correlation = consensus_correlation,
    fit_donor = fit_donor, fit_donor_ebayes = fit_donor_ebayes, fit_mouse = fit_mouse, fit_mouse_ebayes = fit_mouse_ebayes,
    results_donor = results_donor, results_duplicate_correlation = results_duplicate_correlation, comparison = comparison,
    diagnostics = diagnostics), project_output(file.path(out_dir, "Serum_limma_donor_and_duplicateCorrelation_objects.rds"),
    name = "serum_limma_donor_and_duplicate_correlation_objects.rds"))
message("\nSerum limma analyses completed.\n", "Recipient mice: ", nrow(features_raw), "\nIndependent donors: ",
    nrow(donor_metadata), "\nMetabolites analyzed: ", ncol(features_log), "\nConsensus within-donor correlation: ",
    signif(consensus_correlation, 4), "\nDonor-level FDR < 0.10: ", sum(results_donor$FDR < 0.1, na.rm = TRUE),
    "\nMouse-level duplicateCorrelation FDR < 0.10: ", sum(results_duplicate_correlation$FDR < 0.1, na.rm = TRUE),
    "\nEffect-direction agreement: ", round(100 * mean(comparison$same_direction, na.rm = TRUE), 1),
    "%\n\nMain output:\n", output_xlsx)
}))
