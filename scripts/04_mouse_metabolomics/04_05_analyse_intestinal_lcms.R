# SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0
# Caecum/faeces limma

source("scripts/helpers/setup.R")
invisible(run_analysis(
  section = "04_mouse_metabolomics",
  directory = "metabolomics",
  script = "04_05_analyse_intestinal_lcms.R",
  restricted = FALSE,
  code = {
required_packages <- c("limma", "statmod", "dplyr", "tidyr", "purrr", "tibble", "stringr", "readr", "writexl")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]
if (length(missing_packages) > 0) {
    stop("Install these packages first: ", paste(missing_packages, collapse = ", "), "\nFor limma, use BiocManager::install('limma').")
}
input_rds <- project_path("metabolomics/20250714_MassSpecData/Cecum_Feces_cleaned_after_collaborator_review/Cecum_Feces_LCMS_limma_input_objects.rds")
out_dir <- project_path("metabolomics/20250714_MassSpecData/Cecum_Feces_cleaned_after_collaborator_review/Limma_analysis")
dir.create(project_output(out_dir, directory = TRUE), recursive = TRUE, showWarnings = FALSE)
organs_to_analyze <- c("cecum", "feces")
representations_to_analyze <- c("rt_resolved", "one_per_metabolite")
group_col <- "Separation"
sample_col <- "Sample_Name"
donor_col <- "Donor"
mouse_col <- "Donor_Mouse_ID"
group_levels <- c("ctrl", "ppi")
min_positive_per_group <- 4L
min_donors_detected_per_group <- 2L
robust_ebayes <- TRUE
if (!project_file_exists(input_rds)) {
    stop("Input RDS not found:\n", input_rds)
}
clean_obj <- readRDS(project_input(input_rds))
required_top_level <- c("rt_resolved", "one_per_metabolite", "metadata")
missing_top_level <- setdiff(required_top_level, names(clean_obj))
if (length(missing_top_level) > 0) {
    stop("The input RDS is missing: ", paste(missing_top_level, collapse = ", "))
}
unknown_organs <- setdiff(organs_to_analyze, names(clean_obj$metadata))
if (length(unknown_organs) > 0) {
    stop("Metadata are missing for organ(s): ", paste(unknown_organs, collapse = ", "))
}
unknown_representations <- setdiff(representations_to_analyze, names(clean_obj))
if (length(unknown_representations) > 0) {
    stop("The input RDS lacks representation(s): ", paste(unknown_representations, collapse = ", "))
}
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
    tibble::tibble(analyte = rownames(fit_ebayes$coefficients), method = method_name, log2_FC = as.numeric(estimate),
        fold_change = 2^log2_FC, moderated_SE = as.numeric(moderated_se), df_total = as.numeric(df_total),
        conf_low = as.numeric(conf_low), conf_high = as.numeric(conf_high), fc_ci_low = 2^conf_low, fc_ci_high = 2^conf_high,
        AveExpr = as.numeric(fit_ebayes$Amean), moderated_t = as.numeric(fit_ebayes$t[, coef_index]),
        p_value = as.numeric(fit_ebayes$p.value[, coef_index]), FDR = stats::p.adjust(p_value, method = "BH"),
        B = as.numeric(fit_ebayes$lods[, coef_index])) %>% dplyr::arrange(p_value)
}
validate_and_prepare <- function(features_raw, metadata, organ, representation) {
    features_raw <- as.matrix(features_raw)
    metadata <- as.data.frame(metadata, check.names = FALSE)
    storage.mode(features_raw) <- "numeric"
    required_metadata <- c(sample_col, group_col, donor_col, mouse_col, "Organ")
    missing_metadata <- setdiff(required_metadata, names(metadata))
    if (length(missing_metadata) > 0) {
        stop(organ, "/", representation, " metadata are missing: ", paste(missing_metadata, collapse = ", "))
    }
    if (is.null(rownames(features_raw)) || anyDuplicated(rownames(features_raw))) {
        stop(organ, "/", representation, ": invalid or duplicated matrix row names.")
    }
    if (is.null(colnames(features_raw)) || anyNA(colnames(features_raw)) || any(colnames(features_raw) ==
        "") || anyDuplicated(colnames(features_raw))) {
        stop(organ, "/", representation, ": invalid or duplicated analyte names.")
    }
    if (anyDuplicated(metadata[[sample_col]])) {
        stop(organ, ": metadata contain duplicated sample names.")
    }
    if (!identical(rownames(features_raw), as.character(metadata[[sample_col]]))) {
        sample_match <- match(rownames(features_raw), as.character(metadata[[sample_col]]))
        if (anyNA(sample_match)) {
            stop(organ, ": matrix rows cannot be matched to metadata.")
        }
        metadata <- metadata[sample_match, , drop = FALSE]
    }
    stopifnot(identical(rownames(features_raw), as.character(metadata[[sample_col]])))
    if (any(features_raw < 0, na.rm = TRUE)) {
        stop(organ, "/", representation, ": negative values detected. Inspect the preprocessing rather than ",
            "silently truncating LC-MS intensities.")
    }
    metadata[[group_col]] <- factor(stringr::str_to_lower(as.character(metadata[[group_col]])), levels = group_levels)
    metadata[[donor_col]] <- factor(as.character(metadata[[donor_col]]))
    if (anyNA(metadata[[group_col]]) || anyNA(metadata[[donor_col]])) {
        stop(organ, ": missing/unexpected group or donor labels.")
    }
    if (!all(stringr::str_to_lower(metadata$Organ) == organ)) {
        stop("Metadata organ mismatch for ", organ, ".")
    }
    donor_group_check <- metadata %>% dplyr::distinct(Donor_check = .data[[donor_col]], Group_check = .data[[group_col]]) %>%
        dplyr::count(Donor_check, name = "n_groups") %>% dplyr::filter(n_groups != 1L)
    if (nrow(donor_group_check) > 0) {
        stop(organ, ": at least one donor maps to more than one group.")
    }
    group_factor <- metadata[[group_col]]
    usable_feature <- apply(features_raw, 2, function(x) {
        x_finite <- x[is.finite(x)]
        variable <- length(x_finite) >= 2 && is.finite(stats::sd(x_finite)) && stats::sd(x_finite) >
            0
        n_positive_per_group <- tapply(is.finite(x) & x > 0, group_factor, sum)
        enough_positive <- all(n_positive_per_group >= min_positive_per_group)
        n_donors_detected_per_group <- vapply(group_levels, function(current_group) {
            donors_in_group <- unique(as.character(metadata[[donor_col]][group_factor == current_group]))
            sum(vapply(donors_in_group, function(current_donor) {
                current_rows <- group_factor == current_group & as.character(metadata[[donor_col]]) ==
                  current_donor
                any(is.finite(x[current_rows]) & x[current_rows] > 0)
            }, FUN.VALUE = logical(1)))
        }, FUN.VALUE = integer(1))
        enough_donor_coverage <- all(n_donors_detected_per_group >= min_donors_detected_per_group)
        isTRUE(variable) && isTRUE(enough_positive) && isTRUE(enough_donor_coverage)
    })
    excluded_unusable <- tibble::tibble(analyte = colnames(features_raw), retained = usable_feature) %>%
        dplyr::filter(!retained)
    features_use <- features_raw[, usable_feature, drop = FALSE]
    if (ncol(features_use) == 0) {
        stop(organ, "/", representation, ": no usable analytes remain.")
    }
    pseudocounts <- apply(features_use, 2, function(x) {
        positive <- x[is.finite(x) & x > 0]
        if (length(positive) == 0)
            1e-12
        else min(positive)/2
    })
    features_log <- log2(sweep(features_use, 2, pseudocounts, "+"))
    rownames(features_log) <- rownames(features_use)
    colnames(features_log) <- colnames(features_use)
    list(features_raw = features_raw, metadata = metadata, features_log = features_log, pseudocounts = pseudocounts,
        excluded_unusable = excluded_unusable)
}
run_one_analysis <- function(features_raw, metadata, organ, representation) {
    prepared <- validate_and_prepare(features_raw = features_raw, metadata = metadata, organ = organ,
        representation = representation)
    features_log <- prepared$features_log
    metadata <- prepared$metadata
    donor_ids <- levels(droplevels(metadata[[donor_col]]))
    donor_by_analyte <- vapply(donor_ids, function(current_donor) {
        donor_rows <- metadata[[donor_col]] == current_donor
        values <- colMeans(features_log[donor_rows, , drop = FALSE], na.rm = TRUE)
        values[is.nan(values)] <- NA_real_
        values
    }, FUN.VALUE = numeric(ncol(features_log)))
    Y_donor <- donor_by_analyte
    rownames(Y_donor) <- colnames(features_log)
    colnames(Y_donor) <- donor_ids
    donor_metadata <- metadata %>% dplyr::transmute(Donor = as.character(.data[[donor_col]]), Separation = .data[[group_col]]) %>%
        dplyr::distinct() %>% dplyr::slice(match(donor_ids, Donor))
    if (!identical(donor_metadata$Donor, colnames(Y_donor))) {
        stop(organ, "/", representation, ": donor profiles and metadata are misaligned.")
    }
    donor_metadata$Separation <- factor(donor_metadata$Separation, levels = group_levels)
    design_donor <- stats::model.matrix(~Separation, data = donor_metadata)
    coef_name <- paste0("Separation", group_levels[2])
    fit_donor <- limma::lmFit(Y_donor, design_donor)
    fit_donor_ebayes <- limma::eBayes(fit_donor, robust = robust_ebayes)
    results_donor <- extract_limma_results(fit_donor_ebayes, coef_name, "Donor-aggregated limma") %>%
        dplyr::mutate(organ = organ, representation = representation, .before = 1)
    Y_mouse <- t(features_log)
    design_mouse <- stats::model.matrix(~Separation, data = data.frame(Separation = metadata[[group_col]]))
    block_mouse <- metadata[[donor_col]]
    corfit <- limma::duplicateCorrelation(object = Y_mouse, design = design_mouse, block = block_mouse)
    consensus_correlation <- corfit$consensus.correlation
    if (is.null(consensus_correlation)) {
        consensus_correlation <- corfit$consensus
    }
    if (length(consensus_correlation) != 1 || !is.finite(consensus_correlation)) {
        stop(organ, "/", representation, ": duplicateCorrelation did not return a finite consensus correlation.")
    }
    fit_mouse <- limma::lmFit(object = Y_mouse, design = design_mouse, block = block_mouse, correlation = consensus_correlation)
    fit_mouse_ebayes <- limma::eBayes(fit_mouse, robust = robust_ebayes)
    results_mouse <- extract_limma_results(fit_mouse_ebayes, coef_name, "Mouse-level duplicateCorrelation limma") %>%
        dplyr::mutate(organ = organ, representation = representation, .before = 1)
    comparison <- results_donor %>% dplyr::select(organ, representation, analyte, donor_log2_FC = log2_FC,
        donor_fold_change = fold_change, donor_conf_low = conf_low, donor_conf_high = conf_high, donor_p_value = p_value,
        donor_FDR = FDR) %>% dplyr::full_join(results_mouse %>% dplyr::select(organ, representation,
        analyte, mouse_log2_FC = log2_FC, mouse_fold_change = fold_change, mouse_conf_low = conf_low,
        mouse_conf_high = conf_high, mouse_p_value = p_value, mouse_FDR = FDR), by = c("organ", "representation",
        "analyte")) %>% dplyr::mutate(same_direction = sign(donor_log2_FC) == sign(mouse_log2_FC), absolute_log2FC_difference = abs(donor_log2_FC -
        mouse_log2_FC), donor_FDR_below_0_10 = donor_FDR < 0.1, mouse_FDR_below_0_10 = mouse_FDR < 0.1) %>%
        dplyr::arrange(donor_p_value)
    diagnostics <- tibble::tibble(organ = organ, representation = representation, item = c("Recipient mice",
        "Independent donors", "Analytes before usability filtering", "Analytes analyzed", "Unusable analytes excluded",
        "duplicateCorrelation consensus correlation", "Donor limma FDR < 0.05", "Donor limma FDR < 0.10",
        "Mouse duplicateCorrelation FDR < 0.05", "Mouse duplicateCorrelation FDR < 0.10", "Effect-direction agreement"),
        value = c(nrow(prepared$features_raw), nrow(donor_metadata), ncol(prepared$features_raw), ncol(features_log),
            nrow(prepared$excluded_unusable), consensus_correlation, sum(results_donor$FDR < 0.05, na.rm = TRUE),
            sum(results_donor$FDR < 0.1, na.rm = TRUE), sum(results_mouse$FDR < 0.05, na.rm = TRUE),
            sum(results_mouse$FDR < 0.1, na.rm = TRUE), mean(comparison$same_direction, na.rm = TRUE)))
    donor_profile_export <- data.frame(analyte = rownames(Y_donor), Y_donor, check.names = FALSE)
    pseudocount_export <- tibble::tibble(analyte = names(prepared$pseudocounts), pseudocount = as.numeric(prepared$pseudocounts))
    list(organ = organ, representation = representation, metadata_mouse = metadata, features_raw = prepared$features_raw,
        features_log = features_log, pseudocounts = prepared$pseudocounts, pseudocount_export = pseudocount_export,
        excluded_unusable = prepared$excluded_unusable, donor_metadata = donor_metadata, Y_donor = Y_donor,
        donor_profile_export = donor_profile_export, Y_mouse = Y_mouse, design_donor = design_donor,
        design_mouse = design_mouse, duplicateCorrelation = corfit, consensus_correlation = consensus_correlation,
        fit_donor = fit_donor, fit_donor_ebayes = fit_donor_ebayes, fit_mouse = fit_mouse, fit_mouse_ebayes = fit_mouse_ebayes,
        results_donor = results_donor, results_mouse = results_mouse, comparison = comparison, diagnostics = diagnostics)
}
make_sheet_prefix <- function(organ, representation) {
    organ_short <- ifelse(organ == "cecum", "Cec", "Fec")
    rep_short <- ifelse(representation == "rt_resolved", "RT", "Met")
    paste0(organ_short, "_", rep_short)
}
analysis_grid <- tidyr::expand_grid(organ = organs_to_analyze, representation = representations_to_analyze)
analysis_results <- purrr::pmap(analysis_grid, function(organ, representation) {
    if (!organ %in% names(clean_obj[[representation]])) {
        stop(representation, " matrix is missing for organ: ", organ)
    }
    run_one_analysis(features_raw = clean_obj[[representation]][[organ]], metadata = clean_obj$metadata[[organ]],
        organ = organ, representation = representation)
})
names(analysis_results) <- paste(analysis_grid$organ, analysis_grid$representation, sep = "__")
combined_donor <- purrr::map_dfr(analysis_results, "results_donor")
combined_mouse <- purrr::map_dfr(analysis_results, "results_mouse")
combined_comparison <- purrr::map_dfr(analysis_results, "comparison")
combined_diagnostics <- purrr::map_dfr(analysis_results, "diagnostics")
cross_organ_donor_comparison <- combined_donor %>% dplyr::select(representation, organ, analyte, log2_FC,
    p_value, FDR) %>% tidyr::pivot_wider(names_from = organ, values_from = c(log2_FC, p_value, FDR),
    names_glue = "{organ}_{.value}") %>% dplyr::mutate(same_direction = sign(cecum_log2_FC) == sign(feces_log2_FC),
    absolute_log2FC_difference = abs(cecum_log2_FC - feces_log2_FC)) %>% dplyr::arrange(representation,
    cecum_p_value)
workbook_sheets <- list(Diagnostics = combined_diagnostics, All_donor_results = combined_donor, All_mouse_results = combined_mouse,
    All_method_comparison = combined_comparison, Cross_organ_donor = cross_organ_donor_comparison)
for (i in seq_along(analysis_results)) {
    current <- analysis_results[[i]]
    prefix <- make_sheet_prefix(current$organ, current$representation)
    workbook_sheets[[paste0(prefix, "_Donor")]] <- current$results_donor
    workbook_sheets[[paste0(prefix, "_Mouse")]] <- current$results_mouse
    workbook_sheets[[paste0(prefix, "_Compare")]] <- current$comparison
    workbook_sheets[[paste0(prefix, "_DonorProf")]] <- current$donor_profile_export
    workbook_sheets[[paste0(prefix, "_Pseudo")]] <- current$pseudocount_export
    workbook_sheets[[paste0(prefix, "_Excluded")]] <- current$excluded_unusable
}
output_xlsx <- file.path(out_dir, "Cecum_Feces_LCMS_limma_donor_and_duplicateCorrelation_results.xlsx")
writexl::write_xlsx(workbook_sheets, path = project_output(output_xlsx))
readr::write_csv(combined_donor, project_output(file.path(out_dir, "Cecum_Feces_LCMS_all_donor_level_results.csv"),
    name = "caecum_faeces_lcms_all_donor_level_results.csv"))
readr::write_csv(combined_mouse, project_output(file.path(out_dir, "Cecum_Feces_LCMS_all_duplicateCorrelation_results.csv"),
    name = "caecum_faeces_lcms_all_duplicate_correlation_results.csv"))
readr::write_csv(combined_comparison, project_output(file.path(out_dir, "Cecum_Feces_LCMS_all_method_comparison.csv"),
    name = "caecum_faeces_lcms_all_method_comparison.csv"))
readr::write_csv(cross_organ_donor_comparison, project_output(file.path(out_dir, "Cecum_Feces_LCMS_cross_organ_donor_comparison.csv"),
    name = "caecum_faeces_lcms_cross_organ_donor_comparison.csv"))
for (i in seq_along(analysis_results)) {
    current <- analysis_results[[i]]
    stem <- paste(current$organ, current$representation, sep = "_")
    readr::write_csv(current$results_donor, project_output(file.path(out_dir, paste0(stem, "_donor_level_results.csv"))))
    readr::write_csv(current$results_mouse, project_output(file.path(out_dir, paste0(stem, "_duplicateCorrelation_results.csv"))))
    readr::write_csv(current$comparison, project_output(file.path(out_dir, paste0(stem, "_method_comparison.csv"))))
}
saveRDS(list(input_rds = input_rds, settings = list(organs_to_analyze = organs_to_analyze, representations_to_analyze = representations_to_analyze,
    group_levels = group_levels, min_positive_per_group = min_positive_per_group, min_donors_detected_per_group = min_donors_detected_per_group,
    robust_ebayes = robust_ebayes), analysis_grid = analysis_grid, analysis_results = analysis_results,
    combined_donor = combined_donor, combined_mouse = combined_mouse, combined_comparison = combined_comparison,
    combined_diagnostics = combined_diagnostics, cross_organ_donor_comparison = cross_organ_donor_comparison),
    project_output(file.path(out_dir, "Cecum_Feces_LCMS_limma_donor_and_duplicateCorrelation_objects.rds"),
        name = "caecum_faeces_lcms_limma_donor_and_duplicate_correlation_objects.rds"))
message("\nCaecum/faeces LC-MS limma analyses completed.\n", "Combinations analysed: ", nrow(analysis_grid),
    "\nOrgans: ", paste(organs_to_analyze, collapse = ", "), "\nRepresentations: ", paste(representations_to_analyze,
        collapse = ", "), "\n\nDonor-level FDR < 0.10 by analysis:\n", paste(combined_donor %>% dplyr::group_by(organ,
        representation) %>% dplyr::summarise(n = sum(FDR < 0.1, na.rm = TRUE), .groups = "drop") %>%
        dplyr::transmute(text = paste0(organ, " / ", representation, ": ", n)) %>% dplyr::pull(text),
        collapse = "\n"), "\n\nMouse duplicateCorrelation FDR < 0.10 by analysis:\n", paste(combined_mouse %>%
        dplyr::group_by(organ, representation) %>% dplyr::summarise(n = sum(FDR < 0.1, na.rm = TRUE),
        .groups = "drop") %>% dplyr::transmute(text = paste0(organ, " / ", representation, ": ", n)) %>%
        dplyr::pull(text), collapse = "\n"), "\n\nMain output:\n", output_xlsx)
run_combined_compartment_analysis <- TRUE
combined_min_positive_per_group_compartment <- 4L
combined_min_donors_detected_per_group_compartment <- 2L
prepare_combined_compartment_data <- function(clean_obj, representation) {
    if (!all(c("cecum", "feces") %in% names(clean_obj[[representation]]))) {
        stop(representation, ": caecum and/or feces matrix is missing.")
    }
    X_cecum <- as.matrix(clean_obj[[representation]]$cecum)
    X_feces <- as.matrix(clean_obj[[representation]]$feces)
    storage.mode(X_cecum) <- "numeric"
    storage.mode(X_feces) <- "numeric"
    common_analytes <- intersect(colnames(X_cecum), colnames(X_feces))
    if (length(common_analytes) == 0L) {
        stop(representation, ": no analytes are shared between compartments.")
    }
    only_cecum <- setdiff(colnames(X_cecum), colnames(X_feces))
    only_feces <- setdiff(colnames(X_feces), colnames(X_cecum))
    if (length(only_cecum) > 0L || length(only_feces) > 0L) {
        warning(representation, ": the combined model uses only analytes present in both matrices. ",
            "Cecum-only = ", length(only_cecum), "; feces-only = ", length(only_feces), ".")
    }
    common_analytes <- colnames(X_cecum)[colnames(X_cecum) %in% common_analytes]
    X_cecum <- X_cecum[, common_analytes, drop = FALSE]
    X_feces <- X_feces[, common_analytes, drop = FALSE]
    metadata_cecum <- as.data.frame(clean_obj$metadata$cecum, check.names = FALSE)
    metadata_feces <- as.data.frame(clean_obj$metadata$feces, check.names = FALSE)
    X_all <- rbind(X_cecum, X_feces)
    metadata_all <- dplyr::bind_rows(metadata_cecum, metadata_feces)
    required_metadata <- c(sample_col, group_col, donor_col, mouse_col, "Organ")
    missing_metadata <- setdiff(required_metadata, names(metadata_all))
    if (length(missing_metadata) > 0L) {
        stop(representation, " combined metadata are missing: ", paste(missing_metadata, collapse = ", "))
    }
    if (anyDuplicated(metadata_all[[sample_col]])) {
        stop(representation, ": combined metadata contain duplicated sample names.")
    }
    sample_match <- match(rownames(X_all), as.character(metadata_all[[sample_col]]))
    if (anyNA(sample_match)) {
        stop(representation, ": combined matrix rows cannot be matched to metadata.")
    }
    metadata_all <- metadata_all[sample_match, , drop = FALSE]
    stopifnot(identical(rownames(X_all), as.character(metadata_all[[sample_col]])))
    if (any(X_all < 0, na.rm = TRUE)) {
        stop(representation, ": negative intensities detected in the combined matrix.")
    }
    metadata_all[[group_col]] <- factor(stringr::str_to_lower(as.character(metadata_all[[group_col]])),
        levels = group_levels)
    metadata_all[[donor_col]] <- factor(as.character(metadata_all[[donor_col]]))
    metadata_all[[mouse_col]] <- factor(as.character(metadata_all[[mouse_col]]))
    metadata_all$Organ <- factor(stringr::str_to_lower(as.character(metadata_all$Organ)), levels = c("cecum",
        "feces"))
    if (anyNA(metadata_all[[group_col]]) || anyNA(metadata_all[[donor_col]]) || anyNA(metadata_all[[mouse_col]]) ||
        anyNA(metadata_all$Organ)) {
        stop(representation, ": unexpected or missing combined metadata labels.")
    }
    pairing_check <- metadata_all %>% dplyr::count(mouse = .data[[mouse_col]], Organ, name = "n") %>%
        tidyr::pivot_wider(names_from = Organ, values_from = n, values_fill = 0)
    if (!all(c("cecum", "feces") %in% names(pairing_check)) || any(pairing_check$cecum != 1L) || any(pairing_check$feces !=
        1L)) {
        stop(representation, ": each mouse must have exactly one caecum and one feces sample.")
    }
    donor_group_check <- metadata_all %>% dplyr::distinct(donor = .data[[donor_col]], group = .data[[group_col]]) %>%
        dplyr::count(donor, name = "n_groups") %>% dplyr::filter(n_groups != 1L)
    if (nrow(donor_group_check) > 0L) {
        stop(representation, ": a donor maps to more than one treatment group.")
    }
    group_compartment <- interaction(metadata_all[[group_col]], metadata_all$Organ, drop = FALSE, sep = "__")
    usable_feature <- apply(X_all, 2, function(x) {
        finite_x <- x[is.finite(x)]
        variable <- length(finite_x) >= 2L && is.finite(stats::sd(finite_x)) && stats::sd(finite_x) >
            0
        positive_by_cell <- tapply(is.finite(x) & x > 0, group_compartment, sum)
        enough_positive <- length(positive_by_cell) == 4L && all(is.finite(positive_by_cell)) && all(positive_by_cell >=
            combined_min_positive_per_group_compartment)
        donor_coverage_by_cell <- vapply(levels(group_compartment), function(current_cell) {
            current_rows <- group_compartment == current_cell
            donors_in_cell <- unique(as.character(metadata_all[[donor_col]][current_rows]))
            sum(vapply(donors_in_cell, function(current_donor) {
                donor_rows <- current_rows & as.character(metadata_all[[donor_col]]) == current_donor
                any(is.finite(x[donor_rows]) & x[donor_rows] > 0)
            }, FUN.VALUE = logical(1)))
        }, FUN.VALUE = integer(1))
        enough_donor_coverage <- length(donor_coverage_by_cell) == 4L && all(donor_coverage_by_cell >=
            combined_min_donors_detected_per_group_compartment)
        isTRUE(variable) && isTRUE(enough_positive) && isTRUE(enough_donor_coverage)
    })
    excluded_unusable <- tibble::tibble(analyte = colnames(X_all), retained = usable_feature, reason = ifelse(usable_feature,
        NA_character_, paste0("Failed variability and/or prevalence rule: >=", combined_min_positive_per_group_compartment,
            " positive mice and >=", combined_min_donors_detected_per_group_compartment, " donors in every group x compartment cell"))) %>%
        dplyr::filter(!retained)
    X_use <- X_all[, usable_feature, drop = FALSE]
    if (ncol(X_use) == 0L) {
        stop(representation, ": no analytes pass the combined prevalence rule.")
    }
    pseudocounts <- apply(X_use, 2, function(x) {
        positive <- x[is.finite(x) & x > 0]
        if (length(positive) == 0L)
            1e-12
        else min(positive)/2
    })
    X_log <- log2(sweep(X_use, 2, pseudocounts, "+"))
    rownames(X_log) <- rownames(X_use)
    colnames(X_log) <- colnames(X_use)
    list(X_raw = X_all, X_log = X_log, metadata = metadata_all, pseudocounts = pseudocounts, excluded_unusable = excluded_unusable,
        only_cecum = only_cecum, only_feces = only_feces)
}
run_combined_one_representation <- function(clean_obj, representation) {
    prepared <- prepare_combined_compartment_data(clean_obj = clean_obj, representation = representation)
    X_log <- prepared$X_log
    metadata <- prepared$metadata
    Y_all <- t(X_log)
    coef_group <- paste0("Separation", group_levels[2])
    coef_compartment <- "Organfeces"
    donor_compartment_metadata <- metadata %>% dplyr::transmute(Donor = as.character(.data[[donor_col]]),
        Separation = .data[[group_col]], Organ = .data$Organ) %>% dplyr::distinct() %>% dplyr::arrange(Donor,
        Organ)
    donor_compartment_metadata$Separation <- factor(donor_compartment_metadata$Separation, levels = group_levels)
    donor_compartment_metadata$Organ <- factor(donor_compartment_metadata$Organ, levels = c("cecum",
        "feces"))
    donor_compartment_keys <- paste(donor_compartment_metadata$Donor, donor_compartment_metadata$Organ,
        sep = "__")
    Y_donor_compartment <- vapply(donor_compartment_keys, function(current_key) {
        current_rows <- paste(as.character(metadata[[donor_col]]), as.character(metadata$Organ), sep = "__") ==
            current_key
        values <- colMeans(X_log[current_rows, , drop = FALSE], na.rm = TRUE)
        values[is.nan(values)] <- NA_real_
        values
    }, FUN.VALUE = numeric(ncol(X_log)))
    rownames(Y_donor_compartment) <- colnames(X_log)
    colnames(Y_donor_compartment) <- donor_compartment_keys
    design_donor_compartment <- stats::model.matrix(~Separation + Organ, data = donor_compartment_metadata)
    block_donor_compartment <- factor(donor_compartment_metadata$Donor)
    corfit_donor_compartment <- limma::duplicateCorrelation(object = Y_donor_compartment, design = design_donor_compartment,
        block = block_donor_compartment)
    donor_compartment_correlation <- corfit_donor_compartment$consensus.correlation
    if (is.null(donor_compartment_correlation)) {
        donor_compartment_correlation <- corfit_donor_compartment$consensus
    }
    if (length(donor_compartment_correlation) != 1L || !is.finite(donor_compartment_correlation)) {
        stop(representation, ": donor-compartment duplicateCorrelation failed.")
    }
    fit_donor_compartment <- limma::lmFit(object = Y_donor_compartment, design = design_donor_compartment,
        block = block_donor_compartment, correlation = donor_compartment_correlation)
    fit_donor_compartment_ebayes <- limma::eBayes(fit_donor_compartment, robust = robust_ebayes)
    donor_compartment_group <- extract_limma_results(fit_donor_compartment_ebayes, coef_group, "Combined donor-compartment limma: PPI effect adjusted for organ") %>%
        dplyr::mutate(representation = representation, model = "donor_compartment_primary", contrast = "PPI_vs_control_adjusted_for_compartment",
            .before = 1)
    donor_compartment_organ <- extract_limma_results(fit_donor_compartment_ebayes, coef_compartment,
        "Combined donor-compartment limma: feces versus caecum") %>% dplyr::mutate(representation = representation,
        model = "donor_compartment_primary", contrast = "feces_vs_cecum_adjusted_for_group", .before = 1)
    design_all <- stats::model.matrix(~Separation + Organ, data = data.frame(Separation = metadata[[group_col]],
        Organ = metadata$Organ))
    block_all <- metadata[[donor_col]]
    corfit_all <- limma::duplicateCorrelation(object = Y_all, design = design_all, block = block_all)
    all_sample_correlation <- corfit_all$consensus.correlation
    if (is.null(all_sample_correlation)) {
        all_sample_correlation <- corfit_all$consensus
    }
    if (length(all_sample_correlation) != 1L || !is.finite(all_sample_correlation)) {
        stop(representation, ": all-sample duplicateCorrelation failed.")
    }
    fit_all <- limma::lmFit(object = Y_all, design = design_all, block = block_all, correlation = all_sample_correlation)
    fit_all_ebayes <- limma::eBayes(fit_all, robust = robust_ebayes)
    all_sample_group <- extract_limma_results(fit_all_ebayes, coef_group, "All-sample donor-blocked limma: PPI effect adjusted for organ") %>%
        dplyr::mutate(representation = representation, model = "all_samples_donor_blocked_sensitivity",
            contrast = "PPI_vs_control_adjusted_for_compartment", .before = 1)
    all_sample_organ <- extract_limma_results(fit_all_ebayes, coef_compartment, "All-sample donor-blocked limma: feces versus caecum") %>%
        dplyr::mutate(representation = representation, model = "all_samples_donor_blocked_sensitivity",
            contrast = "feces_vs_cecum_adjusted_for_group", .before = 1)
    design_interaction <- stats::model.matrix(~Separation * Organ, data = data.frame(Separation = metadata[[group_col]],
        Organ = metadata$Organ))
    interaction_coef <- paste0("Separation", group_levels[2], ":Organfeces")
    corfit_interaction <- limma::duplicateCorrelation(object = Y_all, design = design_interaction, block = block_all)
    interaction_correlation <- corfit_interaction$consensus.correlation
    if (is.null(interaction_correlation)) {
        interaction_correlation <- corfit_interaction$consensus
    }
    if (length(interaction_correlation) != 1L || !is.finite(interaction_correlation)) {
        stop(representation, ": interaction duplicateCorrelation failed.")
    }
    fit_interaction <- limma::lmFit(object = Y_all, design = design_interaction, block = block_all, correlation = interaction_correlation)
    contrast_matrix <- matrix(0, nrow = ncol(design_interaction), ncol = 3, dimnames = list(colnames(design_interaction),
        c("PPI_effect_in_cecum", "PPI_effect_in_feces", "PPI_effect_difference_feces_minus_cecum")))
    if (!all(c(coef_group, interaction_coef) %in% rownames(contrast_matrix))) {
        stop(representation, ": expected interaction coefficients were not found. Available: ", paste(rownames(contrast_matrix),
            collapse = ", "))
    }
    contrast_matrix[coef_group, "PPI_effect_in_cecum"] <- 1
    contrast_matrix[coef_group, "PPI_effect_in_feces"] <- 1
    contrast_matrix[interaction_coef, "PPI_effect_in_feces"] <- 1
    contrast_matrix[interaction_coef, "PPI_effect_difference_feces_minus_cecum"] <- 1
    fit_interaction_contrasts <- limma::contrasts.fit(fit_interaction, contrasts = contrast_matrix)
    fit_interaction_ebayes <- limma::eBayes(fit_interaction_contrasts, robust = robust_ebayes)
    interaction_results <- purrr::map_dfr(colnames(contrast_matrix), function(current_contrast) {
        extract_limma_results(fit_interaction_ebayes, current_contrast, paste0("All-sample interaction model: ",
            current_contrast)) %>% dplyr::mutate(representation = representation, model = "all_samples_group_by_compartment_interaction",
            contrast = current_contrast, .before = 1)
    })
    group_model_comparison <- donor_compartment_group %>% dplyr::select(representation, analyte, donor_compartment_log2_FC = log2_FC,
        donor_compartment_p_value = p_value, donor_compartment_FDR = FDR) %>% dplyr::full_join(all_sample_group %>%
        dplyr::select(representation, analyte, all_sample_log2_FC = log2_FC, all_sample_p_value = p_value,
            all_sample_FDR = FDR), by = c("representation", "analyte")) %>% dplyr::mutate(same_direction = sign(donor_compartment_log2_FC) ==
        sign(all_sample_log2_FC), absolute_log2FC_difference = abs(donor_compartment_log2_FC - all_sample_log2_FC)) %>%
        dplyr::arrange(donor_compartment_p_value)
    pseudocount_export <- tibble::tibble(analyte = names(prepared$pseudocounts), pseudocount = as.numeric(prepared$pseudocounts))
    donor_compartment_profile <- data.frame(analyte = rownames(Y_donor_compartment), Y_donor_compartment,
        check.names = FALSE)
    diagnostics <- tibble::tibble(representation = representation, item = c("Organ samples", "Recipient mice",
        "Independent donors", "Donor-compartment observations in primary model", "Analytes shared before combined filtering",
        "Analytes retained in combined model", "Analytes excluded by combined prevalence rule", "Cecum-only analytes excluded from combined model",
        "Feces-only analytes excluded from combined model", "Primary donor-compartment consensus correlation",
        "All-sample donor-blocked consensus correlation", "Interaction-model donor-blocked consensus correlation",
        "Primary combined PPI-effect FDR < 0.05", "Primary combined PPI-effect FDR < 0.10", "All-sample PPI-effect FDR < 0.05",
        "All-sample PPI-effect FDR < 0.10", "Interaction FDR < 0.10"), value = c(nrow(metadata), dplyr::n_distinct(metadata[[mouse_col]]),
        dplyr::n_distinct(metadata[[donor_col]]), nrow(donor_compartment_metadata), ncol(prepared$X_raw),
        ncol(prepared$X_log), nrow(prepared$excluded_unusable), length(prepared$only_cecum), length(prepared$only_feces),
        donor_compartment_correlation, all_sample_correlation, interaction_correlation, sum(donor_compartment_group$FDR <
            0.05, na.rm = TRUE), sum(donor_compartment_group$FDR < 0.1, na.rm = TRUE), sum(all_sample_group$FDR <
            0.05, na.rm = TRUE), sum(all_sample_group$FDR < 0.1, na.rm = TRUE), sum(interaction_results$contrast ==
            "PPI_effect_difference_feces_minus_cecum" & interaction_results$FDR < 0.1, na.rm = TRUE)))
    list(representation = representation, metadata = metadata, X_raw = prepared$X_raw, X_log = prepared$X_log,
        pseudocounts = prepared$pseudocounts, pseudocount_export = pseudocount_export, excluded_unusable = prepared$excluded_unusable,
        donor_compartment_metadata = donor_compartment_metadata, Y_donor_compartment = Y_donor_compartment,
        donor_compartment_profile = donor_compartment_profile, design_donor_compartment = design_donor_compartment,
        corfit_donor_compartment = corfit_donor_compartment, fit_donor_compartment = fit_donor_compartment,
        fit_donor_compartment_ebayes = fit_donor_compartment_ebayes, donor_compartment_group = donor_compartment_group,
        donor_compartment_organ = donor_compartment_organ, design_all = design_all, corfit_all = corfit_all,
        fit_all = fit_all, fit_all_ebayes = fit_all_ebayes, all_sample_group = all_sample_group, all_sample_organ = all_sample_organ,
        design_interaction = design_interaction, contrast_matrix = contrast_matrix, corfit_interaction = corfit_interaction,
        fit_interaction = fit_interaction, fit_interaction_ebayes = fit_interaction_ebayes, interaction_results = interaction_results,
        group_model_comparison = group_model_comparison, diagnostics = diagnostics)
}
if (isTRUE(run_combined_compartment_analysis)) {
    combined_compartment_results <- purrr::map(representations_to_analyze, ~run_combined_one_representation(clean_obj = clean_obj,
        representation = .x))
    names(combined_compartment_results) <- representations_to_analyze
    combined_primary_group <- purrr::map_dfr(combined_compartment_results, "donor_compartment_group")
    combined_primary_organ <- purrr::map_dfr(combined_compartment_results, "donor_compartment_organ")
    combined_all_sample_group <- purrr::map_dfr(combined_compartment_results, "all_sample_group")
    combined_all_sample_organ <- purrr::map_dfr(combined_compartment_results, "all_sample_organ")
    combined_interaction_results <- purrr::map_dfr(combined_compartment_results, "interaction_results")
    combined_group_model_comparison <- purrr::map_dfr(combined_compartment_results, "group_model_comparison")
    combined_compartment_diagnostics <- purrr::map_dfr(combined_compartment_results, "diagnostics")
    combined_workbook_sheets <- list(Diagnostics = combined_compartment_diagnostics, Primary_PPI_effect = combined_primary_group,
        Primary_compartment = combined_primary_organ, Allsample_PPI_effect = combined_all_sample_group,
        Allsample_compartment = combined_all_sample_organ, Group_by_compartment = combined_interaction_results,
        Primary_vs_allsample = combined_group_model_comparison)
    for (representation in names(combined_compartment_results)) {
        current <- combined_compartment_results[[representation]]
        rep_short <- ifelse(representation == "rt_resolved", "RT", "Met")
        combined_workbook_sheets[[paste0(rep_short, "_PrimaryGroup")]] <- current$donor_compartment_group
        combined_workbook_sheets[[paste0(rep_short, "_AllGroup")]] <- current$all_sample_group
        combined_workbook_sheets[[paste0(rep_short, "_Interaction")]] <- current$interaction_results
        combined_workbook_sheets[[paste0(rep_short, "_Compare")]] <- current$group_model_comparison
        combined_workbook_sheets[[paste0(rep_short, "_DonorCompProf")]] <- current$donor_compartment_profile
        combined_workbook_sheets[[paste0(rep_short, "_Pseudo")]] <- current$pseudocount_export
        combined_workbook_sheets[[paste0(rep_short, "_Excluded")]] <- current$excluded_unusable
    }
    combined_output_xlsx <- file.path(out_dir, "Cecum_Feces_LCMS_combined_compartment_covariate_results.xlsx")
    writexl::write_xlsx(combined_workbook_sheets, path = project_output(combined_output_xlsx))
    readr::write_csv(combined_primary_group, project_output(file.path(out_dir, "Cecum_Feces_LCMS_combined_primary_PPI_effect.csv"),
        name = "caecum_faeces_lcms_combined_primary_ppi_effect.csv"))
    readr::write_csv(combined_all_sample_group, project_output(file.path(out_dir, "Cecum_Feces_LCMS_combined_all_samples_PPI_effect.csv"),
        name = "caecum_faeces_lcms_combined_all_samples_ppi_effect.csv"))
    readr::write_csv(combined_all_sample_organ, project_output(file.path(out_dir, "Cecum_Feces_LCMS_combined_all_samples_compartment_effect.csv"),
        name = "caecum_faeces_lcms_combined_all_samples_compartment_effect.csv"))
    readr::write_csv(combined_interaction_results, project_output(file.path(out_dir, "Cecum_Feces_LCMS_combined_group_by_compartment_results.csv"),
        name = "caecum_faeces_lcms_combined_group_by_compartment_results.csv"))
    readr::write_csv(combined_group_model_comparison, project_output(file.path(out_dir, "Cecum_Feces_LCMS_combined_primary_vs_all_samples.csv"),
        name = "caecum_faeces_lcms_combined_primary_vs_all_samples.csv"))
    saveRDS(list(settings = list(representations_to_analyze = representations_to_analyze, group_levels = group_levels,
        combined_min_positive_per_group_compartment = combined_min_positive_per_group_compartment, combined_min_donors_detected_per_group_compartment = combined_min_donors_detected_per_group_compartment,
        robust_ebayes = robust_ebayes, interpretation = c("Primary combined model uses 12 donor x compartment means and blocks by donor.",
            "All-sample model uses all 32 organ samples and blocks by donor.", "The all-sample donor block is an equicorrelation approximation and does not separately estimate same-mouse pairing.",
            "The interaction tests whether the PPI effect differs between feces and caecum.")), combined_compartment_results = combined_compartment_results,
        combined_primary_group = combined_primary_group, combined_primary_organ = combined_primary_organ,
        combined_all_sample_group = combined_all_sample_group, combined_all_sample_organ = combined_all_sample_organ,
        combined_interaction_results = combined_interaction_results, combined_group_model_comparison = combined_group_model_comparison,
        combined_compartment_diagnostics = combined_compartment_diagnostics), project_output(file.path(out_dir,
        "Cecum_Feces_LCMS_combined_compartment_covariate_objects.rds"), name = "caecum_faeces_lcms_combined_compartment_covariate_objects.rds"))
    message("\nCombined caecum + feces compartment-adjusted analyses completed.", "\nPrimary model: donor x compartment means, blocked by donor.",
        "\nSensitivity: all 32 samples, blocked by donor.", "\nInteraction: PPI effect difference between feces and caecum.",
        "\n\nCombined PPI-effect FDR < 0.10 by representation:\n", paste(combined_primary_group %>% dplyr::group_by(representation) %>%
            dplyr::summarise(n = sum(FDR < 0.1, na.rm = TRUE), .groups = "drop") %>% dplyr::transmute(text = paste0(representation,
            ": ", n)) %>% dplyr::pull(text), collapse = "\n"), "\n\nCombined output:\n", combined_output_xlsx)
}
}))
