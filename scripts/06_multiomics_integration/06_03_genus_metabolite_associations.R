# SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0
# Genus-metabolite associations and donor exclusion
# Run through run.R so input and output paths resolve within the package.
# Inputs: see data/input_manifest.csv and docs/RUN_ORDER.md.
# Writes generated files under results/; bundled inputs are read-only.

message("Loaded 06_03_genus_metabolite_associations.R")
suppressPackageStartupMessages({
    library(dplyr)
    library(tidyr)
    library(purrr)
    library(readr)
    library(stringr)
    library(tibble)
    library(limma)
    library(writexl)
})
script_path <- tryCatch(ppi_normalize_path(sys.frame(1)$ofile, winslash = "/", mustWork = FALSE), error = function(e) NA_character_)
script_dir <- if (!is.na(script_path) && nzchar(script_path)) dirname(script_path) else getwd()
config_path <- ppi_code("scripts/06_multiomics_integration/06_01_configure_multiomics.R")
if (!ppi_file_exists(config_path)) config_path <- "06_01_configure_multiomics.R"
ppi_source(config_path)
prep_file <- file.path(association_root, "01_prepared_inputs", "PPI_multiomics_prepared_inputs.rds")
if (!ppi_file_exists(prep_file)) {
    stop("Prepared input bundle not found: ", prep_file, "\nRun 06_02_prepare_multiomics_inputs.R first.")
}
bundle <- readRDS(ppi_input(prep_file))
settings <- bundle$settings
results_dir <- file.path(association_root, "02_association_results")
family_dir <- file.path(results_dir, "complete_family_tables")
dir.create(ppi_output(family_dir, directory = TRUE), recursive = TRUE, showWarnings = FALSE)
rank_z_vector <- function(x) {
    x <- as.numeric(x)
    r <- rank(x, ties.method = "average", na.last = "keep")
    z <- as.numeric(scale(r))
    z[!is.finite(z)] <- 0
    z
}
rank_z_matrix <- function(X) {
    X <- as.matrix(X)
    out <- vapply(seq_len(ncol(X)), function(j) rank_z_vector(X[, j]), numeric(nrow(X)))
    out <- as.matrix(out)
    rownames(out) <- rownames(X)
    colnames(out) <- colnames(X)
    out
}
collapse_by_key <- function(X, metadata, key_col) {
    X <- as.matrix(X)
    metadata <- as.data.frame(metadata)
    keys <- as.character(metadata[[key_col]])
    valid <- !is.na(keys) & keys != ""
    X <- X[valid, , drop = FALSE]
    metadata <- metadata[valid, , drop = FALSE]
    keys <- keys[valid]
    if (anyDuplicated(keys)) {
        X_sum <- rowsum(X, group = keys, reorder = FALSE)
        key_counts <- as.numeric(table(keys)[rownames(X_sum)])
        X <- sweep(X_sum, 1, key_counts, "/")
        key_order <- rownames(X)
        metadata <- metadata %>% filter(.data[[key_col]] %in% key_order) %>% arrange(match(.data[[key_col]], 
            key_order)) %>% distinct(.data[[key_col]], .keep_all = TRUE)
    }
    else {
        rownames(X) <- keys
        metadata <- metadata[match(rownames(X), metadata[[key_col]]), , drop = FALSE]
    }
    rownames(metadata) <- as.character(metadata[[key_col]])
    list(matrix = X, metadata = metadata)
}
align_predictor_outcome <- function(predictor_dataset, outcome_dataset, key_col, minimum_samples, minimum_donors, 
    require_day = FALSE) {
    pred <- collapse_by_key(predictor_dataset$matrix, predictor_dataset$metadata, key_col)
    out <- collapse_by_key(outcome_dataset$matrix, outcome_dataset$metadata, key_col)
    common <- intersect(rownames(pred$matrix), rownames(out$matrix))
    common <- sort(common)
    if (length(common) < minimum_samples) {
        return(list(valid = FALSE, reason = paste0("Only ", length(common), " matched samples; minimum is ", 
            minimum_samples), n_samples = length(common)))
    }
    X <- pred$matrix[common, , drop = FALSE]
    Y <- out$matrix[common, , drop = FALSE]
    pred_meta <- pred$metadata[common, , drop = FALSE]
    out_meta <- out$metadata[common, , drop = FALSE]
    pred_group <- as.character(pred_meta$Group)
    out_group <- as.character(out_meta$Group)
    pred_donor <- as.character(pred_meta$Donor)
    out_donor <- as.character(out_meta$Donor)
    inconsistent_group <- !is.na(pred_group) & !is.na(out_group) & pred_group != out_group
    inconsistent_donor <- !is.na(pred_donor) & !is.na(out_donor) & pred_donor != out_donor
    if (any(inconsistent_group)) {
        stop("PPI-group mismatch among matched samples: ", paste(common[inconsistent_group], collapse = ", "))
    }
    if (any(inconsistent_donor)) {
        stop("Donor mismatch among matched samples: ", paste(common[inconsistent_donor], collapse = ", "))
    }
    meta <- tibble(Sample_Key = common, Group = factor(coalesce(pred_group, out_group), levels = c("ctrl", 
        "ppi")), Donor = factor(coalesce(pred_donor, out_donor)))
    if (anyNA(meta$Group)) {
        stop("Missing PPI group after matching: ", paste(meta$Sample_Key[is.na(meta$Group)], collapse = ", "))
    }
    if (anyNA(meta$Donor)) {
        stop("Missing donor after matching: ", paste(meta$Sample_Key[is.na(meta$Donor)], collapse = ", "))
    }
    if (isTRUE(require_day)) {
        pred_day <- if ("Day" %in% names(pred_meta)) {
            suppressWarnings(as.integer(as.character(pred_meta$Day)))
        }
        else {
            rep(NA_integer_, length(common))
        }
        out_day <- if ("Day" %in% names(out_meta)) {
            suppressWarnings(as.integer(as.character(out_meta$Day)))
        }
        else {
            rep(NA_integer_, length(common))
        }
        inconsistent_day <- !is.na(pred_day) & !is.na(out_day) & pred_day != out_day
        if (any(inconsistent_day)) {
            stop("Day mismatch among matched samples: ", paste(common[inconsistent_day], collapse = ", "))
        }
        meta$Day <- dplyr::coalesce(pred_day, out_day)
        meta$DayF <- factor(meta$Day, levels = settings$exvivo_days)
        if (anyNA(meta$DayF)) {
            stop("Missing or out-of-range ex vivo day after matching: ", paste(meta$Sample_Key[is.na(meta$DayF)], 
                collapse = ", "))
        }
    }
    n_donors <- n_distinct(meta$Donor)
    if (n_donors < minimum_donors) {
        return(list(valid = FALSE, reason = paste0("Only ", n_donors, " matched donors; minimum is ", 
            minimum_donors), n_samples = length(common), n_donors = n_donors))
    }
    list(valid = TRUE, X = X, Y = Y, metadata = meta, n_samples = nrow(X), n_donors = n_donors, matched_keys = common)
}
estimate_consensus_correlation <- function(Y, design, block) {
    if (length(unique(block)) == length(block)) {
        return(list(correlation = 0, status = "no_repeated_blocks"))
    }
    result <- tryCatch(suppressWarnings(limma::duplicateCorrelation(t(Y), design = design, block = block)), 
        error = function(e) e)
    if (inherits(result, "error") || !is.finite(result$consensus.correlation)) {
        status <- if (inherits(result, "error")) 
            conditionMessage(result)
        else "non_finite_consensus"
        return(list(correlation = 0, status = paste0("fallback_zero: ", status)))
    }
    correlation <- max(min(as.numeric(result$consensus.correlation), 0.95), -0.95)
    list(correlation = correlation, status = "estimated")
}
raw_spearman_cross_matrix <- function(X, Y) {
    stats::cor(X, Y, method = "spearman", use = "pairwise.complete.obs")
}
run_association_family <- function(aligned, family_id, layer, predictor_dataset_id, outcome_dataset_id, 
    predictor_type, outcome_type, organ = NA_character_, library = NA_character_, include_day = FALSE, 
    predictor_annotation = NULL, outcome_annotation = NULL) {
    if (!isTRUE(aligned$valid)) {
        return(list(results = tibble(), diagnostics = tibble(family_id = family_id, layer = layer, predictor_dataset = predictor_dataset_id, 
            outcome_dataset = outcome_dataset_id, status = "skipped", message = aligned$reason, n_samples = aligned$n_samples %||% 
                NA_integer_, n_donors = aligned$n_donors %||% NA_integer_), model_data = NULL))
    }
    X_original <- aligned$X
    Y_original <- aligned$Y
    meta <- aligned$metadata
    keep_x <- vapply(seq_len(ncol(X_original)), function(j) {
        is.finite(var(X_original[, j], na.rm = TRUE)) && var(X_original[, j], na.rm = TRUE) > 0
    }, logical(1))
    keep_y <- vapply(seq_len(ncol(Y_original)), function(j) {
        is.finite(var(Y_original[, j], na.rm = TRUE)) && var(Y_original[, j], na.rm = TRUE) > 0
    }, logical(1))
    X_original <- X_original[, keep_x, drop = FALSE]
    Y_original <- Y_original[, keep_y, drop = FALSE]
    if (ncol(X_original) == 0L || ncol(Y_original) == 0L) {
        return(list(results = tibble(), diagnostics = tibble(family_id = family_id, layer = layer, predictor_dataset = predictor_dataset_id, 
            outcome_dataset = outcome_dataset_id, status = "skipped", message = "No non-invariant predictors or outcomes after matching", 
            n_samples = nrow(X_original), n_donors = n_distinct(meta$Donor)), model_data = NULL))
    }
    X_model <- if (isTRUE(settings$rank_standardize)) 
        rank_z_matrix(X_original)
    else scale(X_original)
    Y_model <- if (isTRUE(settings$rank_standardize)) 
        rank_z_matrix(Y_original)
    else scale(Y_original)
    X_model <- as.matrix(X_model)
    Y_model <- as.matrix(Y_model)
    raw_rho <- raw_spearman_cross_matrix(X_original, Y_original)
    null_formula <- if (include_day) 
        ~Group + DayF
    else ~Group
    null_design <- model.matrix(null_formula, data = meta)
    cor_info <- estimate_consensus_correlation(Y = Y_model, design = null_design, block = meta$Donor)
    result_list <- vector("list", ncol(X_model))
    for (j in seq_len(ncol(X_model))) {
        dat <- meta
        dat$Predictor <- X_model[, j]
        full_formula <- if (include_day) 
            ~Group + DayF + Predictor
        else ~Group + Predictor
        design <- model.matrix(full_formula, data = dat)
        coef_name <- "Predictor"
        fit <- limma::lmFit(t(Y_model), design = design, block = dat$Donor, correlation = cor_info$correlation)
        fit <- limma::eBayes(fit, robust = settings$robust_ebayes, trend = settings$trend_ebayes)
        coef_index <- match(coef_name, colnames(design))
        beta <- fit$coefficients[, coef_index]
        t_stat <- fit$t[, coef_index]
        p_value <- fit$p.value[, coef_index]
        moderated_se <- fit$stdev.unscaled[, coef_index] * sqrt(fit$s2.post)
        df_total <- fit$df.total
        critical <- qt(0.975, df = df_total)
        outcome_names <- rownames(fit$coefficients)
        predictor_name <- colnames(X_model)[[j]]
        result_list[[j]] <- tibble(family_id = family_id, layer = layer, predictor_dataset = predictor_dataset_id, 
            outcome_dataset = outcome_dataset_id, organ = organ, library = library, predictor_type = predictor_type, 
            outcome_type = outcome_type, predictor = predictor_name, outcome = outcome_names, raw_spearman_rho = as.numeric(raw_rho[predictor_name, 
                outcome_names]), adjusted_rank_beta = as.numeric(beta), moderated_SE = as.numeric(moderated_se), 
            conf_low = as.numeric(beta - critical * moderated_se), conf_high = as.numeric(beta + critical * 
                moderated_se), t_statistic = as.numeric(t_stat), p_value = as.numeric(p_value), df_total = as.numeric(df_total), 
            n_samples = aligned$n_samples, n_donors = aligned$n_donors, n_control_samples = sum(meta$Group == 
                "ctrl"), n_ppi_samples = sum(meta$Group == "ppi"), consensus_donor_correlation = cor_info$correlation, 
            consensus_status = cor_info$status, adjusted_for_group = TRUE, adjusted_for_day = include_day, 
            repeated_measure_block = "Donor")
    }
    results <- bind_rows(result_list) %>% mutate(FDR = p.adjust(p_value, method = "BH"), significant_FDR0p1 = is.finite(FDR) & 
        FDR < settings$association_fdr_cutoff, effect_direction = case_when(adjusted_rank_beta > 0 ~ 
        "positive", adjusted_rank_beta < 0 ~ "negative", TRUE ~ "zero")) %>% arrange(FDR, p_value)
    if (!is.null(predictor_annotation) && nrow(predictor_annotation) > 0L) {
        results <- results %>% left_join(predictor_annotation, by = c(predictor = "variable"))
    }
    if (!is.null(outcome_annotation) && nrow(outcome_annotation) > 0L) {
        results <- results %>% left_join(outcome_annotation, by = c(outcome = "variable"), suffix = c("_predictor", 
            "_outcome"))
    }
    diagnostics <- tibble(family_id = family_id, layer = layer, predictor_dataset = predictor_dataset_id, 
        outcome_dataset = outcome_dataset_id, organ = organ, library = library, status = "completed", 
        message = NA_character_, n_samples = aligned$n_samples, n_donors = aligned$n_donors, n_control_samples = sum(meta$Group == 
            "ctrl"), n_ppi_samples = sum(meta$Group == "ppi"), n_predictors = ncol(X_model), n_outcomes = ncol(Y_model), 
        n_tests = nrow(results), n_FDR0p1 = sum(results$significant_FDR0p1, na.rm = TRUE), consensus_donor_correlation = cor_info$correlation, 
        consensus_status = cor_info$status)
    model_data <- list(family_id = family_id, X_model = X_model, Y_model = Y_model, X_original = X_original, 
        Y_original = Y_original, metadata = meta, include_day = include_day, consensus_correlation = cor_info$correlation)
    list(results = results, diagnostics = diagnostics, model_data = model_data)
}
`%||%` <- function(x, y) {
    if (is.null(x) || length(x) == 0L) 
        y
    else x
}
safe_family_filename <- function(family_id) {
    paste0(str_replace_all(family_id, "[^A-Za-z0-9._-]+", "_"), ".csv.gz")
}
save_family_result <- function(family_result) {
    if (nrow(family_result$results) == 0L) 
        return(invisible(NULL))
    write_csv(family_result$results, ppi_output(file.path(family_dir, safe_family_filename(unique(family_result$results$family_id)))))
}
all_family_results <- list()
all_family_diagnostics <- list()
family_model_data <- list()
register_family <- function(result) {
    id <- if (nrow(result$diagnostics) > 0L) 
        result$diagnostics$family_id[[1]]
    else paste0("family_", length(all_family_results) + 1L)
    all_family_results[[id]] <<- result$results
    all_family_diagnostics[[id]] <<- result$diagnostics
    family_model_data[[id]] <<- result$model_data
    save_family_result(result)
}
if (isTRUE(bundle$run_layers$metabolite_pathway)) {
    for (met_id in bundle$mouse_metabolomics_datasets) {
        met_ds <- bundle$metabolomics[[met_id]]
        for (path_id in names(bundle$pathways)) {
            path_ds <- bundle$pathways[[path_id]]
            family_id <- paste("metabolite_pathway", met_id, path_ds$organ, path_ds$library, sep = "__")
            message("Running ", family_id)
            aligned <- align_predictor_outcome(predictor_dataset = met_ds, outcome_dataset = path_ds, 
                key_col = "Mouse_Key", minimum_samples = settings$minimum_samples_mouse, minimum_donors = settings$minimum_donors, 
                require_day = FALSE)
            pathway_annotation <- path_ds$selection %>% transmute(variable = Pathway, pathway_GSEA_NES = GSEA_NES, 
                pathway_GSEA_P = GSEA_P, pathway_GSEA_FDR = GSEA_FDR)
            register_family(run_association_family(aligned = aligned, family_id = family_id, layer = "metabolite_pathway", 
                predictor_dataset_id = met_id, outcome_dataset_id = path_id, predictor_type = "metabolite", 
                outcome_type = "pathway", organ = path_ds$organ, library = path_ds$library, include_day = FALSE, 
                outcome_annotation = pathway_annotation))
        }
    }
}
if (isTRUE(bundle$run_layers$metabolite_genus)) {
    for (met_id in bundle$mouse_metabolomics_datasets) {
        family_id <- paste("metabolite_genus", met_id, "mouse_genus", sep = "__")
        message("Running ", family_id)
        aligned <- align_predictor_outcome(predictor_dataset = bundle$metabolomics[[met_id]], outcome_dataset = bundle$microbiome$mouse_genus, 
            key_col = "Mouse_Key", minimum_samples = settings$minimum_samples_mouse, minimum_donors = settings$minimum_donors, 
            require_day = FALSE)
        register_family(run_association_family(aligned = aligned, family_id = family_id, layer = "metabolite_genus", 
            predictor_dataset_id = met_id, outcome_dataset_id = "mouse_genus", predictor_type = "metabolite", 
            outcome_type = "genus", include_day = FALSE))
    }
    for (met_id in bundle$exvivo_metabolomics_datasets) {
        family_id <- paste("metabolite_genus", met_id, "culture_genus", sep = "__")
        message("Running ", family_id)
        aligned <- align_predictor_outcome(predictor_dataset = bundle$metabolomics[[met_id]], outcome_dataset = bundle$microbiome$culture_genus, 
            key_col = "Culture_Key", minimum_samples = settings$minimum_samples_exvivo, minimum_donors = settings$minimum_donors, 
            require_day = TRUE)
        register_family(run_association_family(aligned = aligned, family_id = family_id, layer = "metabolite_genus", 
            predictor_dataset_id = met_id, outcome_dataset_id = "culture_genus", predictor_type = "metabolite", 
            outcome_type = "genus", include_day = TRUE))
    }
}
if (isTRUE(bundle$run_layers$genus_pathway)) {
    for (path_id in names(bundle$pathways)) {
        path_ds <- bundle$pathways[[path_id]]
        family_id <- paste("genus_pathway", "mouse_genus", path_ds$organ, path_ds$library, sep = "__")
        message("Running ", family_id)
        aligned <- align_predictor_outcome(predictor_dataset = bundle$microbiome$mouse_genus, outcome_dataset = path_ds, 
            key_col = "Mouse_Key", minimum_samples = settings$minimum_samples_mouse, minimum_donors = settings$minimum_donors, 
            require_day = FALSE)
        pathway_annotation <- path_ds$selection %>% transmute(variable = Pathway, pathway_GSEA_NES = GSEA_NES, 
            pathway_GSEA_P = GSEA_P, pathway_GSEA_FDR = GSEA_FDR)
        register_family(run_association_family(aligned = aligned, family_id = family_id, layer = "genus_pathway", 
            predictor_dataset_id = "mouse_genus", outcome_dataset_id = path_id, predictor_type = "genus", 
            outcome_type = "pathway", organ = path_ds$organ, library = path_ds$library, include_day = FALSE, 
            outcome_annotation = pathway_annotation))
    }
}
resolve_matrix_column <- function(target, available_names) {
    target <- as.character(target)[1]
    available_names <- as.character(available_names)
    if (length(target) == 0L || is.na(target) || target == "" || length(available_names) == 0L) {
        return(NA_integer_)
    }
    idx <- match(target, available_names)
    if (!is.na(idx)) 
        return(idx)
    target_trim <- stringr::str_squish(target)
    available_trim <- stringr::str_squish(available_names)
    idx <- match(target_trim, available_trim)
    if (!is.na(idx)) 
        return(idx)
    target_ascii <- iconv(target_trim, from = "", to = "ASCII//TRANSLIT")
    available_ascii <- iconv(available_trim, from = "", to = "ASCII//TRANSLIT")
    idx <- match(target_ascii, available_ascii)
    if (!is.na(idx)) 
        return(idx)
    target_safe <- make.names(target_ascii, unique = FALSE)
    available_safe <- make.names(available_ascii, unique = FALSE)
    candidates <- which(available_safe == target_safe)
    if (length(candidates) == 1L) 
        return(candidates[[1]])
    NA_integer_
}
loo_direction_stability <- function(model_data, significant_results) {
    empty_result <- tibble(family_id = character(), predictor = character(), outcome = character(), full_adjusted_rank_beta = numeric(), 
        n_donors_tested = integer(), n_successful_refits = integer(), proportion_same_direction = numeric(), 
        minimum_LOO_beta = numeric(), maximum_LOO_beta = numeric(), LOO_coefficients = character(), LOO_status = character())
    if (is.null(model_data) || nrow(significant_results) == 0L) {
        return(empty_result)
    }
    X_model <- as.matrix(model_data$X_model)
    Y_model <- as.matrix(model_data$Y_model)
    meta <- as.data.frame(model_data$metadata)
    if (is.null(colnames(X_model)) && !is.null(colnames(model_data$X_original))) {
        colnames(X_model) <- colnames(model_data$X_original)
    }
    if (is.null(colnames(Y_model)) && !is.null(colnames(model_data$Y_original))) {
        colnames(Y_model) <- colnames(model_data$Y_original)
    }
    if (is.null(colnames(X_model)) || is.null(colnames(Y_model))) {
        stop("LOO refitting cannot proceed for family ", model_data$family_id, ": X_model or Y_model has no column names.")
    }
    if (nrow(X_model) != nrow(meta) || nrow(Y_model) != nrow(meta)) {
        stop("LOO refitting cannot proceed for family ", model_data$family_id, ": model matrices and metadata have inconsistent row counts.")
    }
    donors <- sort(unique(as.character(meta$Donor)))
    donors <- donors[!is.na(donors) & donors != ""]
    if (length(donors) < 2L) {
        return(significant_results %>% transmute(family_id = as.character(.data$family_id), predictor = as.character(.data$predictor), 
            outcome = as.character(.data$outcome), full_adjusted_rank_beta = as.numeric(.data$adjusted_rank_beta), 
            n_donors_tested = length(donors), n_successful_refits = 0L, proportion_same_direction = NA_real_, 
            minimum_LOO_beta = NA_real_, maximum_LOO_beta = NA_real_, LOO_coefficients = NA_character_, 
            LOO_status = "fewer_than_two_donors"))
    }
    donor_cache <- setNames(vector("list", length(donors)), donors)
    for (donor in donors) {
        keep <- as.character(meta$Donor) != donor
        dat <- droplevels(meta[keep, , drop = FALSE])
        null_formula <- if (isTRUE(model_data$include_day)) {
            ~Group + DayF
        }
        else {
            ~Group
        }
        null_design <- model.matrix(null_formula, data = dat)
        cor_result <- tryCatch(estimate_consensus_correlation(Y = Y_model[keep, , drop = FALSE], design = null_design, 
            block = dat$Donor), error = function(e) list(correlation = model_data$consensus_correlation))
        correlation <- as.numeric(cor_result$correlation)[1]
        if (!is.finite(correlation)) 
            correlation <- 0
        correlation <- max(min(correlation, 0.99), -0.99)
        donor_cache[[donor]] <- list(keep = keep, metadata = dat, correlation = correlation)
    }
    out <- vector("list", nrow(significant_results))
    for (i in seq_len(nrow(significant_results))) {
        hit <- significant_results[i, , drop = FALSE]
        predictor <- as.character(hit$predictor[[1]])
        outcome <- as.character(hit$outcome[[1]])
        full_beta <- as.numeric(hit$adjusted_rank_beta[[1]])
        predictor_index <- resolve_matrix_column(predictor, colnames(X_model))
        outcome_index <- resolve_matrix_column(outcome, colnames(Y_model))
        if (is.na(predictor_index) || is.na(outcome_index)) {
            missing_part <- paste(c(if (is.na(predictor_index)) paste0("predictor_not_found:", predictor), 
                if (is.na(outcome_index)) paste0("outcome_not_found:", outcome)), collapse = " | ")
            out[[i]] <- tibble(family_id = as.character(hit$family_id[[1]]), predictor = predictor, outcome = outcome, 
                full_adjusted_rank_beta = full_beta, n_donors_tested = length(donors), n_successful_refits = 0L, 
                proportion_same_direction = NA_real_, minimum_LOO_beta = NA_real_, maximum_LOO_beta = NA_real_, 
                LOO_coefficients = NA_character_, LOO_status = missing_part)
            next
        }
        coefficients <- rep(NA_real_, length(donors))
        names(coefficients) <- donors
        for (d in seq_along(donors)) {
            donor <- donors[[d]]
            cached <- donor_cache[[donor]]
            keep <- cached$keep
            dat <- cached$metadata
            dat$Predictor <- as.numeric(X_model[keep, predictor_index])
            full_formula <- if (isTRUE(model_data$include_day)) {
                ~Group + DayF + Predictor
            }
            else {
                ~Group + Predictor
            }
            design <- model.matrix(full_formula, data = dat)
            predictor_coef_index <- match("Predictor", colnames(design))
            if (is.na(predictor_coef_index)) 
                next
            y_reduced <- Y_model[keep, outcome_index, drop = FALSE]
            fit <- tryCatch(suppressWarnings(limma::lmFit(t(y_reduced), design = design, block = dat$Donor, 
                correlation = cached$correlation)), error = function(e) NULL)
            if (!is.null(fit)) {
                coefficients[[d]] <- as.numeric(fit$coefficients[1, predictor_coef_index])
            }
        }
        finite_coef <- coefficients[is.finite(coefficients)]
        same_direction <- if (length(finite_coef) > 0L && is.finite(full_beta) && full_beta != 0) {
            mean(sign(finite_coef) == sign(full_beta))
        }
        else {
            NA_real_
        }
        out[[i]] <- tibble(family_id = as.character(hit$family_id[[1]]), predictor = predictor, outcome = outcome, 
            full_adjusted_rank_beta = full_beta, n_donors_tested = length(donors), n_successful_refits = length(finite_coef), 
            proportion_same_direction = same_direction, minimum_LOO_beta = if (length(finite_coef) > 
                0L) 
                min(finite_coef)
            else NA_real_, maximum_LOO_beta = if (length(finite_coef) > 0L) 
                max(finite_coef)
            else NA_real_, LOO_coefficients = paste(paste0(names(coefficients), "=", signif(coefficients, 
                4)), collapse = ";"), LOO_status = if (length(finite_coef) == length(donors)) {
                "all_refits_successful"
            }
            else if (length(finite_coef) > 0L) {
                "partial_refit_failure"
            }
            else {
                "all_refits_failed"
            })
    }
    bind_rows(out)
}
all_results <- bind_rows(all_family_results)
diagnostics <- bind_rows(all_family_diagnostics)
if (nrow(all_results) == 0L) {
    stop("No association family produced results. Inspect association_family_diagnostics.csv after rerunning.")
}
significant_results <- all_results %>% filter(.data$significant_FDR0p1) %>% arrange(layer, family_id, 
    FDR)
loo_results <- tibble(family_id = character(), predictor = character(), outcome = character(), full_adjusted_rank_beta = numeric(), 
    n_donors_tested = integer(), n_successful_refits = integer(), proportion_same_direction = numeric(), 
    minimum_LOO_beta = numeric(), maximum_LOO_beta = numeric(), LOO_coefficients = character(), LOO_status = character())
if (isTRUE(settings$run_leave_one_donor_out) && nrow(significant_results) > 0L) {
    message("Running leave-one-donor-out direction checks for FDR-significant associations...")
    loo_results <- imap_dfr(family_model_data, function(model_data, family_id) {
        hits <- significant_results %>% filter(.data$family_id == .env$family_id)
        if (nrow(hits) == 0L || is.null(model_data)) {
            return(tibble())
        }
        loo_direction_stability(model_data, hits)
    })
    loo_key_counts <- loo_results %>% count(.data$family_id, .data$predictor, .data$outcome, name = "n_key_rows")
    duplicated_loo_keys <- loo_key_counts %>% filter(.data$n_key_rows != 1L)
    if (nrow(duplicated_loo_keys) > 0L) {
        stop("LOO results contain duplicated or missing association keys. ", "This indicates that hits were not restricted to their own family.")
    }
    if (nrow(loo_results) != nrow(significant_results)) {
        stop("LOO row count (", nrow(loo_results), ") does not equal the number of significant associations (", 
            nrow(significant_results), ").")
    }
    significant_results <- significant_results %>% left_join(loo_results, by = c("family_id", "predictor", 
        "outcome"), relationship = "one-to-one")
}
if (any(diagnostics$status != "completed")) {
    warning(sum(diagnostics$status != "completed"), " association families were skipped. Review association_family_diagnostics.csv.")
}
write_csv(diagnostics, ppi_output(file.path(results_dir, "association_family_diagnostics.csv"), name = "association_family_diagnostics.csv"))
write_csv(significant_results, ppi_output(file.path(results_dir, "ALL_significant_associations_FDR0p1.csv.gz")))
write_csv(loo_results, ppi_output(file.path(results_dir, "leave_one_donor_out_direction_stability.csv"), 
    name = "leave_one_donor_out_direction_stability.csv"))
layer_summary <- diagnostics %>% group_by(layer) %>% summarise(n_families = n(), n_completed = sum(status == 
    "completed"), n_tests = sum(n_tests, na.rm = TRUE), n_FDR0p1 = sum(n_FDR0p1, na.rm = TRUE), .groups = "drop")
family_summary <- diagnostics %>% arrange(layer, family_id)
write_csv(layer_summary, ppi_output(file.path(results_dir, "association_layer_summary.csv"), name = "association_layer_summary.csv"))
write_csv(family_summary, ppi_output(file.path(results_dir, "association_family_summary.csv"), name = "association_family_summary.csv"))
writexl::write_xlsx(list(Layer_summary = layer_summary, Family_summary = family_summary, Significant_all = significant_results, 
    LOO_direction_stability = loo_results), path = ppi_output(file.path(results_dir, "PPI_multiomics_association_summary.xlsx"), 
    name = "ppi_multiomics_association_summary.xlsx"))
saveRDS(list(all_results = all_results, significant_results = significant_results, diagnostics = diagnostics, 
    loo_results = loo_results, family_model_data = family_model_data, settings = settings), ppi_output(file.path(results_dir, 
    "PPI_multiomics_association_results.rds"), name = "ppi_multiomics_association_results.rds"), compress = "xz")
message("Association analysis complete. Outputs written to: ", results_dir)
