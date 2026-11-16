# SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0
# DESeq2 covariate sensitivity
# Run through run.R so input and output paths resolve within the package.
# Inputs: see data/input_manifest.csv and docs/RUN_ORDER.md.
# Writes generated files under results/; bundled inputs are read-only.

suppressPackageStartupMessages({
    library(DESeq2)
    library(SummarizedExperiment)
    library(dplyr)
    library(stringr)
    library(readxl)
    library(tibble)
    library(writexl)
})
load_rnaseq_script <- file.path(ppi_code("scripts/05_mouse_transcriptomics/05_02_load_primary_deseq2_cache.R"))
ppi_source(load_rnaseq_script)
if (!exists("dds_list")) stop("dds_list not found.")
if (!exists("meta")) stop("meta not found.")
out_dir <- file.path(getwd(), "RNAseq_DESeq2_covariate_sensitivity")
dir.create(ppi_output(out_dir, directory = TRUE), recursive = TRUE, showWarnings = FALSE)
exclude_samples <- c("25Jul1866-C02", "25Jul1866-G04")
alpha <- 0.1
mouse_table_path <- ppi_path("Validation/PPI_mouse_table.xlsx")
mouse_cov <- readxl::read_xlsx(ppi_input(mouse_table_path)) %>% dplyr::rename(mouse_number = Mouse, isolator = Isolator, 
    cage_position = Käfigplatz, line = Linie, sex_male = `♂`, dob = DOB, source_zk = `aus ZK`, age_weeks = `Alter (w)`, 
    weaned = Abgesetzt, inoculum = Inoculum, cage = Käfig) %>% dplyr::mutate(mouse_number = as.integer(mouse_number), 
    isolator = factor(isolator), cage_position = factor(cage_position), source_zk = factor(source_zk), 
    inoculum = factor(inoculum), cage = factor(cage), age_weeks = as.numeric(age_weeks))
meta_mouse <- meta %>% as.data.frame() %>% tibble::rownames_to_column("sampleID")
external_name_col <- c("external_name", "ExternalName", "externalName", "External_name", "external.name", 
    "External.Name", "External name", "external name")
external_name_col <- external_name_col[external_name_col %in% colnames(meta_mouse)][1]
if (is.na(external_name_col)) {
    stop("Could not find external name column in meta. Available columns are:\n", paste(colnames(meta_mouse), 
        collapse = ", "))
}
meta_mouse <- meta_mouse %>% dplyr::mutate(external_name = as.character(.data[[external_name_col]]), 
    mouse_number = as.integer(stringr::str_match(external_name, "_m_(\\d+)_")[, 2]), donor_id = stringr::str_extract(external_name, 
        "^(Ctrl|PPI)_\\d+"), group = factor(stringr::str_extract(external_name, "^(Ctrl|PPI)"), levels = c("Ctrl", 
        "PPI")), organ_from_external_name = stringr::str_extract(external_name, "(Brain|Colon|Liver)$")) %>% 
    dplyr::left_join(mouse_cov, by = "mouse_number")
message("Checking mouse metadata matching...")
failed_mouse_matches <- meta_mouse %>% dplyr::filter(is.na(age_weeks) | is.na(isolator)) %>% dplyr::select(sampleID, 
    name, external_name, mouse_number, donor_id, group)
if (nrow(failed_mouse_matches) > 0) {
    warning("Some samples did not match the mouse covariate table.")
    print(failed_mouse_matches)
}
donor_mismatches <- meta_mouse %>% dplyr::filter(!is.na(inoculum), !is.na(donor_id), as.character(inoculum) != 
    as.character(donor_id)) %>% dplyr::select(sampleID, name, external_name, donor_id, inoculum, mouse_number)
if (nrow(donor_mismatches) > 0) {
    warning("Some donor IDs differ between external name and mouse table.")
    print(donor_mismatches)
}
write.csv(meta_mouse, ppi_output(file.path(out_dir, "meta_with_mouse_covariates.csv"), name = "meta_with_mouse_covariates.csv"), 
    row.names = FALSE)
add_mouse_covariates_to_dds <- function(dds, meta_mouse) {
    m <- match(colnames(dds), meta_mouse$name)
    if (any(is.na(m))) {
        stop("Some dds samples could not be matched to meta_mouse$name:\n", paste(colnames(dds)[is.na(m)], 
            collapse = ", "))
    }
    colData(dds)$rnaseq_sample_id <- meta_mouse$name[m]
    colData(dds)$external_name <- meta_mouse$external_name[m]
    colData(dds)$mouse_number <- meta_mouse$mouse_number[m]
    colData(dds)$donor_id <- factor(meta_mouse$donor_id[m])
    colData(dds)$age_weeks <- meta_mouse$age_weeks[m]
    colData(dds)$isolator <- factor(meta_mouse$isolator[m])
    colData(dds)$cage_position <- factor(meta_mouse$cage_position[m])
    colData(dds)$cage <- factor(meta_mouse$cage[m])
    colData(dds)$group <- factor(meta_mouse$group[m], levels = c("Ctrl", "PPI"))
    dds
}
run_deseq2_sensitivity_one <- function(dds, organ_name, model_label, design_formula, contrast = c("group", 
    "PPI", "Ctrl"), result_name = NULL, add_mouse_covariates = TRUE, meta_mouse, exclude_samples, alpha = 0.1, 
    min_count = 10, min_prop_samples = 0.2) {
    message("------------------------------------------------------------")
    message("Running ", organ_name, " | ", model_label)
    message("Design: ", paste(deparse(design_formula), collapse = " "))
    message("------------------------------------------------------------")
    keep_samples <- !colnames(dds) %in% exclude_samples
    dds_sub <- dds[, keep_samples]
    if (add_mouse_covariates) {
        dds_sub <- add_mouse_covariates_to_dds(dds_sub, meta_mouse)
    }
    dds_sub$group <- factor(dds_sub$group, levels = c("Ctrl", "PPI"))
    dds_sub$group <- relevel(dds_sub$group, ref = "Ctrl")
    if ("isolator" %in% colnames(as.data.frame(colData(dds_sub)))) {
        dds_sub$isolator <- droplevels(factor(dds_sub$isolator))
    }
    if ("cage_position" %in% colnames(as.data.frame(colData(dds_sub)))) {
        dds_sub$cage_position <- droplevels(factor(dds_sub$cage_position))
    }
    if ("age_weeks" %in% colnames(as.data.frame(colData(dds_sub)))) {
        dds_sub$age_z <- as.numeric(scale(dds_sub$age_weeks))
    }
    design_vars <- all.vars(design_formula)
    missing_design_vars <- setdiff(design_vars, colnames(as.data.frame(colData(dds_sub))))
    if (length(missing_design_vars) > 0) {
        warning("Skipping ", organ_name, " | ", model_label, ": missing design variables: ", paste(missing_design_vars, 
            collapse = ", "))
        return(list(failed = TRUE, reason = paste0("missing_design_vars: ", paste(missing_design_vars, 
            collapse = ", ")), organ = organ_name, model = model_label, dds = NULL, res = NULL, res_df = NULL))
    }
    cov_complete <- stats::complete.cases(as.data.frame(colData(dds_sub))[, design_vars, drop = FALSE])
    if (!all(cov_complete)) {
        warning("Removing ", sum(!cov_complete), " samples with missing covariates in ", organ_name, 
            " | ", model_label)
        dds_sub <- dds_sub[, cov_complete]
    }
    design(dds_sub) <- design_formula
    design_matrix <- model.matrix(design(dds_sub), data = as.data.frame(colData(dds_sub)))
    design_rank <- qr(design_matrix)$rank
    n_coef <- ncol(design_matrix)
    if (design_rank < n_coef) {
        warning("Skipping ", organ_name, " | ", model_label, ": design is not full rank.")
        return(list(failed = TRUE, reason = "design_not_full_rank", organ = organ_name, model = model_label, 
            dds = NULL, res = NULL, res_df = NULL, design_matrix = design_matrix, design_rank = design_rank, 
            n_coef = n_coef))
    }
    keep_genes <- rowSums(DESeq2::counts(dds_sub) >= min_count) >= ceiling(min_prop_samples * ncol(dds_sub))
    dds_sub <- dds_sub[keep_genes, ]
    dds_sub <- DESeq2::DESeq(dds_sub)
    rn <- DESeq2::resultsNames(dds_sub)
    message("resultsNames:")
    print(rn)
    if (!is.null(result_name)) {
        if (!result_name %in% rn) {
            warning("Skipping ", organ_name, " | ", model_label, ": result_name ", result_name, " not found.")
            return(list(failed = TRUE, reason = paste0("result_name_not_found: ", result_name), organ = organ_name, 
                model = model_label, dds = dds_sub, res = NULL, res_df = NULL, results_names = rn))
        }
        res <- DESeq2::results(dds_sub, name = result_name, alpha = alpha)
    }
    else {
        res <- DESeq2::results(dds_sub, contrast = contrast, alpha = alpha)
    }
    res_df <- as.data.frame(res) %>% tibble::rownames_to_column("gene") %>% dplyr::mutate(organ = organ_name, 
        model = model_label, sig_0.1 = !is.na(padj) & padj < 0.1, sig_0.05 = !is.na(padj) & padj < 0.05) %>% 
        dplyr::arrange(padj)
    list(failed = FALSE, reason = NA_character_, organ = organ_name, model = model_label, design_formula = design_formula, 
        design_rank = design_rank, n_coef = n_coef, results_names = rn, dds = dds_sub, res = res, res_df = res_df)
}
model_specs <- list(M0_unadjusted = list(design_formula = ~group, contrast = c("group", "PPI", "Ctrl"), 
    result_name = NULL), M1_age_adjusted = list(design_formula = ~age_z + group, contrast = c("group", 
    "PPI", "Ctrl"), result_name = NULL), M2_isolator_adjusted = list(design_formula = ~isolator + group, 
    contrast = c("group", "PPI", "Ctrl"), result_name = NULL), M3_age_isolator_adjusted = list(design_formula = ~age_z + 
    isolator + group, contrast = c("group", "PPI", "Ctrl"), result_name = NULL))
deseq2_sensitivity <- list()
for (model_name in names(model_specs)) {
    spec <- model_specs[[model_name]]
    deseq2_sensitivity[[model_name]] <- lapply(names(dds_list), function(org) {
        run_deseq2_sensitivity_one(dds = dds_list[[org]], organ_name = org, model_label = model_name, 
            design_formula = spec$design_formula, contrast = spec$contrast, result_name = spec$result_name, 
            add_mouse_covariates = TRUE, meta_mouse = meta_mouse, exclude_samples = exclude_samples, 
            alpha = alpha)
    })
    names(deseq2_sensitivity[[model_name]]) <- names(dds_list)
}
for (model_name in names(deseq2_sensitivity)) {
    for (org in names(deseq2_sensitivity[[model_name]])) {
        obj <- deseq2_sensitivity[[model_name]][[org]]
        if (!isTRUE(obj$failed) && !is.null(obj$res_df)) {
            write.csv(obj$res_df, ppi_output(file.path(out_dir, paste0("DESeq2_", org, "_", model_name, 
                ".csv"))), row.names = FALSE)
        }
    }
}
direct_age_effect <- lapply(names(dds_list), function(org) {
    run_deseq2_sensitivity_one(dds = dds_list[[org]], organ_name = org, model_label = "AGE_direct_effect_adjusted_for_group", 
        design_formula = ~group + age_z, contrast = NULL, result_name = "age_z", add_mouse_covariates = TRUE, 
        meta_mouse = meta_mouse, exclude_samples = exclude_samples, alpha = alpha)
})
names(direct_age_effect) <- names(dds_list)
for (org in names(direct_age_effect)) {
    obj <- direct_age_effect[[org]]
    if (!isTRUE(obj$failed) && !is.null(obj$res_df)) {
        write.csv(obj$res_df, ppi_output(file.path(out_dir, paste0("DESeq2_", org, "_direct_age_effect_adjusted_for_group.csv"))), 
            row.names = FALSE)
    }
}
summarise_deseq2_models <- function(model_list) {
    dplyr::bind_rows(lapply(names(model_list), function(model_name) {
        dplyr::bind_rows(lapply(names(model_list[[model_name]]), function(org) {
            obj <- model_list[[model_name]][[org]]
            if (isTRUE(obj$failed) || is.null(obj$res_df)) {
                return(tibble::tibble(model = model_name, organ = org, failed = TRUE, reason = obj$reason, 
                  n_sig_0.1 = NA_integer_, n_sig_0.05 = NA_integer_, min_padj = NA_real_, median_padj = NA_real_))
            }
            res <- obj$res_df
            tibble::tibble(model = model_name, organ = org, failed = FALSE, reason = NA_character_, n_sig_0.1 = sum(res$padj < 
                0.1, na.rm = TRUE), n_sig_0.05 = sum(res$padj < 0.05, na.rm = TRUE), min_padj = min(res$padj, 
                na.rm = TRUE), median_padj = median(res$padj, na.rm = TRUE))
        }))
    }))
}
sensitivity_summary <- summarise_deseq2_models(deseq2_sensitivity)
print(sensitivity_summary)
write.csv(sensitivity_summary, ppi_output(file.path(out_dir, "DESeq2_sensitivity_DEG_count_summary.csv"), 
    name = "de_seq2_sensitivity_deg_count_summary.csv"), row.names = FALSE)
direct_age_summary <- dplyr::bind_rows(lapply(names(direct_age_effect), function(org) {
    obj <- direct_age_effect[[org]]
    if (isTRUE(obj$failed) || is.null(obj$res_df)) {
        return(tibble::tibble(organ = org, failed = TRUE, reason = obj$reason, n_age_sig_0.1 = NA_integer_, 
            n_age_sig_0.05 = NA_integer_, min_padj = NA_real_, median_padj = NA_real_))
    }
    res <- obj$res_df
    tibble::tibble(organ = org, failed = FALSE, reason = NA_character_, n_age_sig_0.1 = sum(res$padj < 
        0.1, na.rm = TRUE), n_age_sig_0.05 = sum(res$padj < 0.05, na.rm = TRUE), min_padj = min(res$padj, 
        na.rm = TRUE), median_padj = median(res$padj, na.rm = TRUE))
}))
print(direct_age_summary)
write.csv(direct_age_summary, ppi_output(file.path(out_dir, "DESeq2_direct_age_effect_summary.csv"), 
    name = "de_seq2_direct_age_effect_summary.csv"), row.names = FALSE)
compare_model_to_unadjusted <- function(deseq2_sensitivity, reference_model = "M0_unadjusted", alpha = 0.1) {
    ref_list <- deseq2_sensitivity[[reference_model]]
    dplyr::bind_rows(lapply(names(deseq2_sensitivity), function(model_name) {
        test_list <- deseq2_sensitivity[[model_name]]
        dplyr::bind_rows(lapply(names(ref_list), function(org) {
            ref_obj <- ref_list[[org]]
            test_obj <- test_list[[org]]
            if (isTRUE(ref_obj$failed) || isTRUE(test_obj$failed)) {
                return(tibble::tibble(model = model_name, reference = reference_model, organ = org, failed = TRUE, 
                  n_sig_ref = NA_integer_, n_sig_test = NA_integer_, overlap_sig = NA_integer_, only_ref = NA_integer_, 
                  only_test = NA_integer_, cor_log2FC = NA_real_, cor_neglog10_p = NA_real_, median_abs_delta_log2FC = NA_real_, 
                  max_abs_delta_log2FC = NA_real_))
            }
            ref <- ref_obj$res_df %>% dplyr::select(gene, ref_log2FC = log2FoldChange, ref_pvalue = pvalue, 
                ref_padj = padj)
            test <- test_obj$res_df %>% dplyr::select(gene, test_log2FC = log2FoldChange, test_pvalue = pvalue, 
                test_padj = padj)
            joined <- ref %>% dplyr::inner_join(test, by = "gene")
            joined %>% dplyr::summarise(model = model_name, reference = reference_model, organ = org, 
                failed = FALSE, n_sig_ref = sum(ref_padj < alpha, na.rm = TRUE), n_sig_test = sum(test_padj < 
                  alpha, na.rm = TRUE), overlap_sig = sum(ref_padj < alpha & test_padj < alpha, na.rm = TRUE), 
                only_ref = sum(ref_padj < alpha & !(test_padj < alpha), na.rm = TRUE), only_test = sum(!(ref_padj < 
                  alpha) & test_padj < alpha, na.rm = TRUE), cor_log2FC = cor(ref_log2FC, test_log2FC, 
                  use = "pairwise.complete.obs"), cor_neglog10_p = cor(-log10(ref_pvalue), -log10(test_pvalue), 
                  use = "pairwise.complete.obs"), median_abs_delta_log2FC = median(abs(test_log2FC - 
                  ref_log2FC), na.rm = TRUE), max_abs_delta_log2FC = max(abs(test_log2FC - ref_log2FC), 
                  na.rm = TRUE), .groups = "drop")
        }))
    }))
}
sensitivity_comparison <- compare_model_to_unadjusted(deseq2_sensitivity = deseq2_sensitivity, reference_model = "M0_unadjusted", 
    alpha = alpha)
print(sensitivity_comparison)
write.csv(sensitivity_comparison, ppi_output(file.path(out_dir, "DESeq2_sensitivity_comparison_vs_unadjusted.csv"), 
    name = "de_seq2_sensitivity_comparison_vs_unadjusted.csv"), row.names = FALSE)
genes_of_interest <- c("ENSMUSG00000026822", "ENSMUSG00000041828", "ENSMUSG00000004328", "ENSMUSG00000016496", 
    "ENSMUSG00000030281", "ENSMUSG00000024067", "ENSMUSG00000026072", "ENSMUSG00000055782", "ENSMUSG00000031765", 
    "ENSMUSG00000031762")
extract_genes_from_deseq2_sensitivity <- function(deseq2_sensitivity, genes_of_interest) {
    dplyr::bind_rows(lapply(names(deseq2_sensitivity), function(model_name) {
        dplyr::bind_rows(lapply(names(deseq2_sensitivity[[model_name]]), function(org) {
            obj <- deseq2_sensitivity[[model_name]][[org]]
            if (isTRUE(obj$failed) || is.null(obj$res_df)) {
                return(NULL)
            }
            obj$res_df %>% dplyr::filter(gene %in% genes_of_interest) %>% dplyr::mutate(model = model_name, 
                organ = org) %>% dplyr::select(model, organ, gene, log2FoldChange, lfcSE, stat, pvalue, 
                padj, sig_0.1, sig_0.05)
        }))
    }))
}
genes_of_interest_sensitivity <- extract_genes_from_deseq2_sensitivity(deseq2_sensitivity = deseq2_sensitivity, 
    genes_of_interest = genes_of_interest)
print(genes_of_interest_sensitivity)
write.csv(genes_of_interest_sensitivity, ppi_output(file.path(out_dir, "DESeq2_genes_of_interest_sensitivity.csv"), 
    name = "de_seq2_genes_of_interest_sensitivity.csv"), row.names = FALSE)
all_sensitivity_results <- dplyr::bind_rows(lapply(names(deseq2_sensitivity), function(model_name) {
    dplyr::bind_rows(lapply(names(deseq2_sensitivity[[model_name]]), function(org) {
        obj <- deseq2_sensitivity[[model_name]][[org]]
        if (isTRUE(obj$failed) || is.null(obj$res_df)) {
            return(NULL)
        }
        obj$res_df
    }))
}))
write.csv(all_sensitivity_results, ppi_output(file.path(out_dir, "DESeq2_all_sensitivity_results_combined.csv"), 
    name = "de_seq2_all_sensitivity_results_combined.csv"), row.names = FALSE)
direct_age_all_results <- dplyr::bind_rows(lapply(names(direct_age_effect), function(org) {
    obj <- direct_age_effect[[org]]
    if (isTRUE(obj$failed) || is.null(obj$res_df)) {
        return(NULL)
    }
    obj$res_df
}))
write.csv(direct_age_all_results, ppi_output(file.path(out_dir, "DESeq2_direct_age_effect_all_results_combined.csv"), 
    name = "de_seq2_direct_age_effect_all_results_combined.csv"), row.names = FALSE)
failed_models <- dplyr::bind_rows(lapply(names(deseq2_sensitivity), function(model_name) {
    dplyr::bind_rows(lapply(names(deseq2_sensitivity[[model_name]]), function(org) {
        obj <- deseq2_sensitivity[[model_name]][[org]]
        if (!isTRUE(obj$failed)) 
            return(NULL)
        tibble::tibble(model = model_name, organ = org, reason = obj$reason)
    }))
}))
write.csv(failed_models, ppi_output(file.path(out_dir, "DESeq2_failed_sensitivity_models.csv"), name = "de_seq2_failed_sensitivity_models.csv"), 
    row.names = FALSE)
writexl::write_xlsx(list(DEG_count_summary = sensitivity_summary, Comparison_vs_unadjusted = sensitivity_comparison, 
    Direct_age_effect = direct_age_summary, Genes_of_interest = genes_of_interest_sensitivity, Failed_models = failed_models), 
    ppi_output(file.path(out_dir, "DESeq2_covariate_sensitivity_summary.xlsx"), name = "de_seq2_covariate_sensitivity_summary.xlsx"))
saveRDS(list(deseq2_sensitivity = deseq2_sensitivity, direct_age_effect = direct_age_effect, sensitivity_summary = sensitivity_summary, 
    sensitivity_comparison = sensitivity_comparison, direct_age_summary = direct_age_summary, genes_of_interest_sensitivity = genes_of_interest_sensitivity, 
    meta_mouse = meta_mouse, mouse_cov = mouse_cov, exclude_samples = exclude_samples, alpha = alpha), 
    ppi_output(file.path(out_dir, "DESeq2_covariate_sensitivity_full_objects.rds"), name = "de_seq2_covariate_sensitivity_full_objects.rds"))
message("============================================================")
message("DESeq2 covariate sensitivity analysis complete.")
message("Results saved to: ", out_dir)
message("============================================================")
