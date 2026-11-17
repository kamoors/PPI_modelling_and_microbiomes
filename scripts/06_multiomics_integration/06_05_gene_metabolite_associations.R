# SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0
# Gene-metabolite DESeq2/GSEA/ORA

source("scripts/helpers/setup.R")
invisible(run_analysis(
  section = "06_multiomics_integration",
  directory = "integration",
  script = "06_05_gene_metabolite_associations.R",
  restricted = FALSE,
  code = {
message("Loaded 06_05_gene_metabolite_associations.R")
message("Model: counts ~ Group + Metabolite only; no donor-aware or LODO components.")
required_packages <- c("DESeq2", "SummarizedExperiment", "fgsea", "clusterProfiler", "msigdbr", "AnnotationDbi",
    "org.Mm.eg.db", "dplyr", "tidyr", "purrr", "readr", "readxl", "stringr", "tibble", "writexl")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0L) {
    stop("Missing required R packages: ", paste(missing_packages, collapse = ", "), "\nInstall missing Bioconductor packages with BiocManager::install() ",
        "and missing CRAN packages with install.packages(), then rerun.")
}
suppressPackageStartupMessages({
    library(DESeq2)
    library(SummarizedExperiment)
    library(fgsea)
    library(clusterProfiler)
    library(msigdbr)
    library(AnnotationDbi)
    library(org.Mm.eg.db)
    library(dplyr)
    library(tidyr)
    library(purrr)
    library(readr)
    library(readxl)
    library(stringr)
    library(tibble)
    library(writexl)
})
options(stringsAsFactors = FALSE)
set.seed(20260717)
script_path <- tryCatch(project_normalize_path(sys.frame(1)$ofile, winslash = "/", mustWork = FALSE), error = function(e) NA_character_)
script_dir <- if (!is.na(script_path) && nzchar(script_path)) {
    dirname(script_path)
} else {
    getwd()
}
config_candidates <- c(project_code("scripts/06_multiomics_integration/06_01_configure_multiomics.R"), file.path(script_dir,
    "PPI_multiomics_correlations_00_config_v2_fixed.R"), "06_01_configure_multiomics.R", "PPI_multiomics_correlations_00_config_v2_fixed.R")
config_path <- config_candidates[project_file_exists(config_candidates)][1]
if (is.na(config_path)) {
    stop("Could not find the multi-omics configuration script. Place this script ", "beside 06_01_configure_multiomics.R or ",
        "PPI_multiomics_correlations_00_config_v2_fixed.R.")
}
project_source(config_path)
prepared_file <- file.path(association_root, "01_prepared_inputs", "PPI_multiomics_prepared_inputs.rds")
if (!project_file_exists(prepared_file)) {
    stop("Prepared multi-omics input bundle not found: ", prepared_file, "\nRun the corrected input-preparation script first.")
}
bundle <- readRDS(project_input(prepared_file))
output_root <- file.path(association_root, "05_gene_metabolite_DESeq2_Group_Metabolite_Hallmark_Reactome")
gene_result_dir <- file.path(output_root, "complete_gene_results")
gsea_result_dir <- file.path(output_root, "gsea_family_tables")
ora_result_dir <- file.path(output_root, "ora_family_tables")
audit_dir <- file.path(output_root, "audits")
for (directory in c(output_root, gene_result_dir, gsea_result_dir, ora_result_dir, audit_dir)) {
    dir.create(project_output(directory, directory = TRUE), recursive = TRUE, showWarnings = FALSE)
}
gm_settings <- list(metabolite_scope = "significant", metabolite_differential_fdr_cutoff = 0.1, gene_fdr_cutoff = 0.1,
    gsea_fdr_cutoff = 0.1, ora_fdr_cutoff = 0.1, ora_gene_fdr_cutoff = 0.1, ora_min_genes_per_direction = 10L,
    minimum_samples = 8L, minimum_samples_per_group = 3L, minimum_metabolite_sd = 1e-08, tissue_map = list(serum_lcms = c("Liver",
        "Colon", "Brain"), cecum_lcms = c("Liver", "Colon", "Brain"), feces_lcms = c("Liver", "Colon",
        "Brain"), cecum_nmr = c("Liver", "Colon", "Brain")), exclude_rna_samples = c("25Jul1866-C02",
        "25Jul1866-G04"), hallmark_min_size = 15L, hallmark_max_size = 300L, reactome_min_size = 15L,
    reactome_max_size = 500L, independent_filtering = TRUE, cooks_cutoff = FALSE)
allowed_metabolite_scopes <- c("significant", "all_tested", "all_matrix")
if (!gm_settings$metabolite_scope %in% allowed_metabolite_scopes) {
    stop("gm_settings$metabolite_scope must be one of: ", paste(allowed_metabolite_scopes, collapse = ", "))
}
setting_value_to_text <- function(x) {
    if (is.null(x))
        return(NA_character_)
    if (inherits(x, "formula"))
        return(paste(deparse(x), collapse = " "))
    if (is.function(x))
        return("<function>")
    if (is.list(x)) {
        if (length(x) == 0L)
            return("")
        values <- vapply(x, setting_value_to_text, character(1))
        if (!is.null(names(x))) {
            item_names <- names(x)
            item_names[is.na(item_names) | !nzchar(item_names)] <- paste0("item", which(is.na(item_names) |
                !nzchar(item_names)))
            values <- paste0(item_names, "=", values)
        }
        return(paste(values, collapse = ";"))
    }
    paste(as.character(x), collapse = ";")
}
settings_audit <- tibble::enframe(gm_settings, name = "setting", value = "raw_value") %>% dplyr::mutate(value = purrr::map_chr(.data$raw_value,
    setting_value_to_text)) %>% dplyr::select(.data$setting, .data$value)
readr::write_csv(settings_audit, project_output(file.path(audit_dir, "gene_metabolite_analysis_settings.csv"),
    name = "gene_metabolite_analysis_settings.csv"))
`%||%` <- function(x, y) {
    if (is.null(x) || length(x) == 0L)
        y
    else x
}
safe_file_id <- function(x) {
    stringr::str_replace_all(as.character(x), "[^A-Za-z0-9._-]+", "_")
}
strip_ensembl_version <- function(x) {
    sub("\\.[0-9]+$", "", as.character(x))
}
normalize_group <- function(x) {
    x <- stringr::str_to_lower(stringr::str_squish(as.character(x)))
    dplyr::case_when(stringr::str_detect(x, "^ctrl|^control|^ctl") ~ "ctrl", stringr::str_detect(x, "^ppi") ~
        "ppi", TRUE ~ NA_character_)
}
normalize_organ_name <- function(x) {
    x <- stringr::str_to_lower(as.character(x))
    dplyr::case_when(stringr::str_detect(x, "liver") ~ "Liver", stringr::str_detect(x, "colon") ~ "Colon",
        stringr::str_detect(x, "brain") ~ "Brain", TRUE ~ stringr::str_to_title(x))
}
first_existing_column <- function(data, candidates, required = TRUE) {
    hit <- candidates[candidates %in% names(data)]
    if (length(hit) == 0L) {
        if (required) {
            stop("None of the expected columns were present: ", paste(candidates, collapse = ", "))
        }
        return(NULL)
    }
    hit[[1]]
}
resolve_sheet <- function(path, candidates) {
    sheets <- readxl::excel_sheets(project_input(path))
    hit <- candidates[candidates %in% sheets]
    if (length(hit) == 0L) {
        stop("None of the expected sheets were found in ", path, ". Expected: ", paste(candidates, collapse = ", "),
            ". Available: ", paste(sheets, collapse = ", "))
    }
    hit[[1]]
}
normalize_file_basename <- function(x) {
    ext <- tolower(tools::file_ext(x))
    stem <- tools::file_path_sans_ext(basename(x))
    stem <- stringr::str_remove(stem, " ?\\([0-9]+\\)$")
    paste0(tolower(stem), ifelse(nzchar(ext), paste0(".", ext), ""))
}
get_prepared_resolved_path <- function(input_name, fallback = NA_character_) {
    audit <- as.data.frame(bundle$resolved_paths)
    idx <- match(input_name, audit$input)
    if (!is.na(idx)) {
        hit <- as.character(audit$resolved_path[[idx]])
        if (!is.na(hit) && nzchar(hit))
            return(hit)
    }
    fallback
}
resolve_near_input <- function(primary_input, filename_candidates, search_root) {
    filename_candidates <- unique(filename_candidates)
    nearby <- file.path(dirname(primary_input), filename_candidates)
    nearby_hit <- nearby[project_file_exists(nearby)]
    if (length(nearby_hit) >= 1L) {
        return(project_normalize_path(nearby_hit[[1]], winslash = "/"))
    }
    if (!project_dir_exists(search_root)) {
        stop("Search root does not exist: ", search_root)
    }
    all_files <- project_list_files(search_root, recursive = TRUE, full.names = TRUE)
    desired <- normalize_file_basename(filename_candidates)
    hits <- all_files[normalize_file_basename(all_files) %in% desired]
    if (length(hits) == 0L) {
        stop("Could not locate any of these files near ", primary_input, " or under ", search_root, ": ",
            paste(filename_candidates, collapse = ", "))
    }
    input_dir <- project_normalize_path(dirname(primary_input), winslash = "/", mustWork = FALSE)
    hit_dirs <- project_normalize_path(dirname(hits), winslash = "/", mustWork = FALSE)
    same_dir <- which(hit_dirs == input_dir)
    if (length(same_dir) >= 1L) {
        return(project_normalize_path(hits[same_dir[[1]]], winslash = "/"))
    }
    if (length(hits) > 1L) {
        stop("Multiple candidate result files were found and none was beside the ", "prepared input. Specify the exact path in the script. Candidates:\n- ",
            paste(hits, collapse = "\n- "))
    }
    project_normalize_path(hits[[1]], winslash = "/")
}
match_feature_name_detailed <- function(target, available_names) {
    target <- as.character(target)[1]
    available_names <- as.character(available_names)
    if (is.na(target) || !nzchar(target) || length(available_names) == 0L) {
        return(list(index = NA_integer_, status = "not_matched"))
    }
    strategies <- list(exact = list(target, available_names), trimmed = list(stringr::str_squish(target),
        stringr::str_squish(available_names)), ascii = list(iconv(stringr::str_squish(target), from = "",
        to = "ASCII//TRANSLIT"), iconv(stringr::str_squish(available_names), from = "", to = "ASCII//TRANSLIT")))
    for (strategy_name in names(strategies)) {
        values <- strategies[[strategy_name]]
        idx <- which(values[[2]] == values[[1]])
        if (length(idx) == 1L) {
            return(list(index = idx[[1]], status = strategy_name))
        }
        if (length(idx) > 1L) {
            return(list(index = NA_integer_, status = paste0("ambiguous_", strategy_name)))
        }
    }
    target_ascii <- iconv(stringr::str_squish(target), from = "", to = "ASCII//TRANSLIT")
    available_ascii <- iconv(stringr::str_squish(available_names), from = "", to = "ASCII//TRANSLIT")
    target_safe <- make.names(tolower(target_ascii), unique = FALSE)
    available_safe <- make.names(tolower(available_ascii), unique = FALSE)
    idx <- which(available_safe == target_safe)
    if (length(idx) == 1L) {
        return(list(index = idx[[1]], status = "normalized"))
    }
    if (length(idx) > 1L) {
        return(list(index = NA_integer_, status = "ambiguous_normalized"))
    }
    list(index = NA_integer_, status = "not_matched")
}
safe_zscore <- function(x) {
    x <- as.numeric(x)
    mu <- mean(x, na.rm = TRUE)
    s <- stats::sd(x, na.rm = TRUE)
    if (!is.finite(s) || s <= gm_settings$minimum_metabolite_sd) {
        return(rep(NA_real_, length(x)))
    }
    (x - mu)/s
}
is_full_rank_formula <- function(formula, data) {
    model_matrix <- tryCatch(stats::model.matrix(formula, data = data), error = function(e) NULL)
    !is.null(model_matrix) && qr(model_matrix)$rank == ncol(model_matrix)
}
design_condition_number <- function(formula, data) {
    model_matrix <- tryCatch(stats::model.matrix(formula, data = data), error = function(e) NULL)
    if (is.null(model_matrix))
        return(NA_real_)
    tryCatch(kappa(model_matrix, exact = TRUE), error = function(e) NA_real_)
}
rna_loader_path <- get_prepared_resolved_path("rna_loader_script", input_paths$rna_loader_script)
if (is.na(rna_loader_path) || !project_file_exists(rna_loader_path)) {
    stop("RNA loader script not found: ", rna_loader_path)
}
rna_env <- new.env(parent = globalenv())
project_sys_source(rna_loader_path, envir = rna_env)
is_dds_list <- function(x) {
    is.list(x) && length(x) > 0L && all(vapply(x, inherits, logical(1), what = "DESeqDataSet"))
}
extract_dds_from_nested <- function(x) {
    if (!is.list(x) || length(x) == 0L)
        return(NULL)
    out <- purrr::map(x, function(z) {
        if (inherits(z, "DESeqDataSet"))
            return(z)
        if (is.list(z)) {
            candidates <- c("dds", "DESeqDataSet", "dds_no_outliers")
            for (candidate in candidates) {
                if (!is.null(z[[candidate]]) && inherits(z[[candidate]], "DESeqDataSet")) {
                  return(z[[candidate]])
                }
            }
        }
        NULL
    })
    keep <- !vapply(out, is.null, logical(1))
    out <- out[keep]
    if (length(out) == 0L)
        NULL
    else out
}
candidate_names <- c("dds_list", "dds_list_no_outliers", "dds_no_outliers", "deseq_no_outliers")
dds_list <- NULL
for (candidate_name in candidate_names) {
    if (!exists(candidate_name, envir = rna_env, inherits = FALSE))
        next
    object <- get(candidate_name, envir = rna_env, inherits = FALSE)
    if (is_dds_list(object)) {
        dds_list <- object
        break
    }
    nested <- extract_dds_from_nested(object)
    if (!is.null(nested) && is_dds_list(nested)) {
        dds_list <- nested
        break
    }
}
if (is.null(dds_list)) {
    for (object_name in ls(rna_env, all.names = TRUE)) {
        object <- get(object_name, envir = rna_env, inherits = FALSE)
        if (is_dds_list(object)) {
            dds_list <- object
            break
        }
        nested <- extract_dds_from_nested(object)
        if (!is.null(nested) && is_dds_list(nested)) {
            dds_list <- nested
            break
        }
    }
}
if (is.null(dds_list) || !is_dds_list(dds_list)) {
    stop("Could not identify an organ-specific list of DESeqDataSet objects after ", "sourcing: ", rna_loader_path)
}
if (is.null(names(dds_list)) || any(names(dds_list) == "")) {
    inferred_names <- purrr::map_chr(dds_list, function(dds) {
        col_data <- as.data.frame(SummarizedExperiment::colData(dds))
        organ_col <- first_existing_column(col_data, c("organ", "Organ", "tissue"), required = FALSE)
        if (is.null(organ_col))
            return(NA_character_)
        normalize_organ_name(unique(col_data[[organ_col]])[[1]])
    })
    names(dds_list) <- inferred_names
} else {
    names(dds_list) <- normalize_organ_name(names(dds_list))
}
required_organs <- unique(unlist(gm_settings$tissue_map, use.names = FALSE))
missing_organs <- setdiff(required_organs, names(dds_list))
if (length(missing_organs) > 0L) {
    stop("RNA loader did not provide DESeqDataSet objects for: ", paste(missing_organs, collapse = ", "))
}
dds_list <- purrr::imap(dds_list, function(dds, organ) {
    keep <- !colnames(dds) %in% gm_settings$exclude_rna_samples
    dds[, keep]
})
rna_object_audit <- purrr::imap_dfr(dds_list, function(dds, organ) {
    tibble::tibble(Organ = organ, n_samples = ncol(dds), n_genes = nrow(dds), excluded_samples_still_present = {
        still_present <- intersect(colnames(dds), gm_settings$exclude_rna_samples)
        if (length(still_present) == 0L) {
            "none"
        }
        else {
            paste(still_present, collapse = ";")
        }
    })
})
readr::write_csv(rna_object_audit, project_output(file.path(audit_dir, "rna_object_audit.csv"), name = "rna_object_audit.csv"))
requested_mouse_layers <- names(gm_settings$tissue_map)
missing_prepared_layers <- setdiff(requested_mouse_layers, names(bundle$metabolomics))
if (length(missing_prepared_layers) > 0L) {
    stop("Prepared bundle is missing requested mouse metabolomics layers: ", paste(missing_prepared_layers,
        collapse = ", "))
}
analysis_datasets <- intersect(requested_mouse_layers, bundle$mouse_metabolomics_datasets)
missing_from_mouse_registry <- setdiff(requested_mouse_layers, analysis_datasets)
if (length(missing_from_mouse_registry) > 0L) {
    stop("These requested layers are not registered as mouse metabolomics datasets ", "in the prepared bundle: ",
        paste(missing_from_mouse_registry, collapse = ", "))
}
all_bundle_layers <- names(bundle$metabolomics)
excluded_nonmouse_layers <- setdiff(all_bundle_layers, analysis_datasets)
nonmouse_layer_audit <- tibble::tibble(dataset_id = excluded_nonmouse_layers, included_in_gene_metabolite_model = rep(FALSE,
    length(excluded_nonmouse_layers)), reason = rep(paste0("No one-to-one recipient-mouse/RNA sample alignment; ",
    "the counts ~ Group + Metabolite model requires matched mouse samples."), length(excluded_nonmouse_layers)))
readr::write_csv(nonmouse_layer_audit, project_output(file.path(audit_dir, "excluded_nonmouse_metabolomics_layers.csv"),
    name = "excluded_nonmouse_metabolomics_layers.csv"))
serum_input <- get_prepared_resolved_path("serum_lcms", input_paths$serum_lcms)
cecum_nmr_input <- get_prepared_resolved_path("cecum_nmr", input_paths$cecum_nmr)
serum_result_file <- resolve_near_input(serum_input, c("Serum_limma_donor_and_duplicateCorrelation_results.xlsx"),
    metabolomics_root)
cecum_feces_combined_result_file <- file.path(metabolomics_root, "20250714_MassSpecData", "Cecum_Feces_cleaned_after_collaborator_review",
    "Limma_analysis", "Cecum_Feces_LCMS_combined_compartment_limma_FINAL_results.xlsx")
if (!project_file_exists(cecum_feces_combined_result_file)) {
    stop("Combined caecum/faeces differential-result workbook not found: ", cecum_feces_combined_result_file)
}
read_metabolite_selection <- function(dataset_id) {
    dataset <- bundle$metabolomics[[dataset_id]]
    if (is.null(dataset)) {
        stop("Prepared metabolomics dataset missing: ", dataset_id)
    }
    available <- colnames(dataset$matrix)
    if (dataset_id == "serum_lcms") {
        sheet <- resolve_sheet(serum_result_file, c("Mouse_duplicateCorrelation", "Donor_level_limma"))
        raw <- readxl::read_excel(project_input(serum_result_file), sheet = sheet)
        feature_col <- first_existing_column(raw, c("metabolite", "analyte"))
        effect_col <- first_existing_column(raw, c("log2_FC", "log2FC", "effect"), required = FALSE)
        p_col <- first_existing_column(raw, c("p_value", "p.value", "P.Value"), required = FALSE)
        fdr_col <- first_existing_column(raw, c("FDR", "adj.P.Val", "padj"))
        source_path <- serum_result_file
        source_component <- sheet
    }
    else if (dataset_id %in% c("cecum_lcms", "feces_lcms")) {
        sheet <- resolve_sheet(cecum_feces_combined_result_file, c("Sheet1"))
        raw_all <- readxl::read_excel(project_input(cecum_feces_combined_result_file), sheet = sheet)
        required_filter_columns <- c("representation", "model", "contrast")
        missing_filter_columns <- setdiff(required_filter_columns, names(raw_all))
        if (length(missing_filter_columns) > 0L) {
            stop("Combined caecum/faeces workbook is missing columns: ", paste(missing_filter_columns,
                collapse = ", "))
        }
        raw <- raw_all %>% dplyr::filter(stringr::str_to_lower(as.character(.data$representation)) ==
            "one_per_metabolite", as.character(.data$model) == "all_samples_donor_blocked_sensitivity",
            as.character(.data$contrast) == "PPI_vs_control_adjusted_for_compartment")
        if (nrow(raw) == 0L) {
            stop("No rows matched the requested combined caecum/faeces model.")
        }
        feature_col <- first_existing_column(raw, c("analyte", "metabolite"))
        effect_col <- first_existing_column(raw, c("log2_FC", "log2FC", "effect"), required = FALSE)
        p_col <- first_existing_column(raw, c("p_value", "p.value", "P.Value"), required = FALSE)
        fdr_col <- first_existing_column(raw, c("FDR", "adj.P.Val", "padj"))
        source_path <- cecum_feces_combined_result_file
        source_component <- paste0(sheet, "; representation=one_per_metabolite; ", "model=all_samples_donor_blocked_sensitivity; ",
            "contrast=PPI_vs_control_adjusted_for_compartment; ", "applied_to=", dataset_id)
    }
    else if (dataset_id == "cecum_nmr") {
        sheet <- resolve_sheet(cecum_nmr_input, c("Mouse_dupCor_sensitivity", "Donor_averaged_primary"))
        raw <- readxl::read_excel(project_input(cecum_nmr_input), sheet = sheet)
        feature_col <- first_existing_column(raw, c("analysis_group", "metabolite"))
        effect_col <- first_existing_column(raw, c("score_difference", "log2_FC", "effect"), required = FALSE)
        p_col <- first_existing_column(raw, c("p_value", "p.value", "P.Value"), required = FALSE)
        fdr_col <- first_existing_column(raw, c("FDR", "adj.P.Val", "padj"))
        source_path <- cecum_nmr_input
        source_component <- sheet
    }
    else {
        stop("No differential-result source configured for: ", dataset_id)
    }
    source_feature <- as.character(raw[[feature_col]])
    details <- purrr::map(source_feature, match_feature_name_detailed, available_names = available)
    matched_index <- purrr::map_int(details, "index")
    matching_status <- purrr::map_chr(details, "status")
    matched_ok <- !is.na(matched_index)
    source_table <- tibble::tibble(dataset_id = dataset_id, source_feature = source_feature, metabolite = dplyr::if_else(matched_ok,
        available[matched_index], NA_character_), matching_status = matching_status, differential_effect = if (!is.null(effect_col)) {
        suppressWarnings(as.numeric(raw[[effect_col]]))
    }
    else {
        NA_real_
    }, differential_p_value = if (!is.null(p_col)) {
        suppressWarnings(as.numeric(raw[[p_col]]))
    }
    else {
        NA_real_
    }, differential_FDR = suppressWarnings(as.numeric(raw[[fdr_col]])), source_path = source_path, source_component = source_component) %>%
        dplyr::mutate(differential_significant = is.finite(.data$differential_FDR) & .data$differential_FDR <
            gm_settings$metabolite_differential_fdr_cutoff)
    best <- source_table %>% dplyr::filter(!is.na(.data$metabolite)) %>% dplyr::arrange(.data$differential_FDR,
        .data$differential_p_value) %>% dplyr::distinct(.data$metabolite, .keep_all = TRUE)
    reverse_audit <- tibble::tibble(dataset_id = dataset_id, metabolite = available) %>% dplyr::left_join(best %>%
        dplyr::select(.data$dataset_id, .data$metabolite, .data$source_feature, .data$matching_status,
            .data$differential_effect, .data$differential_p_value, .data$differential_FDR, .data$differential_significant,
            .data$source_path, .data$source_component), by = c("dataset_id", "metabolite")) %>% dplyr::mutate(differential_test_status = dplyr::case_when(is.na(.data$source_feature) ~
        "not_tested", .data$differential_significant %in% TRUE ~ "tested_significant", TRUE ~ "tested_not_significant"),
        selected_for_gene_analysis = dplyr::case_when(gm_settings$metabolite_scope == "significant" ~
            .data$differential_test_status == "tested_significant", gm_settings$metabolite_scope == "all_tested" ~
            .data$differential_test_status != "not_tested", gm_settings$metabolite_scope == "all_matrix" ~
            TRUE, TRUE ~ FALSE))
    selected <- reverse_audit %>% dplyr::filter(.data$selected_for_gene_analysis) %>% dplyr::pull("metabolite")
    list(source_table = source_table, reverse_audit = reverse_audit, selected = selected)
}
metabolite_selection <- purrr::map(analysis_datasets, read_metabolite_selection)
names(metabolite_selection) <- analysis_datasets
metabolite_source_audit <- dplyr::bind_rows(purrr::map(metabolite_selection, "source_table"))
metabolite_reverse_audit <- dplyr::bind_rows(purrr::map(metabolite_selection, "reverse_audit"))
readr::write_csv(metabolite_source_audit, project_output(file.path(audit_dir, "selected_metabolite_source_matching.csv"),
    name = "selected_metabolite_source_matching.csv"))
readr::write_csv(metabolite_reverse_audit, project_output(file.path(audit_dir, "selected_metabolite_reverse_audit.csv"),
    name = "selected_metabolite_reverse_audit.csv"))
selection_summary <- metabolite_reverse_audit %>% dplyr::count(.data$dataset_id, .data$differential_test_status,
    .data$selected_for_gene_analysis, name = "n_metabolites")
readr::write_csv(selection_summary, project_output(file.path(audit_dir, "selected_metabolite_summary.csv"),
    name = "selected_metabolite_summary.csv"))
if (!any(metabolite_reverse_audit$selected_for_gene_analysis)) {
    stop("No metabolites were selected for gene-level analysis under scope: ", gm_settings$metabolite_scope)
}
rna_bridge <- bundle$rna_bridge %>% dplyr::mutate(Organ = normalize_organ_name(.data$Organ), Group = factor(normalize_group(.data$Group),
    levels = c("ctrl", "ppi")))
prepare_family_alignment <- function(dataset_id, organ) {
    dataset <- bundle$metabolomics[[dataset_id]]
    dds <- dds_list[[organ]]
    selected <- metabolite_selection[[dataset_id]]$selected
    if (length(selected) == 0L) {
        return(list(valid = FALSE, reason = "no_selected_metabolites"))
    }
    metabolite_matrix <- as.matrix(dataset$matrix)
    metabolite_metadata <- as.data.frame(dataset$metadata)
    if (is.null(rownames(metabolite_matrix))) {
        stop("Metabolite matrix has no row names: ", dataset_id)
    }
    if (!"Sample_ID" %in% names(metabolite_metadata)) {
        stop("Metabolite metadata lacks Sample_ID: ", dataset_id)
    }
    if (!"Mouse_Key" %in% names(metabolite_metadata)) {
        stop("Metabolite metadata lacks Mouse_Key: ", dataset_id)
    }
    metabolite_metadata <- metabolite_metadata[match(rownames(metabolite_matrix), metabolite_metadata$Sample_ID),
        , drop = FALSE]
    if (anyNA(metabolite_metadata$Mouse_Key)) {
        stop("Missing Mouse_Key after aligning metabolite metadata for ", dataset_id)
    }
    if (anyDuplicated(metabolite_metadata$Mouse_Key)) {
        stop("Duplicated Mouse_Key in metabolite dataset: ", dataset_id)
    }
    rownames(metabolite_matrix) <- as.character(metabolite_metadata$Mouse_Key)
    keep_metabolites <- intersect(selected, colnames(metabolite_matrix))
    if (length(keep_metabolites) == 0L) {
        return(list(valid = FALSE, reason = "selected_metabolites_absent_from_matrix"))
    }
    metabolite_matrix <- metabolite_matrix[, keep_metabolites, drop = FALSE]
    bridge <- rna_bridge %>% dplyr::filter(.data$Organ == .env$organ) %>% dplyr::filter(.data$RNA_Sample %in%
        colnames(dds)) %>% dplyr::distinct(.data$RNA_Sample, .keep_all = TRUE)
    bridge <- bridge[match(colnames(dds), bridge$RNA_Sample), , drop = FALSE]
    bridge <- bridge[!is.na(bridge$Mouse_Key), , drop = FALSE]
    dds <- dds[, bridge$RNA_Sample]
    common_mouse <- intersect(as.character(bridge$Mouse_Key), rownames(metabolite_matrix))
    keep_bridge <- as.character(bridge$Mouse_Key) %in% common_mouse
    bridge <- bridge[keep_bridge, , drop = FALSE]
    dds <- dds[, bridge$RNA_Sample]
    metabolite_matrix <- metabolite_matrix[as.character(bridge$Mouse_Key), , drop = FALSE]
    rownames(metabolite_matrix) <- bridge$RNA_Sample
    model_metadata <- data.frame(RNA_Sample = as.character(bridge$RNA_Sample), Mouse_Key = as.character(bridge$Mouse_Key),
        Group = factor(normalize_group(bridge$Group), levels = c("ctrl", "ppi")), stringsAsFactors = FALSE,
        row.names = as.character(bridge$RNA_Sample))
    if (!identical(colnames(dds), model_metadata$RNA_Sample) || !identical(rownames(metabolite_matrix),
        model_metadata$RNA_Sample)) {
        stop("RNA/metabolite alignment order failed for ", dataset_id, " / ", organ)
    }
    if (anyNA(model_metadata$Group)) {
        stop("Missing Group after alignment for ", dataset_id, " / ", organ)
    }
    group_counts <- table(factor(model_metadata$Group, levels = c("ctrl", "ppi")))
    valid <- nrow(model_metadata) >= gm_settings$minimum_samples && all(group_counts >= gm_settings$minimum_samples_per_group)
    list(valid = valid, reason = if (valid) {
        "valid"
    } else {
        "insufficient_samples_or_group_size"
    }, dataset_id = dataset_id, organ = organ, dds = dds, metabolite_matrix = metabolite_matrix, metadata = model_metadata)
}
family_alignment_cache <- list()
family_alignment_audit_list <- list()
metabolite_alignment_audit_list <- list()
for (dataset_id in analysis_datasets) {
    organs <- gm_settings$tissue_map[[dataset_id]]
    for (organ in organs) {
        family_key <- paste(dataset_id, organ, sep = "__")
        aligned <- prepare_family_alignment(dataset_id, organ)
        family_alignment_cache[[family_key]] <- aligned
        if (!isTRUE(aligned$valid)) {
            family_alignment_audit_list[[family_key]] <- tibble::tibble(dataset_id = dataset_id, Organ = organ,
                valid = FALSE, reason = aligned$reason, n_samples = NA_integer_, n_ctrl = NA_integer_,
                n_ppi = NA_integer_, n_selected_metabolites = length(metabolite_selection[[dataset_id]]$selected))
            next
        }
        family_alignment_audit_list[[family_key]] <- tibble::tibble(dataset_id = dataset_id, Organ = organ,
            valid = TRUE, reason = "valid", n_samples = nrow(aligned$metadata), n_ctrl = sum(aligned$metadata$Group ==
                "ctrl"), n_ppi = sum(aligned$metadata$Group == "ppi"), n_selected_metabolites = ncol(aligned$metabolite_matrix))
        metabolite_alignment_audit_list[[family_key]] <- purrr::map_dfr(colnames(aligned$metabolite_matrix),
            function(metabolite) {
                values <- aligned$metabolite_matrix[, metabolite]
                finite <- is.finite(values)
                metadata_finite <- aligned$metadata[finite, , drop = FALSE]
                group_numeric <- as.numeric(metadata_finite$Group == "ppi")
                point_biserial_r <- if (sum(finite) >= 3L && stats::sd(values[finite], na.rm = TRUE) >
                  gm_settings$minimum_metabolite_sd) {
                  suppressWarnings(stats::cor(values[finite], group_numeric, method = "pearson", use = "complete.obs"))
                }
                else {
                  NA_real_
                }
                model_frame <- data.frame(Group = metadata_finite$Group, Metabolite = safe_zscore(values[finite]))
                tibble::tibble(dataset_id = dataset_id, Organ = organ, metabolite = metabolite, n_samples = sum(finite),
                  n_ctrl = sum(metadata_finite$Group == "ctrl"), n_ppi = sum(metadata_finite$Group ==
                    "ppi"), metabolite_sd = stats::sd(values[finite], na.rm = TRUE), group_metabolite_point_biserial_r = point_biserial_r,
                  design_condition_number = design_condition_number(~Group + Metabolite, model_frame))
            })
    }
}
family_alignment_audit <- dplyr::bind_rows(family_alignment_audit_list)
metabolite_alignment_audit <- dplyr::bind_rows(metabolite_alignment_audit_list)
readr::write_csv(family_alignment_audit, project_output(file.path(audit_dir, "family_alignment_audit.csv"),
    name = "family_alignment_audit.csv"))
readr::write_csv(metabolite_alignment_audit, project_output(file.path(audit_dir, "metabolite_alignment_audit.csv"),
    name = "metabolite_alignment_audit.csv"))
sample_name_alignment_audit <- purrr::imap_dfr(family_alignment_cache, function(aligned, family_key) {
    if (!isTRUE(aligned$valid)) {
        return(tibble::tibble(family_key = family_key, valid_family = FALSE, expression_metadata_names_identical = NA,
            metabolite_expression_names_identical = NA))
    }
    tibble::tibble(family_key = family_key, valid_family = TRUE, expression_metadata_names_identical = identical(colnames(aligned$dds),
        as.character(aligned$metadata$RNA_Sample)), metabolite_expression_names_identical = identical(rownames(aligned$metabolite_matrix),
        colnames(aligned$dds)))
})
readr::write_csv(sample_name_alignment_audit, project_output(file.path(audit_dir, "sample_name_alignment_audit.csv"),
    name = "sample_name_alignment_audit.csv"))
if (any(sample_name_alignment_audit$valid_family & (!sample_name_alignment_audit$expression_metadata_names_identical |
    !sample_name_alignment_audit$metabolite_expression_names_identical), na.rm = TRUE)) {
    stop("At least one valid family has mismatched sample names. Inspect ", file.path(audit_dir, "sample_name_alignment_audit.csv"))
}
if (!any(family_alignment_audit$valid)) {
    stop("No valid RNA–metabolomics families remained after sample alignment.")
}
all_ensembl <- unique(strip_ensembl_version(unlist(purrr::map(dds_list, rownames))))
symbol_map_vector <- AnnotationDbi::mapIds(org.Mm.eg.db, keys = all_ensembl, column = "SYMBOL", keytype = "ENSEMBL",
    multiVals = "first")
gene_map <- tibble::tibble(ENSEMBL = names(symbol_map_vector), SYMBOL = unname(as.character(symbol_map_vector))) %>%
    dplyr::filter(!is.na(.data$ENSEMBL), .data$ENSEMBL != "") %>% dplyr::distinct(.data$ENSEMBL, .keep_all = TRUE)
readr::write_csv(gene_map, project_output(file.path(audit_dir, "ensembl_to_symbol_display_mapping.csv"),
    name = "ensembl_to_symbol_display_mapping.csv"))
available_msig_collections <- msigdbr::msigdbr_collections(db_species = "MM")
readr::write_csv(available_msig_collections, project_output(file.path(audit_dir, "available_mouse_msigdb_collections.csv"),
    name = "available_mouse_msigdb_collections.csv"))
collection_col <- first_existing_column(available_msig_collections, c("gs_collection", "collection"))
subcollection_col <- first_existing_column(available_msig_collections, c("gs_subcollection", "subcollection"))
hallmark_rows <- available_msig_collections %>% dplyr::filter(toupper(as.character(.data[[collection_col]])) ==
    "MH")
reactome_rows <- available_msig_collections %>% dplyr::filter(toupper(as.character(.data[[collection_col]])) ==
    "M2", toupper(as.character(.data[[subcollection_col]])) == "CP:REACTOME")
if (nrow(hallmark_rows) == 0L || nrow(reactome_rows) == 0L) {
    stop("The installed msigdbr database does not expose the required ", "mouse-native collections MH and M2:CP:REACTOME. Inspect ",
        file.path(audit_dir, "available_mouse_msigdb_collections.csv"), " and update msigdbr if necessary.")
}
mouse_hallmark_collection <- as.character(hallmark_rows[[collection_col]][[1]])
mouse_reactome_collection <- as.character(reactome_rows[[collection_col]][[1]])
mouse_reactome_subcollection <- as.character(reactome_rows[[subcollection_col]][[1]])
message("Loading mouse-native MSigDB Hallmark collection: ", mouse_hallmark_collection)
msig_hallmark <- msigdbr::msigdbr(db_species = "MM", species = "Mus musculus", collection = mouse_hallmark_collection)
message("Loading mouse-native MSigDB Reactome collection: ", mouse_reactome_collection, " / ", mouse_reactome_subcollection)
msig_reactome <- msigdbr::msigdbr(db_species = "MM", species = "Mus musculus", collection = mouse_reactome_collection,
    subcollection = mouse_reactome_subcollection)
prepare_msig_enrichment_table <- function(msig_table, library) {
    ensembl_col <- first_existing_column(msig_table, c("ensembl_gene", "db_ensembl_gene"))
    out <- msig_table %>% dplyr::transmute(Pathway = as.character(.data$gs_name), ENSEMBL = strip_ensembl_version(.data[[ensembl_col]])) %>%
        dplyr::filter(!is.na(.data$Pathway), .data$Pathway != "", !is.na(.data$ENSEMBL), .data$ENSEMBL !=
            "") %>% dplyr::distinct()
    if (nrow(out) == 0L) {
        stop("No usable mouse Ensembl identifiers were recovered for ", library)
    }
    out
}
hallmark_term2gene <- prepare_msig_enrichment_table(msig_hallmark, "HALLMARK")
reactome_term2gene <- prepare_msig_enrichment_table(msig_reactome, "REACTOME")
pathway_collections <- list(HALLMARK = list(sets = split(hallmark_term2gene$ENSEMBL, hallmark_term2gene$Pathway),
    term2gene = hallmark_term2gene %>% dplyr::transmute(term = .data$Pathway, gene = .data$ENSEMBL),
    min_size = gm_settings$hallmark_min_size, max_size = gm_settings$hallmark_max_size), REACTOME = list(sets = split(reactome_term2gene$ENSEMBL,
    reactome_term2gene$Pathway), term2gene = reactome_term2gene %>% dplyr::transmute(term = .data$Pathway,
    gene = .data$ENSEMBL), min_size = gm_settings$reactome_min_size, max_size = gm_settings$reactome_max_size))
pathway_collection_audit <- purrr::imap_dfr(pathway_collections, function(collection, library) {
    sizes <- lengths(collection$sets)
    tibble::tibble(Library = library, n_pathways_total = length(sizes), n_pathways_within_size_limits = sum(sizes >=
        collection$min_size & sizes <= collection$max_size), minimum_size = collection$min_size, maximum_size = collection$max_size,
        msigdb_version = paste(unique(c(msig_hallmark$db_version, msig_reactome$db_version)), collapse = ";"))
})
readr::write_csv(pathway_collection_audit, project_output(file.path(audit_dir, "pathway_collection_audit.csv"),
    name = "pathway_collection_audit.csv"))
add_gene_annotation <- function(gene_table) {
    gene_table %>% dplyr::mutate(ENSEMBL = strip_ensembl_version(.data$gene)) %>% dplyr::left_join(gene_map,
        by = "ENSEMBL")
}
make_ranked_gene_vector <- function(gene_table) {
    ranked <- gene_table %>% dplyr::mutate(ENSEMBL = strip_ensembl_version(.data$gene)) %>% dplyr::filter(!is.na(.data$ENSEMBL),
        .data$ENSEMBL != "", is.finite(.data$rank_stat)) %>% dplyr::group_by(.data$ENSEMBL) %>% dplyr::slice_max(order_by = abs(.data$rank_stat),
        n = 1L, with_ties = FALSE) %>% dplyr::ungroup() %>% dplyr::arrange(dplyr::desc(.data$rank_stat))
    stats <- ranked$rank_stat
    names(stats) <- ranked$ENSEMBL
    stats
}
ensembl_ids_to_symbols <- function(ids) {
    ids <- strip_ensembl_version(ids)
    symbols <- gene_map$SYMBOL[match(ids, gene_map$ENSEMBL)]
    ifelse(is.na(symbols) | symbols == "", ids, symbols)
}
prepare_model_input <- function(aligned, metabolite) {
    values_all <- aligned$metabolite_matrix[, metabolite]
    names(values_all) <- rownames(aligned$metabolite_matrix)
    metadata <- as.data.frame(aligned$metadata, stringsAsFactors = FALSE)
    keep <- is.finite(values_all)
    kept_samples <- names(values_all)[keep]
    if (length(kept_samples) == 0L) {
        return(list(valid = FALSE, reason = "no_finite_metabolite_samples"))
    }
    dds <- aligned$dds[, kept_samples]
    metadata <- metadata[match(colnames(dds), metadata$RNA_Sample), , drop = FALSE]
    values <- values_all[colnames(dds)]
    if (anyNA(metadata$RNA_Sample) || !identical(as.character(metadata$RNA_Sample), colnames(dds)) ||
        !identical(names(values), colnames(dds))) {
        stop("RNA/metabolite ordering failed for ", aligned$dataset_id, " / ", aligned$organ, " / ",
            metabolite)
    }
    z_values <- safe_zscore(values)
    group_counts <- table(factor(as.character(metadata$Group), levels = c("ctrl", "ppi")))
    valid <- length(z_values) >= gm_settings$minimum_samples && all(is.finite(z_values)) && all(group_counts >=
        gm_settings$minimum_samples_per_group)
    if (!valid) {
        return(list(valid = FALSE, reason = "insufficient_samples_group_size_or_constant_metabolite"))
    }
    metadata$Metabolite <- as.numeric(z_values)
    metadata$Group <- droplevels(factor(metadata$Group, levels = c("ctrl", "ppi")))
    rownames(metadata) <- colnames(dds)
    if (!identical(colnames(dds), rownames(metadata))) {
        stop("Expression and metadata sample names are not identical ", "after model-input alignment.")
    }
    formula <- ~Group + Metabolite
    if (!is_full_rank_formula(formula, metadata)) {
        return(list(valid = FALSE, reason = "design_matrix_not_full_rank"))
    }
    list(valid = TRUE, dds = dds, metadata = metadata, metabolite = metabolite, design_condition_number = design_condition_number(formula,
        metadata))
}
fit_group_metabolite_deseq2 <- function(model_input) {
    metadata <- as.data.frame(model_input$metadata, stringsAsFactors = FALSE)
    dds <- model_input$dds
    metadata <- metadata[colnames(dds), , drop = FALSE]
    rownames(metadata) <- colnames(dds)
    if (!identical(colnames(dds), rownames(metadata))) {
        stop("DESeq2 input sample names do not match metadata row names.")
    }
    SummarizedExperiment::colData(dds)$Group <- metadata$Group
    SummarizedExperiment::colData(dds)$Metabolite <- metadata$Metabolite
    formula <- ~Group + Metabolite
    if (!is_full_rank_formula(formula, as.data.frame(SummarizedExperiment::colData(dds)))) {
        stop("DESeq2 design matrix is not full rank for ", model_input$metabolite)
    }
    DESeq2::design(dds) <- formula
    dds <- DESeq2::DESeq(dds, quiet = TRUE, betaPrior = FALSE)
    result_names <- DESeq2::resultsNames(dds)
    coefficient_name <- result_names[result_names == "Metabolite"]
    if (length(coefficient_name) != 1L) {
        stop("Could not identify exactly one DESeq2 coefficient named ", "'Metabolite'. Available coefficients: ",
            paste(result_names, collapse = ", "))
    }
    result <- DESeq2::results(dds, name = coefficient_name, independentFiltering = gm_settings$independent_filtering,
        cooksCutoff = gm_settings$cooks_cutoff) %>% as.data.frame() %>% tibble::rownames_to_column("gene")
    result %>% dplyr::transmute(gene = .data$gene, effect_per_SD = as.numeric(.data$log2FoldChange),
        standard_error = as.numeric(.data$lfcSE), statistic = as.numeric(.data$stat), rank_stat = as.numeric(.data$stat),
        rank_stat_type = "DESeq2_Wald_statistic", pvalue = as.numeric(.data$pvalue), gene_FDR_within_metabolite = as.numeric(.data$padj),
        average_expression = as.numeric(.data$baseMean))
}
fit_one_model <- function(model_input) {
    tryCatch({
        gene_table <- fit_group_metabolite_deseq2(model_input)
        list(success = TRUE, gene_table = gene_table, error = NA_character_)
    }, error = function(e) {
        list(success = FALSE, gene_table = NULL, error = conditionMessage(e))
    })
}
run_gsea_for_gene_table <- function(gene_table, dataset_id, organ, metabolite) {
    stats <- make_ranked_gene_vector(gene_table)
    if (length(stats) < 100L) {
        return(tibble::tibble())
    }
    purrr::imap_dfr(pathway_collections, function(collection, library) {
        fgsea_result <- tryCatch(fgsea::fgseaMultilevel(pathways = collection$sets, stats = stats, minSize = collection$min_size,
            maxSize = collection$max_size, eps = 0), error = function(e) {
            warning("GSEA failed for ", dataset_id, " / ", organ, " / ", metabolite, " / ", library,
                ": ", conditionMessage(e))
            NULL
        })
        if (is.null(fgsea_result) || nrow(fgsea_result) == 0L) {
            return(tibble::tibble())
        }
        tibble::as_tibble(fgsea_result) %>% dplyr::transmute(dataset_id = dataset_id, Organ = organ,
            metabolite = metabolite, model = "counts ~ Group + Metabolite", Library = library, Pathway = as.character(.data$pathway),
            ES = as.numeric(.data$ES), NES = as.numeric(.data$NES), pvalue = as.numeric(.data$pval),
            GSEA_FDR_within_metabolite_library = as.numeric(.data$padj), pathway_size = as.integer(.data$size),
            leading_edge_ensembl = purrr::map_chr(.data$leadingEdge, ~paste(as.character(.x), collapse = ";")),
            leading_edge_symbols = purrr::map_chr(.data$leadingEdge, ~paste(ensembl_ids_to_symbols(as.character(.x)),
                collapse = ";")), n_ranked_genes = length(stats))
    })
}
run_ora_for_gene_table <- function(gene_table, dataset_id, organ, metabolite) {
    annotated <- gene_table %>% dplyr::mutate(ENSEMBL = strip_ensembl_version(.data$gene)) %>% dplyr::filter(!is.na(.data$ENSEMBL),
        .data$ENSEMBL != "") %>% dplyr::group_by(.data$ENSEMBL) %>% dplyr::slice_min(order_by = dplyr::coalesce(.data$gene_FDR_within_metabolite,
        Inf), n = 1L, with_ties = FALSE) %>% dplyr::ungroup()
    universe <- unique(annotated$ENSEMBL)
    if (length(universe) < 100L) {
        return(tibble::tibble())
    }
    purrr::map_dfr(c("positive", "negative"), function(direction) {
        selected_table <- annotated %>% dplyr::filter(is.finite(.data$gene_FDR_within_metabolite), .data$gene_FDR_within_metabolite <
            gm_settings$ora_gene_fdr_cutoff)
        if (direction == "positive") {
            selected_table <- selected_table %>% dplyr::filter(.data$effect_per_SD > 0)
        }
        else {
            selected_table <- selected_table %>% dplyr::filter(.data$effect_per_SD < 0)
        }
        selected <- selected_table %>% dplyr::pull("ENSEMBL") %>% unique()
        if (length(selected) < gm_settings$ora_min_genes_per_direction) {
            return(tibble::tibble())
        }
        purrr::imap_dfr(pathway_collections, function(collection, library) {
            enrichment <- tryCatch(clusterProfiler::enricher(gene = selected, universe = universe, TERM2GENE = collection$term2gene,
                pvalueCutoff = 1, qvalueCutoff = 1, pAdjustMethod = "BH", minGSSize = collection$min_size,
                maxGSSize = collection$max_size), error = function(e) {
                warning("ORA failed for ", dataset_id, " / ", organ, " / ", metabolite, " / ", direction,
                  " / ", library, ": ", conditionMessage(e))
                NULL
            })
            if (is.null(enrichment)) {
                return(tibble::tibble())
            }
            enrichment_table <- as.data.frame(enrichment)
            if (nrow(enrichment_table) == 0L) {
                return(tibble::tibble())
            }
            tibble::as_tibble(enrichment_table) %>% dplyr::mutate(gene_ensembl = as.character(.data$geneID),
                gene_symbols = purrr::map_chr(strsplit(as.character(.data$geneID), "/"), ~paste(ensembl_ids_to_symbols(.x),
                  collapse = "/"))) %>% dplyr::transmute(dataset_id = dataset_id, Organ = organ, metabolite = metabolite,
                model = "counts ~ Group + Metabolite", direction = direction, Library = library, Pathway = as.character(.data$ID),
                Description = as.character(.data$Description), GeneRatio = as.character(.data$GeneRatio),
                BgRatio = as.character(.data$BgRatio), pvalue = as.numeric(.data$pvalue), ORA_FDR_within_metabolite_library_direction = as.numeric(.data$p.adjust),
                qvalue = as.numeric(.data$qvalue), gene_ensembl = .data$gene_ensembl, gene_symbols = .data$gene_symbols,
                Count = as.integer(.data$Count), n_selected_genes = length(selected), n_universe_genes = length(universe))
        })
    })
}
model_audit_list <- list()
gene_summary_list <- list()
significant_gene_list <- list()
all_gsea_list <- list()
all_ora_list <- list()
for (family_key in names(family_alignment_cache)) {
    aligned <- family_alignment_cache[[family_key]]
    if (!isTRUE(aligned$valid))
        next
    dataset_id <- aligned$dataset_id
    organ <- aligned$organ
    metabolites <- colnames(aligned$metabolite_matrix)
    message("Running counts ~ Group + Metabolite for ", dataset_id, " / ", organ, " (", length(metabolites),
        " metabolites)")
    family_gene_tables <- list()
    family_gsea_tables <- list()
    family_ora_tables <- list()
    for (metabolite in metabolites) {
        message("  ", metabolite)
        model_input <- prepare_model_input(aligned, metabolite)
        if (!isTRUE(model_input$valid)) {
            model_audit_list[[paste(family_key, metabolite, "input", sep = "__")]] <- tibble::tibble(dataset_id = dataset_id,
                Organ = organ, metabolite = metabolite, model = "counts ~ Group + Metabolite", success = FALSE,
                n_samples = NA_integer_, n_ctrl = NA_integer_, n_ppi = NA_integer_, design_condition_number = NA_real_,
                error = model_input$reason)
            next
        }
        fit <- fit_one_model(model_input)
        model_audit_list[[paste(family_key, metabolite, sep = "__")]] <- tibble::tibble(dataset_id = dataset_id,
            Organ = organ, metabolite = metabolite, model = "counts ~ Group + Metabolite", success = fit$success,
            n_samples = nrow(model_input$metadata), n_ctrl = sum(model_input$metadata$Group == "ctrl"),
            n_ppi = sum(model_input$metadata$Group == "ppi"), design_condition_number = model_input$design_condition_number,
            error = fit$error)
        if (!isTRUE(fit$success))
            next
        gene_table <- fit$gene_table %>% dplyr::mutate(dataset_id = dataset_id, Organ = organ, metabolite = metabolite,
            model = "counts ~ Group + Metabolite", n_samples = nrow(model_input$metadata), n_ctrl = sum(model_input$metadata$Group ==
                "ctrl"), n_ppi = sum(model_input$metadata$Group == "ppi"), .before = 1)
        family_gene_tables[[metabolite]] <- gene_table
        family_gsea_tables[[metabolite]] <- run_gsea_for_gene_table(gene_table, dataset_id, organ, metabolite)
        family_ora_tables[[metabolite]] <- run_ora_for_gene_table(gene_table, dataset_id, organ, metabolite)
    }
    gene_family <- dplyr::bind_rows(family_gene_tables)
    if (nrow(gene_family) > 0L) {
        gene_family <- gene_family %>% dplyr::mutate(gene_FDR_family = p.adjust(.data$pvalue, method = "BH"))
        readr::write_csv(gene_family, project_output(file.path(gene_result_dir, paste0("gene_results__",
            safe_file_id(dataset_id), "__", safe_file_id(organ), ".csv.gz"))))
        significant_gene_rows <- gene_family %>% dplyr::filter(.data$gene_FDR_within_metabolite < gm_settings$gene_fdr_cutoff |
            .data$gene_FDR_family < gm_settings$gene_fdr_cutoff)
        if (nrow(significant_gene_rows) > 0L) {
            significant_gene_list[[family_key]] <- significant_gene_rows
        }
        gene_summary_list[[family_key]] <- gene_family %>% dplyr::group_by(.data$dataset_id, .data$Organ,
            .data$metabolite, .data$model, .data$n_samples, .data$n_ctrl, .data$n_ppi) %>% dplyr::summarise(n_genes_tested = sum(is.finite(.data$pvalue)),
            n_gene_FDR_within_metabolite_0p1 = sum(.data$gene_FDR_within_metabolite < gm_settings$gene_fdr_cutoff,
                na.rm = TRUE), n_gene_FDR_family_0p1 = sum(.data$gene_FDR_family < gm_settings$gene_fdr_cutoff,
                na.rm = TRUE), .groups = "drop")
    }
    gsea_family <- dplyr::bind_rows(family_gsea_tables)
    if (nrow(gsea_family) > 0L) {
        all_gsea_list[[family_key]] <- gsea_family
        readr::write_csv(gsea_family, project_output(file.path(gsea_result_dir, paste0("GSEA__", safe_file_id(family_key),
            ".csv.gz"))))
    }
    ora_family <- dplyr::bind_rows(family_ora_tables)
    if (nrow(ora_family) > 0L) {
        all_ora_list[[family_key]] <- ora_family
        readr::write_csv(ora_family, project_output(file.path(ora_result_dir, paste0("ORA__", safe_file_id(family_key),
            ".csv.gz"))))
    }
}
model_audit <- dplyr::bind_rows(model_audit_list)
gene_summary <- dplyr::bind_rows(gene_summary_list)
significant_genes <- dplyr::bind_rows(significant_gene_list)
gsea_results <- dplyr::bind_rows(all_gsea_list)
ora_results <- dplyr::bind_rows(all_ora_list)
readr::write_csv(model_audit, project_output(file.path(audit_dir, "gene_model_audit.csv"), name = "gene_model_audit.csv"))
readr::write_csv(gene_summary, project_output(file.path(output_root, "gene_model_summary.csv"), name = "gene_model_summary.csv"))
if (nrow(significant_genes) > 0L) {
    readr::write_csv(significant_genes, project_output(file.path(output_root, "significant_gene_metabolite_associations.csv.gz")))
}
if (nrow(model_audit) == 0L || !any(model_audit$success %in% TRUE)) {
    stop("All gene–metabolite models failed. Inspect ", file.path(audit_dir, "gene_model_audit.csv"))
}
if (nrow(gsea_results) > 0L) {
    gsea_results <- gsea_results %>% dplyr::group_by(.data$dataset_id, .data$Organ, .data$Library) %>%
        dplyr::mutate(GSEA_FDR_family = p.adjust(.data$pvalue, method = "BH")) %>% dplyr::ungroup() %>%
        dplyr::group_by(.data$Library) %>% dplyr::mutate(GSEA_FDR_library_global = p.adjust(.data$pvalue,
        method = "BH")) %>% dplyr::ungroup() %>% dplyr::mutate(GSEA_FDR_global = p.adjust(.data$pvalue,
        method = "BH"))
    readr::write_csv(gsea_results, project_output(file.path(output_root, "ALL_gene_metabolite_GSEA_results.csv.gz")))
    significant_gsea <- gsea_results %>% dplyr::filter(.data$GSEA_FDR_within_metabolite_library < gm_settings$gsea_fdr_cutoff |
        .data$GSEA_FDR_family < gm_settings$gsea_fdr_cutoff) %>% dplyr::arrange(.data$GSEA_FDR_family,
        .data$pvalue)
    readr::write_csv(significant_gsea, project_output(file.path(output_root, "significant_gene_metabolite_GSEA_results.csv"),
        name = "significant_gene_metabolite_gsea_results.csv"))
} else {
    significant_gsea <- tibble::tibble()
}
if (nrow(ora_results) > 0L) {
    ora_results <- ora_results %>% dplyr::group_by(.data$dataset_id, .data$Organ, .data$Library, .data$direction) %>%
        dplyr::mutate(ORA_FDR_family = p.adjust(.data$pvalue, method = "BH")) %>% dplyr::ungroup() %>%
        dplyr::group_by(.data$Library, .data$direction) %>% dplyr::mutate(ORA_FDR_library_direction_global = p.adjust(.data$pvalue,
        method = "BH")) %>% dplyr::ungroup() %>% dplyr::mutate(ORA_FDR_global = p.adjust(.data$pvalue,
        method = "BH"))
    readr::write_csv(ora_results, project_output(file.path(output_root, "ALL_gene_metabolite_ORA_results.csv.gz")))
    significant_ora <- ora_results %>% dplyr::filter(.data$ORA_FDR_within_metabolite_library_direction <
        gm_settings$ora_fdr_cutoff | .data$ORA_FDR_family < gm_settings$ora_fdr_cutoff) %>% dplyr::arrange(.data$ORA_FDR_family,
        .data$pvalue)
    readr::write_csv(significant_ora, project_output(file.path(output_root, "significant_gene_metabolite_ORA_results.csv"),
        name = "significant_gene_metabolite_ora_results.csv"))
} else {
    significant_ora <- tibble::tibble()
}
model_completion_summary <- model_audit %>% dplyr::count(.data$success, name = "n_models")
gsea_summary <- if (nrow(gsea_results) > 0L) {
    gsea_results %>% dplyr::group_by(.data$dataset_id, .data$Organ, .data$Library) %>% dplyr::summarise(n_tests = dplyr::n(),
        n_metabolites = dplyr::n_distinct(.data$metabolite), n_within_metabolite_FDR_0p1 = sum(.data$GSEA_FDR_within_metabolite_library <
            gm_settings$gsea_fdr_cutoff, na.rm = TRUE), n_family_FDR_0p1 = sum(.data$GSEA_FDR_family <
            gm_settings$gsea_fdr_cutoff, na.rm = TRUE), n_library_global_FDR_0p1 = sum(.data$GSEA_FDR_library_global <
            gm_settings$gsea_fdr_cutoff, na.rm = TRUE), n_global_FDR_0p1 = sum(.data$GSEA_FDR_global <
            gm_settings$gsea_fdr_cutoff, na.rm = TRUE), .groups = "drop")
} else {
    tibble::tibble()
}
ora_summary <- if (nrow(ora_results) > 0L) {
    ora_results %>% dplyr::group_by(.data$dataset_id, .data$Organ, .data$Library, .data$direction) %>%
        dplyr::summarise(n_tests = dplyr::n(), n_metabolites = dplyr::n_distinct(.data$metabolite), n_within_metabolite_FDR_0p1 = sum(.data$ORA_FDR_within_metabolite_library_direction <
            gm_settings$ora_fdr_cutoff, na.rm = TRUE), n_family_FDR_0p1 = sum(.data$ORA_FDR_family <
            gm_settings$ora_fdr_cutoff, na.rm = TRUE), n_global_FDR_0p1 = sum(.data$ORA_FDR_global <
            gm_settings$ora_fdr_cutoff, na.rm = TRUE), .groups = "drop")
} else {
    tibble::tibble()
}
readr::write_csv(model_completion_summary, project_output(file.path(output_root, "model_completion_summary.csv"),
    name = "model_completion_summary.csv"))
readr::write_csv(gsea_summary, project_output(file.path(output_root, "GSEA_analysis_summary.csv"), name = "gsea_analysis_summary.csv"))
readr::write_csv(ora_summary, project_output(file.path(output_root, "ORA_analysis_summary.csv"), name = "ora_analysis_summary.csv"))
xlsx_sheets <- list(model_completion = model_completion_summary, alignment_families = family_alignment_audit,
    metabolite_selection = selection_summary, gene_summary = gene_summary, GSEA_summary = gsea_summary,
    ORA_summary = ora_summary)
if (nrow(significant_gsea) > 0L) {
    xlsx_sheets$GSEA_significant <- significant_gsea %>% dplyr::select(-dplyr::any_of(c("leading_edge_ensembl",
        "leading_edge_symbols"))) %>% dplyr::slice_head(n = 50000L)
}
if (nrow(significant_ora) > 0L) {
    xlsx_sheets$ORA_significant <- significant_ora %>% dplyr::slice_head(n = 50000L)
}
writexl::write_xlsx(xlsx_sheets, project_output(file.path(output_root, "PPI_gene_metabolite_DESeq2_GSEA_ORA_summary.xlsx"),
    name = "ppi_gene_metabolite_de_seq2_gsea_ora_summary.xlsx"))
saveRDS(list(settings = gm_settings, model_formula = ~Group + Metabolite, metabolite_selection = metabolite_selection,
    family_alignment_audit = family_alignment_audit, metabolite_alignment_audit = metabolite_alignment_audit,
    sample_name_alignment_audit = sample_name_alignment_audit, pathway_collection_audit = pathway_collection_audit,
    model_audit = model_audit, gene_summary = gene_summary, gsea_results = gsea_results, ora_results = ora_results,
    prepared_input_file = prepared_file, rna_loader_path = rna_loader_path, serum_result_file = serum_result_file,
    cecum_feces_combined_result_file = cecum_feces_combined_result_file, cecum_nmr_input = cecum_nmr_input,
    session_info = capture.output(sessionInfo())), project_output(file.path(output_root, "PPI_gene_metabolite_DESeq2_GSEA_ORA_objects.rds"),
    name = "ppi_gene_metabolite_de_seq2_gsea_ora_objects.rds"), compress = "xz")
writeLines(capture.output(sessionInfo()), project_output(file.path(output_root, "sessionInfo.txt"), name = "session_info.txt"))
message("DESeq2 gene–metabolite modelling, Hallmark/Reactome GSEA, ", "and Hallmark/Reactome ORA complete.")
message("Outputs written to: ", output_root)
}))
