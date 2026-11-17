# SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0
# MetaCardis blood-flux analysis; unreported sex-stratified/interaction branches omitted.
# Primary fitting, sex adjustment, host-model QC, clustering and enrichment expressions are unchanged.

source("scripts/helpers/setup.R")
invisible(run_analysis(
  section = "02_human_metabolism",
  directory = "human_cohorts/blood_flux/metacardis",
  script = "02_01_metacardis_blood_flux_associations.R",
  restricted = TRUE,
  code = {
library(R.matlab)
library(readxl)
library(dplyr)
library(purrr)
library(tibble)
library(tidyr)
library(stringr)
library(readr)
library(data.table)
library(matrixStats)
library(irlba)
library(ggplot2)
library(lmtest)
library(sandwich)
library(igraph)
library(grid)
library(scales)
mat_dir <- project_path("human_cohorts/blood_flux/metacardis/Blood_FBA_Metacardis/resultFlux_Metacardis_v2/output_constraintTest")
meta_path <- project_path("human_cohorts/metacardis/Merged_Metacardis-Study.xlsx")
dict_path <- project_path("human_cohorts/kiel/recon-store-metabolites-1 (1).tsv")
map_file <- project_path("human_cohorts/results/figs_PPI_Paper/figs/data/VMH_HMDB_map.csv")
hmdb_file <- project_path("human_cohorts/results/figs_PPI_Paper/figs/data/hmdb_class_annotations_collapsed.tsv")
blood_out_dir <- project_path("human_cohorts/blood_flux/metacardis/Blood_FBA_Metacardis/resultFlux_Metacardis_v2/analysis_outputs")
dir.create(project_output(blood_out_dir, directory = TRUE), showWarnings = FALSE, recursive = TRUE)
analysis_version <- "MetaCardis_Blood_FBA_v2_new_models"
high_cor_cutoff <- 0.9999
result_fdr_cutoff <- 0.05
enrichment_min_class_size <- 5
enrichment_fdr_cutoff <- 0.1
extract_mat_string <- function(x) {
    while (is.list(x) && length(x) == 1) {
        x <- x[[1]]
    }
    out <- as.character(x)
    if (length(out) == 0) {
        return(NA_character_)
    }
    out[1]
}
extract_rxns <- function(rxns_obj) {
    rxns_list <- as.vector(rxns_obj)
    vapply(rxns_list, extract_mat_string, character(1))
}
extract_one_fba_mat <- function(file) {
    mat_object <- R.matlab::readMat(project_input(file))
    if (!"rxns" %in% names(mat_object)) {
        stop("No 'rxns' object in: ", basename(file))
    }
    if (!"f" %in% names(mat_object)) {
        stop("No 'f' object in: ", basename(file))
    }
    rxns <- extract_rxns(mat_object$rxns)
    flux <- as.numeric(mat_object$f)
    if (length(rxns) != length(flux)) {
        stop("Length mismatch in file: ", basename(file), "\nrxns = ", length(rxns), "\nflux = ", length(flux))
    }
    if (any(is.na(rxns) | !nzchar(rxns))) {
        stop("Empty reaction names in: ", basename(file))
    }
    if (anyDuplicated(rxns)) {
        stop("Duplicated reaction names in: ", basename(file))
    }
    mat_sampleID <- if ("ID" %in% names(mat_object)) {
        extract_mat_string(mat_object$ID)
    }
    else {
        NA_character_
    }
    mat_sex <- if ("sex" %in% names(mat_object)) {
        extract_mat_string(mat_object$sex)
    }
    else {
        NA_character_
    }
    sampleID <- stringr::str_extract(basename(file), "x[0-9]+MCx[0-9]+")
    model_sex <- stringr::str_extract(basename(file), "(?<=mWBM_)(Female|Male)")
    if (is.na(sampleID)) {
        stop("Could not extract MetaCardis ID from filename: ", basename(file))
    }
    if (is.na(model_sex)) {
        stop("Could not extract model sex from filename: ", basename(file))
    }
    out <- as.data.frame(matrix(flux, nrow = 1), check.names = FALSE)
    names(out) <- rxns
    out %>% mutate(sampleID = sampleID, mat_sampleID = mat_sampleID, model_sex = model_sex, mat_sex = mat_sex,
        file = basename(file), .before = 1)
}
mat_files <- project_list_files(mat_dir, pattern = "\\.mat$", full.names = TRUE)
message("MAT files found: ", length(mat_files))
if (length(mat_files) == 0) {
    stop("No .mat files found in: ", mat_dir)
}
file_check <- tibble(file = basename(mat_files), sampleID = stringr::str_extract(basename(mat_files),
    "x[0-9]+MCx[0-9]+"), model_sex = stringr::str_extract(basename(mat_files), "(?<=mWBM_)(Female|Male)"))
message("File-name sample ID extraction:")
print(table(is.na(file_check$sampleID), useNA = "ifany"))
message("File-name model sex extraction:")
print(table(is.na(file_check$model_sex), useNA = "ifany"))
message("MetaCardis ID family prefixes:")
print(table(stringr::str_extract(file_check$sampleID, "^x[0-9]+MCx"), useNA = "ifany"))
message("Model sex in file names:")
print(table(file_check$model_sex, useNA = "ifany"))
if (any(is.na(file_check$sampleID)) || any(is.na(file_check$model_sex))) {
    write.csv(file_check, project_output(file.path(blood_out_dir, "MetaCardis_Blood_FBA_v2_file_name_check.csv"),
        name = "meta_cardis_blood_fba_file_name_check.csv"), row.names = FALSE)
    stop("Some MAT filenames could not be parsed. See MetaCardis_Blood_FBA_v2_file_name_check.csv.")
}
safe_extract_one_fba_mat <- purrr::safely(extract_one_fba_mat)
mat_read_results <- purrr::map(mat_files, safe_extract_one_fba_mat)
mat_read_errors <- tibble(file = basename(mat_files), error = vapply(mat_read_results, function(x) {
    if (is.null(x$error)) {
        NA_character_
    }
    else {
        conditionMessage(x$error)
    }
}, character(1))) %>% filter(!is.na(error))
message("MAT files that could not be read: ", nrow(mat_read_errors))
if (nrow(mat_read_errors) > 0) {
    print(mat_read_errors, n = Inf)
    write.csv(mat_read_errors, project_output(file.path(blood_out_dir, "MetaCardis_Blood_FBA_v2_MAT_read_errors.csv"),
        name = "meta_cardis_blood_fba_mat_read_errors.csv"), row.names = FALSE)
}
metac_blood_fba_flux_wide <- mat_read_results %>% purrr::map("result") %>% purrr::compact() %>% bind_rows()
if (nrow(metac_blood_fba_flux_wide) == 0) {
    stop("No MetaCardis Blood FBA MAT files could be loaded.")
}
message("Loaded MetaCardis Blood FBA dimensions:")
print(dim(metac_blood_fba_flux_wide))
print(metac_blood_fba_flux_wide[1:min(5, nrow(metac_blood_fba_flux_wide)), 1:min(10, ncol(metac_blood_fba_flux_wide))])
message("Model sex from filenames:")
print(table(metac_blood_fba_flux_wide$model_sex, useNA = "ifany"))
message("Duplicated sample IDs: ", sum(duplicated(metac_blood_fba_flux_wide$sampleID)))
if (anyDuplicated(metac_blood_fba_flux_wide$sampleID)) {
    duplicated_ids <- metac_blood_fba_flux_wide$sampleID[duplicated(metac_blood_fba_flux_wide$sampleID)]
    stop("Duplicated sample IDs found: ", paste(unique(duplicated_ids), collapse = ", "))
}
rxn_cols_raw <- grep("^DM_", names(metac_blood_fba_flux_wide), value = TRUE)
message("Raw DM reaction columns: ", length(rxn_cols_raw))
if (length(rxn_cols_raw) == 0) {
    stop("No reaction columns beginning with 'DM_' were found.")
}
metac_blood_fba_flux_wide$n_missing <- rowSums(is.na(metac_blood_fba_flux_wide[, rxn_cols_raw, drop = FALSE]))
metac_blood_fba_flux_wide$all_nan <- metac_blood_fba_flux_wide$n_missing == length(rxn_cols_raw)
message("All-NaN samples:")
print(table(metac_blood_fba_flux_wide$all_nan, useNA = "ifany"))
metac_blood_fba <- metac_blood_fba_flux_wide %>% filter(!all_nan)
message("Samples before removing all-NaN: ", nrow(metac_blood_fba_flux_wide))
message("Samples after removing all-NaN: ", nrow(metac_blood_fba))
cardis_meta <- readxl::read_excel(project_input(meta_path)) %>% as.data.frame()
stopifnot("ID" %in% names(cardis_meta))
cardis_meta$ID <- trimws(as.character(cardis_meta$ID))
if (anyDuplicated(cardis_meta$ID)) {
    stop("Duplicated IDs are present in the MetaCardis metadata.")
}
convert_yes_no_numeric <- function(x) {
    x <- trimws(tolower(as.character(x)))
    dplyr::case_when(x %in% c("yes", "1", "true") ~ 1, x %in% c("no", "0", "false") ~ 0, TRUE ~ NA_real_)
}
normalize_yes_no_factor <- function(x) {
    x <- trimws(tolower(as.character(x)))
    out <- dplyr::case_when(x %in% c("yes", "1", "true") ~ "Yes", x %in% c("no", "0", "false") ~ "No",
        TRUE ~ NA_character_)
    factor(out, levels = c("No", "Yes"))
}
med_vars <- c("Antidiabetic.drugs", "Insulin", "Metformin", "Sulfonurea", "DPP.IV.inhibitors", "GLP1.agonists",
    "Anti-hypertensive.drugs", "Beta-blockers", "ACE", "AT2.inhibitor", "Calcium.blockers", "Anse.diuretics",
    "Potassium-.sparing.diuretic", "Thiazides", "Lipid.lowering.drugs", "Statin", "Ezetimib", "Fibrate",
    "Nitrate", "Amiodarone", "Digitalic", "Aspirine", "Clopidogrel", "Coumarin")
med_vars <- intersect(med_vars, names(cardis_meta))
cardis_meta[med_vars] <- lapply(cardis_meta[med_vars], convert_yes_no_numeric)
med_broad <- intersect(c("Antidiabetic.drugs", "Anti-hypertensive.drugs", "Lipid.lowering.drugs", "Nitrate",
    "Amiodarone", "Digitalic", "Aspirine", "Clopidogrel", "Coumarin"), names(cardis_meta))
cardis_meta$medication_count <- rowSums(cardis_meta[, med_broad, drop = FALSE], na.rm = TRUE)
cardis_meta$Proton.pump.inhibitor <- normalize_yes_no_factor(cardis_meta$Proton.pump.inhibitor)
cardis_meta$`Age.(years)` <- as.numeric(cardis_meta$`Age.(years)`)
cardis_meta$`BMI.(kg/m²)` <- as.numeric(cardis_meta$`BMI.(kg/m²)`)
cardis_meta$creatinine <- as.numeric(cardis_meta$creatinine)
cardis_meta$medication_count <- as.numeric(cardis_meta$medication_count)
cardis_meta$Gender <- factor(cardis_meta$Gender)
cardis_meta$Nationality <- factor(cardis_meta$Nationality)
cardis_meta$Diabetic.status <- factor(cardis_meta$Diabetic.status)
for (v in intersect(c("Metformin", "Statin", "Aspirine"), names(cardis_meta))) {
    cardis_meta[[v]] <- factor(cardis_meta[[v]], levels = c(0, 1), labels = c("No", "Yes"))
}
matched_idx <- match(metac_blood_fba$sampleID, cardis_meta$ID)
message("Blood FBA samples: ", nrow(metac_blood_fba))
message("Matched metadata samples: ", sum(!is.na(matched_idx)))
message("Missing from metadata: ", sum(is.na(matched_idx)))
missing_from_metadata <- metac_blood_fba$sampleID[is.na(matched_idx)]
if (length(missing_from_metadata) > 0) {
    message("Unmatched IDs: ", paste(missing_from_metadata, collapse = ", "))
}
keep_matched <- !is.na(matched_idx)
metac_blood_fba <- metac_blood_fba[keep_matched, , drop = FALSE]
metac_blood_meta <- cardis_meta[matched_idx[keep_matched], , drop = FALSE]
rownames(metac_blood_fba) <- NULL
rownames(metac_blood_meta) <- NULL
stopifnot(all(metac_blood_meta$ID == metac_blood_fba$sampleID))
message("Final aligned sample number: ", nrow(metac_blood_meta))
message("PPI counts:")
print(table(metac_blood_meta$Proton.pump.inhibitor, useNA = "ifany"))
message("Participant Gender counts:")
print(table(metac_blood_meta$Gender, useNA = "ifany"))
message("PPI by participant Gender:")
print(table(metac_blood_meta$Proton.pump.inhibitor, metac_blood_meta$Gender, useNA = "ifany"))
normalize_sex_label <- function(x) {
    x <- trimws(tolower(as.character(x)))
    dplyr::case_when(x %in% c("female", "f", "woman", "2") ~ "Female", x %in% c("male", "m", "man", "1") ~
        "Male", TRUE ~ NA_character_)
}
model_sex_check <- tibble(sampleID = metac_blood_fba$sampleID, participant_gender_raw = as.character(metac_blood_meta$Gender),
    participant_gender = normalize_sex_label(metac_blood_meta$Gender), model_sex_raw = as.character(metac_blood_fba$model_sex),
    model_sex = normalize_sex_label(metac_blood_fba$model_sex), mat_sex_raw = as.character(metac_blood_fba$mat_sex),
    mat_sex = normalize_sex_label(metac_blood_fba$mat_sex), file = metac_blood_fba$file) %>% mutate(participant_vs_model_match = participant_gender ==
    model_sex, participant_vs_mat_match = if_else(is.na(mat_sex), NA, participant_gender == mat_sex))
message("Participant Gender versus WBM model sex:")
print(table(model_sex_check$participant_gender, model_sex_check$model_sex, useNA = "ifany"))
message("Participant Gender versus MAT stored sex:")
print(table(model_sex_check$participant_gender, model_sex_check$mat_sex, useNA = "ifany"))
message("Participant/model sex match:")
print(table(model_sex_check$participant_vs_model_match, useNA = "ifany"))
write.csv(model_sex_check, project_output(file.path(blood_out_dir, "MetaCardis_Blood_FBA_v2_model_sex_check.csv"),
    name = "meta_cardis_blood_fba_model_sex_check.csv"), row.names = FALSE)
if (any(model_sex_check$participant_vs_model_match == FALSE, na.rm = TRUE)) {
    warning("Some MetaCardis v2 WBM model-sex labels do not match participant Gender.")
}
rxn_cols <- grep("^DM_", names(metac_blood_fba), value = TRUE)
mat_raw <- as.matrix(metac_blood_fba[, rxn_cols, drop = FALSE])
storage.mode(mat_raw) <- "numeric"
message("Raw reaction matrix dimensions:")
print(dim(mat_raw))
message("Raw flux summary:")
print(summary(as.vector(mat_raw)))
col_missing <- colMeans(is.na(mat_raw))
message("Reactions with any missing values: ", sum(col_missing > 0))
message("Reactions with >10% missing values: ", sum(col_missing > 0.1))
keep_complete <- col_missing == 0
mat <- mat_raw[, keep_complete, drop = FALSE]
vars <- matrixStats::colVars(mat, na.rm = TRUE)
keep_var <- is.finite(vars) & vars > 0
mat <- mat[, keep_var, drop = FALSE]
message("Reactions after missing-value and variance filters: ", ncol(mat))
if (ncol(mat) < 2) {
    stop("Too few usable reactions remain.")
}
mat_log <- sign(mat) * log1p(abs(mat))
stopifnot(nrow(mat_log) == nrow(metac_blood_meta))
stopifnot(all(metac_blood_meta$ID == metac_blood_fba$sampleID))
flux_summary <- data.frame(statistic = c("min", "q25", "median", "mean", "q75", "max"), value = c(min(mat,
    na.rm = TRUE), quantile(mat, 0.25, na.rm = TRUE), median(mat, na.rm = TRUE), mean(mat, na.rm = TRUE),
    quantile(mat, 0.75, na.rm = TRUE), max(mat, na.rm = TRUE)))
print(flux_summary)
message("Absolute raw-flux quantiles:")
print(quantile(abs(as.vector(mat)), probs = c(0.5, 0.9, 0.99, 0.999), na.rm = TRUE))
message("Absolute signed-log flux quantiles:")
print(quantile(abs(as.vector(mat_log)), probs = c(0.5, 0.9, 0.99, 0.999), na.rm = TRUE))
hist(as.vector(mat), breaks = 100, main = "MetaCardis Blood FBA raw flux distribution", xlab = "Flux")
hist(as.vector(mat_log), breaks = 100, main = "MetaCardis Blood FBA signed log1p flux distribution",
    xlab = "signed log1p(flux)")
n_pcs <- min(10, nrow(mat_log) - 1, ncol(mat_log) - 1)
if (n_pcs < 2) {
    stop("Too few dimensions for PCA.")
}
metac_blood_pca_unscaled <- irlba::prcomp_irlba(mat_log, n = n_pcs, center = TRUE, scale. = FALSE)
metac_blood_pca_scaled <- irlba::prcomp_irlba(mat_log, n = n_pcs, center = TRUE, scale. = TRUE)
total_variance_unscaled <- sum(matrixStats::colVars(mat_log), na.rm = TRUE)
total_variance_scaled <- ncol(mat_log)
var_unscaled <- metac_blood_pca_unscaled$sdev^2/total_variance_unscaled
var_scaled <- metac_blood_pca_scaled$sdev^2/total_variance_scaled
metac_blood_pca_var_summary <- bind_rows(data.frame(model = "MetaCardis_Blood_FBA_v2_unscaled", PC = paste0("PC",
    seq_along(var_unscaled)), variance_explained = var_unscaled, percent_variance_explained = 100 * var_unscaled,
    cumulative_percent = 100 * cumsum(var_unscaled)), data.frame(model = "MetaCardis_Blood_FBA_v2_scaled",
    PC = paste0("PC", seq_along(var_scaled)), variance_explained = var_scaled, percent_variance_explained = 100 *
        var_scaled, cumulative_percent = 100 * cumsum(var_scaled)))
print(metac_blood_pca_var_summary)
message("MetaCardis Blood FBA unscaled PC1: ", round(100 * var_unscaled[1], 2), "%")
message("MetaCardis Blood FBA scaled PC1: ", round(100 * var_scaled[1], 2), "%")
metac_blood_scores_scaled <- as.data.frame(metac_blood_pca_scaled$x)
metac_blood_scores_scaled$sampleID <- metac_blood_fba$sampleID
metac_blood_scores_scaled <- metac_blood_scores_scaled %>% left_join(metac_blood_meta %>% select(ID,
    Proton.pump.inhibitor, `Age.(years)`, `BMI.(kg/m²)`, Gender, Nationality, Diabetic.status, medication_count),
    by = c(sampleID = "ID"))
p_pca_var <- ggplot(metac_blood_pca_var_summary, aes(x = as.numeric(gsub("PC", "", PC)), y = percent_variance_explained)) +
    geom_col() + facet_wrap(~model, ncol = 1, scales = "free_y") + theme_bw() + labs(title = "MetaCardis Blood FBA PCA variance explained",
    x = "Principal component", y = "% total variance explained")
p_pca_ppi <- ggplot(metac_blood_scores_scaled, aes(x = PC1, y = PC2, color = Proton.pump.inhibitor)) +
    geom_point(alpha = 0.75) + theme_bw() + labs(title = "MetaCardis Blood FBA PCA by PPI", x = "PC1",
    y = "PC2", color = "PPI")
p_pca_gender <- ggplot(metac_blood_scores_scaled, aes(x = PC1, y = PC2, color = Gender)) + geom_point(alpha = 0.75) +
    theme_bw() + labs(title = "MetaCardis Blood FBA PCA by participant sex", x = "PC1", y = "PC2", color = "Gender")
print(p_pca_var)
print(p_pca_ppi)
print(p_pca_gender)
metac_blood_meta_pc1 <- metac_blood_meta
metac_blood_meta_pc1$PC1_scaled <- metac_blood_pca_scaled$x[, 1]
metac_blood_meta_pc1$PC2_scaled <- metac_blood_pca_scaled$x[, 2]
pc1_fit <- lm(PC1_scaled ~ Proton.pump.inhibitor + `Age.(years)` + `BMI.(kg/m²)` + Gender + Nationality +
    Diabetic.status + medication_count, data = metac_blood_meta_pc1)
print(summary(pc1_fit))
metac_blood_loadings <- metac_blood_pca_scaled$rotation[, 1]
names(metac_blood_loadings) <- colnames(mat_log)
print(summary(metac_blood_loadings))
message("Proportion of positive PC1 loadings: ", mean(metac_blood_loadings > 0))
print(table(cut(metac_blood_loadings, breaks = c(-Inf, -0.002, 0.002, Inf))))
top_metac_blood_pos <- data.frame(reaction = names(sort(metac_blood_loadings, decreasing = TRUE))[1:min(30,
    length(metac_blood_loadings))], loading = sort(metac_blood_loadings, decreasing = TRUE)[1:min(30,
    length(metac_blood_loadings))])
top_metac_blood_neg <- data.frame(reaction = names(sort(metac_blood_loadings))[1:min(30, length(metac_blood_loadings))],
    loading = sort(metac_blood_loadings)[1:min(30, length(metac_blood_loadings))])
run_hc3_scan <- function(mat_log, metadata, formula, target_term, progress_every = 50) {
    covars <- setdiff(all.vars(formula), "met")
    missing_covars <- setdiff(covars, names(metadata))
    if (length(missing_covars) > 0) {
        stop("Missing covariates: ", paste(missing_covars, collapse = ", "))
    }
    keep <- complete.cases(metadata[, covars, drop = FALSE])
    mat_use <- mat_log[keep, , drop = FALSE]
    meta_use <- metadata[keep, , drop = FALSE]
    message("Complete cases: ", nrow(meta_use))
    message("Reactions to test: ", ncol(mat_use))
    fit_one <- function(y) {
        df <- meta_use
        df$met <- y
        if (length(unique(df$met[is.finite(df$met)])) < 2 || sd(df$met, na.rm = TRUE) == 0) {
            return(c(beta = NA_real_, se = NA_real_, t = NA_real_, p = NA_real_, n = nrow(df)))
        }
        fit <- tryCatch(lm(formula, data = df), error = function(e) NULL)
        if (is.null(fit)) {
            return(c(beta = NA_real_, se = NA_real_, t = NA_real_, p = NA_real_, n = nrow(df)))
        }
        ct <- tryCatch(lmtest::coeftest(fit, vcov. = sandwich::vcovHC(fit, type = "HC3")), error = function(e) NULL)
        if (is.null(ct) || !(target_term %in% rownames(ct))) {
            return(c(beta = NA_real_, se = NA_real_, t = NA_real_, p = NA_real_, n = nrow(df)))
        }
        c(beta = ct[target_term, "Estimate"], se = ct[target_term, "Std. Error"], t = ct[target_term,
            "t value"], p = ct[target_term, "Pr(>|t|)"], n = nrow(df))
    }
    result_list <- lapply(seq_len(ncol(mat_use)), function(i) {
        if (progress_every > 0 && i%%progress_every == 0) {
            message("Finished ", i, " / ", ncol(mat_use))
        }
        fit_one(mat_use[, i])
    })
    result_df <- as.data.frame(do.call(rbind, result_list))
    tibble(reaction = colnames(mat_use), beta = as.numeric(result_df$beta), SE = as.numeric(result_df$se),
        t_value = as.numeric(result_df$t), pvalue = as.numeric(result_df$p), n = as.numeric(result_df$n)) %>%
        mutate(FDR = p.adjust(pvalue, method = "fdr")) %>% arrange(FDR, pvalue)
}
fml_metac_blood <- met ~ Proton.pump.inhibitor + `Age.(years)` + `BMI.(kg/m²)` + Gender + Nationality +
    Diabetic.status + medication_count + Metformin + Statin + Aspirine + creatinine
ppi_term <- "Proton.pump.inhibitorYes"
metac_blood_fba_results <- run_hc3_scan(mat_log = mat_log, metadata = metac_blood_meta, formula = fml_metac_blood,
    target_term = ppi_term)
metac_blood_result_summary <- metac_blood_fba_results %>% summarise(n_reactions = n(), nominal_p_0.05 = sum(pvalue <
    0.05, na.rm = TRUE), fdr_0.05 = sum(FDR < result_fdr_cutoff, na.rm = TRUE), min_pvalue = min(pvalue,
    na.rm = TRUE), min_FDR = min(FDR, na.rm = TRUE))
print(metac_blood_result_summary)
print(metac_blood_fba_results %>% arrange(FDR, pvalue) %>% head(30))
load_metabolite_dict <- function(dict_path) {
    md <- readr::read_tsv(project_input(dict_path), show_col_types = FALSE)
    stopifnot(all(c("abbreviation", "fullName") %in% names(md)))
    md %>% select(abbreviation, fullName) %>% mutate(abbreviation = as.character(abbreviation), fullName = as.character(fullName)) %>%
        distinct(abbreviation, .keep_all = TRUE)
}
extract_blood_fba_met_abbr <- function(x) {
    x %>% as.character() %>% str_replace("\\[.*\\]$", "") %>% str_replace("^DM_", "")
}
annotate_blood_fba_results <- function(res_tbl, dict_path) {
    md <- load_metabolite_dict(dict_path)
    res_tbl %>% mutate(reaction = as.character(reaction), met_abbr = extract_blood_fba_met_abbr(reaction)) %>%
        left_join(md, by = c(met_abbr = "abbreviation")) %>% mutate(metabolite_name = coalesce(fullName,
        met_abbr), label = paste0(metabolite_name, " (", met_abbr, ")"))
}
metac_blood_fba_results_annot <- annotate_blood_fba_results(metac_blood_fba_results, dict_path)
print(metac_blood_fba_results_annot %>% select(reaction, met_abbr, metabolite_name, beta, pvalue, FDR) %>%
    arrange(FDR, pvalue) %>% head(30))
cor_mat <- cor(mat_log, use = "pairwise.complete.obs")
cor_values <- cor_mat[upper.tri(cor_mat)]
cor_values <- cor_values[is.finite(cor_values)]
cor_summary <- data.frame(metric = c("median_abs_cor", "p90_abs_cor", "p95_abs_cor", "p99_abs_cor", "p999_abs_cor",
    "prop_abs_cor_gt_0.90", "prop_abs_cor_gt_0.95", "prop_abs_cor_gt_0.99"), value = c(median(abs(cor_values),
    na.rm = TRUE), quantile(abs(cor_values), 0.9, na.rm = TRUE), quantile(abs(cor_values), 0.95, na.rm = TRUE),
    quantile(abs(cor_values), 0.99, na.rm = TRUE), quantile(abs(cor_values), 0.999, na.rm = TRUE), mean(abs(cor_values) >
        0.9, na.rm = TRUE), mean(abs(cor_values) > 0.95, na.rm = TRUE), mean(abs(cor_values) > 0.99,
        na.rm = TRUE)))
print(cor_summary)
hist(cor_values, breaks = 100, main = "MetaCardis Blood FBA pairwise correlations", xlab = "Pearson correlation")
positive_high_cor_pairs <- which(cor_mat > high_cor_cutoff & upper.tri(cor_mat), arr.ind = TRUE)
positive_high_cor_df <- data.frame(rxn1 = colnames(cor_mat)[positive_high_cor_pairs[, 1]], rxn2 = colnames(cor_mat)[positive_high_cor_pairs[,
    2]], cor = cor_mat[positive_high_cor_pairs])
negative_high_cor_pairs <- which(cor_mat < -high_cor_cutoff & upper.tri(cor_mat), arr.ind = TRUE)
negative_high_cor_df <- data.frame(rxn1 = colnames(cor_mat)[negative_high_cor_pairs[, 1]], rxn2 = colnames(cor_mat)[negative_high_cor_pairs[,
    2]], cor = cor_mat[negative_high_cor_pairs])
message("Pairs with cor > ", high_cor_cutoff, ": ", nrow(positive_high_cor_df))
message("Pairs with cor < -", high_cor_cutoff, ": ", nrow(negative_high_cor_df))
if (nrow(positive_high_cor_df) > 0) {
    cluster_graph <- igraph::graph_from_data_frame(positive_high_cor_df[, c("rxn1", "rxn2")], directed = FALSE)
    graph_clusters <- igraph::components(cluster_graph)
    cluster_df <- data.frame(reaction = names(graph_clusters$membership), cluster_id = as.integer(graph_clusters$membership))
} else {
    cluster_df <- data.frame(reaction = character(), cluster_id = integer())
}
all_reactions <- colnames(mat_log)
singleton_reactions <- setdiff(all_reactions, cluster_df$reaction)
next_cluster_id <- if (nrow(cluster_df) == 0) {
    1L
} else {
    max(cluster_df$cluster_id) + 1L
}
singleton_df <- data.frame(reaction = singleton_reactions, cluster_id = seq.int(from = next_cluster_id,
    length.out = length(singleton_reactions)))
cluster_df <- bind_rows(cluster_df, singleton_df) %>% mutate(cluster = paste0("cluster_", cluster_id))
cluster_summary <- cluster_df %>% dplyr::count(cluster, sort = TRUE, name = "n_reactions_total")
print(cluster_summary %>% head(30))
metac_blood_fba_results_clustered <- metac_blood_fba_results_annot %>% left_join(cluster_df %>% select(reaction,
    cluster, cluster_id), by = "reaction")
metac_blood_sig_cluster_summary <- metac_blood_fba_results_clustered %>% filter(FDR < result_fdr_cutoff) %>%
    group_by(cluster) %>% arrange(pvalue, .by_group = TRUE) %>% summarise(n_sig_reactions = n(), representative_reaction = first(reaction),
    representative_metabolite = first(metabolite_name), min_FDR = min(FDR, na.rm = TRUE), min_pvalue = min(pvalue,
        na.rm = TRUE), mean_beta = mean(beta, na.rm = TRUE), median_beta = median(beta, na.rm = TRUE),
    reactions = paste(reaction, collapse = "; "), metabolites = paste(metabolite_name, collapse = "; "),
    .groups = "drop") %>% arrange(min_FDR, min_pvalue)
print(metac_blood_sig_cluster_summary)
message("Raw significant reactions: ", sum(metac_blood_fba_results$FDR < result_fdr_cutoff, na.rm = TRUE))
message("Collapsed significant clusters: ", nrow(metac_blood_sig_cluster_summary))
duplicate_beta_p_summary <- metac_blood_fba_results %>% mutate(beta_round = round(beta, 6), p_round = round(pvalue,
    6)) %>% dplyr::count(beta_round, p_round, sort = TRUE)
print(duplicate_beta_p_summary %>% head(30))
plot_cluster_boxplot <- function(cluster_name, cluster_df, flux_data, metadata, max_reactions = 12, log_transform = TRUE) {
    rxns_to_plot <- cluster_df %>% filter(cluster == cluster_name) %>% pull(reaction)
    rxns_to_plot <- intersect(rxns_to_plot, names(flux_data))
    if (length(rxns_to_plot) == 0) {
        stop("No reactions found for cluster: ", cluster_name)
    }
    rxns_to_plot <- rxns_to_plot[seq_len(min(max_reactions, length(rxns_to_plot)))]
    plot_df <- flux_data %>% select(sampleID, all_of(rxns_to_plot)) %>% left_join(metadata %>% select(ID,
        Proton.pump.inhibitor, Gender), by = c(sampleID = "ID")) %>% pivot_longer(cols = all_of(rxns_to_plot),
        names_to = "reaction", values_to = "flux") %>% mutate(value = if (log_transform) {
        sign(flux) * log1p(abs(flux))
    }
    else {
        flux
    })
    ggplot(plot_df, aes(x = Proton.pump.inhibitor, y = value, fill = Proton.pump.inhibitor)) + geom_boxplot(outlier.alpha = 0.4) +
        facet_wrap(~reaction, scales = "free_y") + theme_bw() + labs(title = paste0("MetaCardis Blood FBA: ",
        cluster_name), x = "PPI", y = ifelse(log_transform, "signed log1p(flux)", "flux")) + theme(legend.position = "none")
}
significant_cluster_names <- metac_blood_sig_cluster_summary$cluster
cluster_plots <- lapply(significant_cluster_names, function(cluster_name) {
    plot_cluster_boxplot(cluster_name = cluster_name, cluster_df = cluster_df, flux_data = metac_blood_fba,
        metadata = metac_blood_meta, max_reactions = 12, log_transform = TRUE)
})
names(cluster_plots) <- significant_cluster_names
if (length(cluster_plots) > 0) {
    print(cluster_plots[[1]])
}
pad_hmdb <- function(x) {
    x <- as.character(x)
    x <- trimws(x)
    x[x == ""] <- NA_character_
    x <- stringr::str_replace(x, "^HMDB", "")
    x <- stringr::str_pad(x, width = 7, side = "left", pad = "0")
    paste0("HMDB", x)
}
build_class_membership <- function(hmdb_annotation) {
    hmdb_annotation <- data.table::as.data.table(hmdb_annotation)
    hmdb_col <- intersect(c("hmdb", "HMDB"), names(hmdb_annotation))[1]
    if (is.na(hmdb_col)) {
        stop("Could not identify HMDB column.")
    }
    class_cols <- intersect(c("kingdom", "super_class", "superclass", "class", "sub_class", "subclass",
        "direct_parent", "class_label"), names(hmdb_annotation))
    class_cols <- setdiff(class_cols, hmdb_col)
    if (length(class_cols) == 0) {
        stop("Could not identify HMDB class columns.")
    }
    hmdb_long_raw <- data.table::melt(hmdb_annotation, id.vars = hmdb_col, measure.vars = class_cols,
        variable.name = "source_column", value.name = "raw_class_label")
    data.table::setnames(hmdb_long_raw, hmdb_col, "hmdb")
    hmdb_long_raw <- hmdb_long_raw[!is.na(hmdb) & hmdb != "" & !is.na(raw_class_label) & raw_class_label !=
        ""]
    hmdb_long <- hmdb_long_raw[, {
        parts <- unlist(strsplit(as.character(raw_class_label), ";", fixed = TRUE))
        parts <- trimws(parts)
        parts <- parts[parts != ""]
        if (length(parts) == 0) {
            NULL
        }
        else {
            data.table::data.table(class_depth = seq_along(parts), class_label = parts)
        }
    }, by = .(hmdb, source_column)]
    hmdb_long[, `:=`(class_level, paste0(as.character(source_column), "_level_", class_depth))]
    hmdb_long <- hmdb_long[!is.na(class_label) & class_label != ""]
    hmdb_long <- hmdb_long[!class_label %in% c("Organic compounds", "Inorganic compounds")]
    unique(hmdb_long[, .(hmdb, class_level, class_label)])
}
performBloodFBAEnrichmentPlot <- function(res_tbl, map_file, hmdb_file, min_class_size = 5, fdr_cutoff = 0.1,
    effect_col = "beta") {
    excr <- data.table::as.data.table(res_tbl)
    if (!"metabr" %in% names(excr)) {
        if (!"met_abbr" %in% names(excr)) {
            excr[, `:=`(met_abbr, extract_blood_fba_met_abbr(reaction))]
        }
        excr[, `:=`(metabr, trimws(as.character(met_abbr)))]
    }
    else {
        excr[, `:=`(metabr, trimws(as.character(metabr)))]
    }
    if (!(effect_col %in% names(excr))) {
        stop("effect_col not found in res_tbl: ", effect_col)
    }
    excr[, `:=`(beta_meta, as.numeric(get(effect_col)))]
    map <- data.table::fread(project_input(map_file), drop = 1, header = TRUE)
    if (!all(c("met_id", "HMDB") %in% names(map))) {
        stop("Map file did not read as expected. Columns are: ", paste(names(map), collapse = ", "))
    }
    map <- map[HMDB != "HMDB0000216"]
    hmdb_annotation <- data.table::fread(project_input(hmdb_file))
    excr[, `:=`(metabr, trimws(metabr))]
    map[, `:=`(met_id, trimws(met_id))]
    map[, `:=`(HMDB, trimws(HMDB))]
    map[HMDB == "", `:=`(HMDB, NA_character_)]
    map <- unique(map, by = c("met_id", "HMDB"))
    data.table::setnames(map, "met_id", "metabr")
    excr <- merge(excr, map, by = "metabr", all.x = TRUE, allow.cartesian = TRUE)
    excr <- excr[!is.na(HMDB)]
    excr[, `:=`(HMDB, pad_hmdb(HMDB))]
    message("Input rows: ", nrow(res_tbl))
    message("Mapped effect-HMDB rows: ", nrow(excr))
    message("Unique mapped metabolites: ", data.table::uniqueN(excr$metabr))
    message("Unique mapped HMDB IDs: ", data.table::uniqueN(excr$HMDB))
    hmdb_long <- build_class_membership(hmdb_annotation)
    hmdb_long[, `:=`(hmdb, pad_hmdb(hmdb))]
    hmdb_long <- hmdb_long[hmdb %in% excr[, HMDB]]
    testClasses <- table(hmdb_long$class_label)
    testClasses <- testClasses[testClasses >= min_class_size]
    if (length(testClasses) == 0) {
        warning("No HMDB classes passed min_class_size.")
        return(list(plot = NULL, enrichment_results = data.table::data.table(), enrichment_results_all = data.table::data.table(),
            mapped_results = excr))
    }
    res <- data.table::rbindlist(lapply(names(testClasses), function(testClass) {
        hmdbs <- unique(hmdb_long[class_label == testClass, hmdb])
        thisClass <- excr[HMDB %in% hmdbs, beta_meta]
        otherClass <- excr[!(HMDB %in% hmdbs), beta_meta]
        if (length(thisClass) < 2 || length(otherClass) < 2) {
            return(NULL)
        }
        data.table::data.table(class = testClass, nThis = length(thisClass), nOther = length(otherClass),
            meanThis = mean(thisClass, na.rm = TRUE), meanOther = mean(otherClass, na.rm = TRUE), medianThis = median(thisClass,
                na.rm = TRUE), medianOther = median(otherClass, na.rm = TRUE), pval = wilcox.test(thisClass,
                otherClass)$p.value)
    }), use.names = TRUE, fill = TRUE)
    res[, `:=`(padj, p.adjust(pval, method = "fdr"))]
    res_all <- copy(res)
    res <- res[padj <= fdr_cutoff]
    res[, `:=`(effect, meanThis - meanOther)]
    data.table::setorder(res, effect)
    if (nrow(res) == 0) {
        warning("No HMDB classes passed enrichment FDR.")
        return(list(plot = NULL, enrichment_results = res, enrichment_results_all = res_all, mapped_results = excr))
    }
    plot_dt <- data.table::rbindlist(lapply(res$class, function(testClass) {
        hmdbs <- unique(hmdb_long[class_label == testClass, hmdb])
        vals <- excr[HMDB %in% hmdbs, .(HMDB, metabr, reaction, beta_meta)]
        vals[, `:=`(class, testClass)]
        vals
    }), use.names = TRUE, fill = TRUE)
    plot_dt[, `:=`(class, factor(class, levels = res$class))]
    plot_dt <- merge(plot_dt, res[, .(class, padj)], by = "class", all.x = TRUE)
    enrichment_plot <- ggplot(plot_dt, aes(x = class, y = beta_meta, fill = padj)) + geom_boxplot(outlier.shape = 16,
        outlier.size = 1.2, outlier.alpha = 0.75, color = "#3A3A3A", linewidth = 0.45) + geom_hline(yintercept = 0,
        linetype = "dashed", color = "#8C8C8C", linewidth = 0.5) + scale_fill_gradient(low = "#5A4E6D",
        high = "#D9D9D9", name = "FDR", guide = guide_colorbar(reverse = TRUE)) + guides(fill = guide_colorbar(title.position = "top")) +
        labs(x = NULL, y = "PPI effect size") + coord_flip(clip = "off") + theme_minimal(base_size = 13) +
        theme(axis.text.y = element_text(size = 10), legend.position = c(0.8, 0.08), legend.direction = "horizontal",
            legend.background = element_rect(fill = ggplot2::alpha("white", 0.75), color = NA), legend.key.width = unit(0.55,
                "cm"), legend.key.height = unit(0.22, "cm"), legend.title = element_text(size = 8), legend.text = element_text(size = 7),
            plot.margin = margin(5.5, 20, 5.5, 5.5))
    list(plot = enrichment_plot, enrichment_results = res, enrichment_results_all = res_all, plot_data = plot_dt,
        mapped_results = excr)
}
metac_blood_enrich <- performBloodFBAEnrichmentPlot(res_tbl = metac_blood_fba_results_annot, map_file = map_file,
    hmdb_file = hmdb_file, min_class_size = enrichment_min_class_size, fdr_cutoff = enrichment_fdr_cutoff,
    effect_col = "beta")
print(metac_blood_enrich$enrichment_results)
if (!is.null(metac_blood_enrich$plot)) {
    print(metac_blood_enrich$plot)
}
write.csv(metac_blood_fba_flux_wide, project_output(file.path(blood_out_dir, "MetaCardis_Blood_FBA_v2_fluxes_mapped_wide.csv"),
    name = "meta_cardis_blood_fba_fluxes_mapped_wide.csv"), row.names = FALSE)
write.csv(metac_blood_pca_var_summary, project_output(file.path(blood_out_dir, "MetaCardis_Blood_FBA_v2_PCA_variance_explained.csv"),
    name = "meta_cardis_blood_fba_pca_variance_explained.csv"), row.names = FALSE)
write.csv(top_metac_blood_pos, project_output(file.path(blood_out_dir, "MetaCardis_Blood_FBA_v2_PC1_top_positive_loadings.csv"),
    name = "meta_cardis_blood_fba_pc1_top_positive_loadings.csv"), row.names = FALSE)
write.csv(top_metac_blood_neg, project_output(file.path(blood_out_dir, "MetaCardis_Blood_FBA_v2_PC1_top_negative_loadings.csv"),
    name = "meta_cardis_blood_fba_pc1_top_negative_loadings.csv"), row.names = FALSE)
write.csv(metac_blood_fba_results, project_output(file.path(blood_out_dir, "MetaCardis_Blood_FBA_v2_HC3_PPI_results.csv"),
    name = "meta_cardis_blood_fba_hc3_ppi_results.csv"), row.names = FALSE)
write.csv(metac_blood_fba_results_annot, project_output(file.path(blood_out_dir, "MetaCardis_Blood_FBA_v2_HC3_PPI_results_annotated.csv"),
    name = "meta_cardis_blood_fba_hc3_ppi_results_annotated.csv"), row.names = FALSE)
write.csv(metac_blood_result_summary, project_output(file.path(blood_out_dir, "MetaCardis_Blood_FBA_v2_HC3_PPI_summary.csv"),
    name = "meta_cardis_blood_fba_hc3_ppi_summary.csv"), row.names = FALSE)
write.csv(cor_summary, project_output(file.path(blood_out_dir, "MetaCardis_Blood_FBA_v2_correlation_summary.csv"),
    name = "meta_cardis_blood_fba_correlation_summary.csv"), row.names = FALSE)
write.csv(positive_high_cor_df, project_output(file.path(blood_out_dir, "MetaCardis_Blood_FBA_v2_high_positive_correlation_pairs.csv"),
    name = "meta_cardis_blood_fba_high_positive_correlation_pairs.csv"), row.names = FALSE)
write.csv(negative_high_cor_df, project_output(file.path(blood_out_dir, "MetaCardis_Blood_FBA_v2_high_negative_correlation_pairs.csv"),
    name = "meta_cardis_blood_fba_high_negative_correlation_pairs.csv"), row.names = FALSE)
write.csv(cluster_df, project_output(file.path(blood_out_dir, "MetaCardis_Blood_FBA_v2_reaction_clusters_cor_0.9999.csv"),
    name = "meta_cardis_blood_fba_reaction_clusters_cor_0_9999.csv"), row.names = FALSE)
write.csv(cluster_summary, project_output(file.path(blood_out_dir, "MetaCardis_Blood_FBA_v2_reaction_cluster_summary.csv"),
    name = "meta_cardis_blood_fba_reaction_cluster_summary.csv"), row.names = FALSE)
write.csv(metac_blood_fba_results_clustered, project_output(file.path(blood_out_dir, "MetaCardis_Blood_FBA_v2_HC3_results_annotated_clustered.csv"),
    name = "meta_cardis_blood_fba_hc3_results_annotated_clustered.csv"), row.names = FALSE)
write.csv(metac_blood_sig_cluster_summary, project_output(file.path(blood_out_dir, "MetaCardis_Blood_FBA_v2_significant_cluster_summary.csv"),
    name = "meta_cardis_blood_fba_significant_cluster_summary.csv"), row.names = FALSE)
write.csv(duplicate_beta_p_summary, project_output(file.path(blood_out_dir, "MetaCardis_Blood_FBA_v2_duplicate_beta_pvalue_summary.csv"),
    name = "meta_cardis_blood_fba_duplicate_beta_pvalue_summary.csv"), row.names = FALSE)
write.csv(metac_blood_enrich$enrichment_results, project_output(file.path(blood_out_dir, "MetaCardis_Blood_FBA_v2_HMDB_enrichment_significant.csv"),
    name = "meta_cardis_blood_fba_hmdb_enrichment_significant.csv"), row.names = FALSE)
write.csv(metac_blood_enrich$enrichment_results_all, project_output(file.path(blood_out_dir, "MetaCardis_Blood_FBA_v2_HMDB_enrichment_all_classes.csv"),
    name = "meta_cardis_blood_fba_hmdb_enrichment_all_classes.csv"), row.names = FALSE)
write.csv(metac_blood_enrich$mapped_results, project_output(file.path(blood_out_dir, "MetaCardis_Blood_FBA_v2_HMDB_mapped_results.csv"),
    name = "meta_cardis_blood_fba_hmdb_mapped_results.csv"), row.names = FALSE)
ggsave(project_output(file.path(blood_out_dir, "MetaCardis_Blood_FBA_v2_PCA_variance_explained.pdf"), name = "meta_cardis_blood_fba_pca_variance_explained.pdf"),
    p_pca_var, width = 7, height = 7)
ggsave(project_output(file.path(blood_out_dir, "MetaCardis_Blood_FBA_v2_PCA_by_PPI.pdf"), name = "meta_cardis_blood_fba_pca_by_ppi.pdf"),
    p_pca_ppi, width = 7, height = 6)
ggsave(project_output(file.path(blood_out_dir, "MetaCardis_Blood_FBA_v2_PCA_by_participant_Gender.pdf"),
    name = "meta_cardis_blood_fba_pca_by_participant_gender.pdf"), p_pca_gender, width = 7, height = 6)
if (!is.null(metac_blood_enrich$plot)) {
    ggsave(project_output(file.path(blood_out_dir, "MetaCardis_Blood_FBA_v2_HMDB_enrichment_plot.pdf"), name = "meta_cardis_blood_fba_hmdb_enrichment_plot.pdf"),
        metac_blood_enrich$plot, width = 8, height = 6)
}
if (length(cluster_plots) > 0) {
    for (cluster_name in names(cluster_plots)) {
        ggsave(project_output(file.path(blood_out_dir, paste0("MetaCardis_Blood_FBA_v2_cluster_boxplot_",
            cluster_name, ".pdf"))), cluster_plots[[cluster_name]], width = 10, height = 7)
    }
}
saveRDS(list(analysis_version = analysis_version, mat_dir = mat_dir, mat_files = mat_files, file_check = file_check,
    mat_read_errors = mat_read_errors, model_sex_check = model_sex_check, metac_blood_fba_flux_wide = metac_blood_fba_flux_wide,
    metac_blood_fba = metac_blood_fba, metac_blood_meta = metac_blood_meta, mat = mat, mat_log = mat_log,
    pca_unscaled = metac_blood_pca_unscaled, pca_scaled = metac_blood_pca_scaled, pca_var_summary = metac_blood_pca_var_summary,
    pca_scores = metac_blood_scores_scaled, top_positive_loadings = top_metac_blood_pos, top_negative_loadings = top_metac_blood_neg,
    full_results = metac_blood_fba_results, annotated_results = metac_blood_fba_results_annot, correlation_summary = cor_summary,
    positive_high_cor_pairs = positive_high_cor_df, negative_high_cor_pairs = negative_high_cor_df, cluster_df = cluster_df,
    cluster_summary = cluster_summary, clustered_results = metac_blood_fba_results_clustered, significant_cluster_summary = metac_blood_sig_cluster_summary),
    project_output(file.path(blood_out_dir, "MetaCardis_Blood_FBA_v2_full_workflow_results.rds"), name = "meta_cardis_blood_fba_full_workflow_results.rds"))
message("MetaCardis Blood FBA v2 workflow complete.")
message("Results saved to: ", blood_out_dir)
}))
