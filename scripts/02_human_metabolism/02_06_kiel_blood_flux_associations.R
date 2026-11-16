# SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0
# Kiel blood-flux primary HC3 regression and PCA diagnostics.
# Run through run.R; inputs and exports use package-local routing.
# Validated from 1,280 original flux files: all 317 primary reactions and 299 downstream Kiel rows match.
# Metadata loading and medication-count preparation are explicit. Participant data/models are not distributed.
# This primary script does not run sex-stratified or PPI-by-sex interaction models.

ppi_source(ppi_code("scripts/00_shared_helpers/00_01_load_kiel_metadata.R"))
focus_data <- cbind(data, phenotype_all_drugs_filt[, 4:16], phenoData_filt[, c(2, 8, 9, 18, 68)])
focus_data <- as.data.frame(focus_data)
focus_data$glucose_mg.dl_mb <- suppressWarnings(as.numeric(focus_data$glucose_mg.dl_mb))
focus_data$diabetes_mb <- as.factor(focus_data$diabetes_mb)
focus_data$BMI.g25.30 <- as.factor(focus_data$BMI.g25.30)
focus_data$Analgesic <- as.factor(focus_data$Analgesic)
focus_data$Analgesic_NSAID <- as.factor(focus_data$Analgesic_NSAID)
focus_data$Antibiotics <- as.factor(focus_data$Antibiotics)
focus_data$Antihypertensives <- as.factor(focus_data$Antihypertensives)
focus_data$Antidiabetic <- as.factor(focus_data$Antidiabetic)
focus_data$Antidepressant <- as.factor(focus_data$Antidepressant)
focus_data$Antiphlogistic <- as.factor(focus_data$Antiphlogistic)
focus_data$Metformin <- as.factor(focus_data$Metformin)
med_broad <- c("Analgesic", "Antibiotics", "Antidepressant", "Antihypertensives", "Antiphlogistic", "Antidiabetic")
stopifnot(nrow(phenotype_all_drugs_filt) == nrow(focus_data))
focus_data[med_broad] <- phenotype_all_drugs_filt[, med_broad]
focus_data[med_broad] <- lapply(focus_data[med_broad], function(x) as.numeric(as.character(x)))
focus_data$medication_count <- rowSums(focus_data[, med_broad, drop = FALSE], na.rm = TRUE)
focus_data$Analgesic <- as.factor(focus_data$Analgesic)
focus_data$Antibiotics <- as.factor(focus_data$Antibiotics)
focus_data$Antidepressant <- as.factor(focus_data$Antidepressant)
focus_data$Antihypertensives <- as.factor(focus_data$Antihypertensives)
focus_data$Antiphlogistic <- as.factor(focus_data$Antiphlogistic)
focus_data$Antidiabetic <- as.factor(focus_data$Antidiabetic)
focus_data$diabetes_mb <- as.factor(focus_data$diabetes_mb)
ppi_term <- "PPI1"
library(R.matlab)
library(dplyr)
library(purrr)
library(tibble)
library(tidyr)
mat_dir <- ppi_path("human_cohorts/blood_flux/kiel_models/resultFlux")
mat_files <- ppi_list_files(mat_dir, pattern = "\\.mat$", full.names = TRUE)
length(mat_files)
extract_rxns <- function(rxns_obj) {
    vapply(seq_along(rxns_obj), function(i) {
        as.character(rxns_obj[[i]][[1]][1, 1])
    }, character(1))
}
extract_one_string <- function(x) {
    as.character(x[[1]][[1]][1, 1])
}
extract_one_fba_mat <- function(file) {
    mat <- R.matlab::readMat(ppi_input(file))
    rxns <- extract_rxns(mat$rxns)
    flux <- as.numeric(mat$f)
    sample_id <- extract_one_string(mat$ID)
    sex <- extract_one_string(mat$sex)
    if (length(rxns) != length(flux)) {
        stop("Length mismatch in file: ", basename(file), "\nrxns = ", length(rxns), "\nf = ", length(flux))
    }
    if (any(!nzchar(rxns))) {
        stop("Empty reaction names in file: ", basename(file))
    }
    out <- as.data.frame(matrix(flux, nrow = 1), check.names = FALSE)
    names(out) <- rxns
    out <- out %>% mutate(sampleID = sample_id, sex = sex, file = basename(file), .before = 1)
    out
}
{
    previous_mat_locale <- Sys.getlocale("LC_CTYPE")
    invisible(Sys.setlocale("LC_CTYPE", "C"))
    tryCatch(blood_fba_flux_wide <- purrr::map_dfr(mat_files, extract_one_fba_mat), finally = invisible(Sys.setlocale("LC_CTYPE", 
        previous_mat_locale)))
}
dim(blood_fba_flux_wide)
bram_dict <- readxl::read_xlsx(ppi_input(ppi_path("Bram/SampleIDTranslateTable.xlsx")))
anon_dict <- read.csv(ppi_input(ppi_path("data/Community_sim/FOCUS/FOCUS/focus_dictionary.csv")))
bram_dict$orig <- anon_dict$x
bram_dict$orig_anon <- anon_dict$anon_out
library(stringr)
library(dplyr)
blood_fba_flux_wide <- blood_fba_flux_wide %>% mutate(bram_sample = str_extract(file, "samp[0-9]+"), 
    original_sampleID = bram_dict$orig[match(bram_sample, bram_dict$`New sample ID`)], original_anon = bram_dict$orig_anon[match(bram_sample, 
        bram_dict$`New sample ID`)]) %>% relocate(original_sampleID, original_anon, bram_sample, .before = 1)
blood_fba_flux_wide <- blood_fba_flux_wide %>% rename(mat_sampleID = sampleID) %>% rename(sampleID = original_sampleID)
library(dplyr)
library(tidyr)
library(matrixStats)
library(irlba)
library(ggplot2)
library(lmtest)
library(sandwich)
library(tibble)
rxn_cols <- setdiff(names(blood_fba_flux_wide), c("sampleID", "mat_sampleID", "original_anon", "bram_sample", 
    "sex", "file"))
blood_fba_flux_wide$n_missing <- rowSums(is.na(blood_fba_flux_wide[, rxn_cols, drop = FALSE]))
blood_fba_flux_wide$all_nan <- blood_fba_flux_wide$n_missing == length(rxn_cols)
blood_fba <- blood_fba_flux_wide %>% filter(!all_nan)
blood_meta <- focus_data[match(blood_fba$sampleID, focus_data$ID), , drop = FALSE]
stopifnot(all(blood_meta$ID == blood_fba$sampleID))
mat <- as.matrix(blood_fba[, rxn_cols, drop = FALSE])
storage.mode(mat) <- "numeric"
col_missing <- colMeans(is.na(mat))
keep_complete <- col_missing == 0
mat <- mat[, keep_complete, drop = FALSE]
vars <- matrixStats::colVars(mat, na.rm = TRUE)
keep_var <- vars > 0
mat <- mat[, keep_var, drop = FALSE]
flux_summary <- data.frame(statistic = c("min", "q25", "median", "mean", "q75", "max"), value = as.numeric(summary(as.vector(mat))))
quantile(abs(as.vector(mat)), probs = c(0.5, 0.9, 0.99, 0.999), na.rm = TRUE)
mat_log <- sign(mat) * log1p(abs(mat))
quantile(abs(as.vector(mat_log)), probs = c(0.5, 0.9, 0.99, 0.999), na.rm = TRUE)
blood_pca_unscaled <- irlba::prcomp_irlba(mat_log, n = 10, center = TRUE, scale. = FALSE)
blood_pca_scaled <- irlba::prcomp_irlba(mat_log, n = 10, center = TRUE, scale. = TRUE)
var_unscaled <- blood_pca_unscaled$sdev^2/sum(blood_pca_unscaled$sdev^2)
var_scaled <- blood_pca_scaled$sdev^2/sum(blood_pca_scaled$sdev^2)
blood_pca_var_summary <- bind_rows(data.frame(model = "Blood_FBA_unscaled", PC = paste0("PC", seq_along(var_unscaled)), 
    variance_explained = var_unscaled, percent_variance_explained = 100 * var_unscaled, cumulative_percent = 100 * 
        cumsum(var_unscaled)), data.frame(model = "Blood_FBA_scaled", PC = paste0("PC", seq_along(var_scaled)), 
    variance_explained = var_scaled, percent_variance_explained = 100 * var_scaled, cumulative_percent = 100 * 
        cumsum(var_scaled)))
blood_scores_scaled <- as.data.frame(blood_pca_scaled$x)
blood_scores_scaled$sampleID <- blood_fba$sampleID
blood_scores_scaled <- blood_scores_scaled %>% left_join(blood_meta %>% dplyr::select(ID, PPI, Age, BMI, 
    Gender, diabetes_mb, medication_count), by = c(sampleID = "ID"))
ggplot(blood_pca_var_summary, aes(x = as.numeric(gsub("PC", "", PC)), y = percent_variance_explained)) + 
    geom_col() + facet_wrap(~model, ncol = 1) + theme_bw() + labs(title = "Blood FBA PCA variance explained", 
    x = "Principal component", y = "% variance explained")
ggplot(blood_scores_scaled, aes(x = PC1, y = PC2, color = PPI)) + geom_point(alpha = 0.75) + theme_bw() + 
    labs(title = "Blood FBA scaled PCA scores", x = "PC1", y = "PC2", color = "PPI")
ggplot(blood_scores_scaled, aes(x = PC1, y = PC2, color = Gender)) + geom_point(alpha = 0.75) + theme_bw() + 
    labs(title = "Blood FBA scaled PCA scores by sex", x = "PC1", y = "PC2", color = "Gender")
blood_meta_pc1 <- blood_meta
blood_meta_pc1$PC1_scaled <- blood_pca_scaled$x[, 1]
blood_meta_pc1$PC2_scaled <- blood_pca_scaled$x[, 2]
blood_loadings <- blood_pca_scaled$rotation[, 1]
names(blood_loadings) <- colnames(mat_log)
mean(blood_loadings > 0)
top_blood_pos <- data.frame(reaction = names(sort(blood_loadings, decreasing = TRUE))[1:30], loading = sort(blood_loadings, 
    decreasing = TRUE)[1:30])
top_blood_neg <- data.frame(reaction = names(sort(blood_loadings))[1:30], loading = sort(blood_loadings)[1:30])
cor_mat <- cor(mat_log, use = "pairwise.complete.obs")
cor_values <- cor_mat[upper.tri(cor_mat)]
fml_blood <- met ~ PPI + Age + BMI + Gender + diabetes_mb + medication_count + Metformin + STATIN + ASPIRIN
ppi_term <- "PPI1"
run_blood_hc3 <- function(mat_log, metadata, formula, ppi_term = "PPI1") {
    covars <- setdiff(all.vars(formula), "met")
    missing_covars <- setdiff(covars, names(metadata))
    if (length(missing_covars) > 0) {
        stop("Missing covariates: ", paste(missing_covars, collapse = ", "))
    }
    keep <- complete.cases(metadata[, covars, drop = FALSE])
    mat_use <- mat_log[keep, , drop = FALSE]
    meta_use <- metadata[keep, , drop = FALSE]
    fit_one <- function(y) {
        df <- meta_use
        df$met <- y
        if (sd(df$met, na.rm = TRUE) == 0) {
            return(c(beta = NA, se = NA, t = NA, p = NA, n = nrow(df)))
        }
        fit <- tryCatch(lm(formula, data = df), error = function(e) NULL)
        if (is.null(fit)) {
            return(c(beta = NA, se = NA, t = NA, p = NA, n = nrow(df)))
        }
        ct <- tryCatch(lmtest::coeftest(fit, vcov. = sandwich::vcovHC(fit, type = "HC3")), error = function(e) NULL)
        if (is.null(ct) || !(ppi_term %in% rownames(ct))) {
            return(c(beta = NA, se = NA, t = NA, p = NA, n = nrow(df)))
        }
        c(beta = ct[ppi_term, "Estimate"], se = ct[ppi_term, "Std. Error"], t = ct[ppi_term, "t value"], 
            p = ct[ppi_term, "Pr(>|t|)"], n = nrow(df))
    }
    res <- lapply(seq_len(ncol(mat_use)), function(i) {
        fit_one(mat_use[, i])
    })
    res_df <- as.data.frame(do.call(rbind, res))
    tibble(reaction = colnames(mat_use), beta = as.numeric(res_df$beta), SE = as.numeric(res_df$se), 
        t_value = as.numeric(res_df$t), pvalue = as.numeric(res_df$p), n = as.numeric(res_df$n)) %>% 
        mutate(FDR = p.adjust(pvalue, method = "fdr")) %>% arrange(FDR, pvalue)
}
blood_fba_results <- run_blood_hc3(mat_log = mat_log, metadata = blood_meta, formula = fml_blood, ppi_term = "PPI1")
sum(blood_fba_results$pvalue < 0.05, na.rm = TRUE)
sum(blood_fba_results$FDR < 0.05, na.rm = TRUE)
blood_out_dir <- ppi_path("human_cohorts/blood_flux/kiel")
dir.create(ppi_output(blood_out_dir, directory = TRUE), showWarnings = FALSE, recursive = TRUE)
write.csv(blood_pca_var_summary, ppi_output(file.path(blood_out_dir, "Blood_FBA_PCA_variance_explained.csv")), 
    row.names = FALSE)
write.csv(blood_fba_results, ppi_output(file.path(blood_out_dir, "Blood_FBA_HC3_PPI_results.csv")), row.names = FALSE)
saveRDS(list(blood_fba = blood_fba, blood_meta = blood_meta, mat_log = mat_log, pca_scaled = blood_pca_scaled, 
    pca_unscaled = blood_pca_unscaled, pca_var_summary = blood_pca_var_summary, blood_fba_results = blood_fba_results, 
    top_positive_loadings = top_blood_pos, top_negative_loadings = top_blood_neg), ppi_output(file.path(blood_out_dir, 
    "Blood_FBA_diagnostics_and_HC3_results.rds")))
