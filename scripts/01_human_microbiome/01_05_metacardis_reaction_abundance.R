# SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0
# MetaCardis reaction abundance: manuscript-matching primary HC3 model.
# The full reaction-abundance model does not include systolic blood pressure.
# Retains exclusion of reported non-users with metabolomics-detected PPI; natural log(x + 1e-6), not CLR.
# Participant matrices and covariates are unavailable publicly. No bootstrap or historical BASE analyses are run.

source("scripts/helpers/setup.R")
invisible(run_analysis(
  section = "01_human_microbiome",
  directory = "human_cohorts/metacardis",
  script = "01_05_metacardis_reaction_abundance.R",
  restricted = TRUE,
  code = {
library(lmtest)
library(sandwich)
library(readxl)
library(dplyr)
library(stringr)
library(MatrixGenerics)
library(tibble)
library(matrixStats)
meta_dir <- project_path("Metacardis/data/MetaCardis")
results_dir <- project_path("human_cohorts/results")
project_setwd(meta_dir)
cardis_meta <- readxl::read_excel(project_input("Merged_Metacardis-Study.xlsx"))
cardis_meta$ID <- trimws(cardis_meta$ID)
met_data_raw <- cardis_meta[, 791:(ncol(cardis_meta) - 3)]
met_data_num <- met_data_raw %>% mutate(across(everything(), ~suppressWarnings(as.numeric(.))))
met_matrix_all <- as.matrix(met_data_num)
met_all_na_rows <- rowSums(is.na(met_matrix_all)) == ncol(met_matrix_all)
cardis_meta_met <- cardis_meta[!met_all_na_rows, , drop = FALSE]
met_matrix_all <- met_matrix_all[!met_all_na_rows, , drop = FALSE]
stopifnot(nrow(cardis_meta_met) == nrow(met_matrix_all))
ppi_cols <- c("omeprazole", "pantoprazole", "dexlansoprazole")
ppi_cols <- ppi_cols[ppi_cols %in% colnames(met_matrix_all)]
if (length(ppi_cols) == 0) {
    stop("None of the expected PPI metabolite columns were found in metabolomics data.")
}
ppi_signal <- met_matrix_all[, ppi_cols, drop = FALSE]
ppi_baseline <- matrixStats::colMins(ppi_signal, na.rm = TRUE)
ppi_detect_matrix <- sweep(ppi_signal, 2, ppi_baseline, FUN = function(x, base) x > base + 0.1)
ppi_detected <- rowSums(ppi_detect_matrix, na.rm = TRUE) > 0
ppi_measured <- ifelse(ppi_detected, "Yes", "No")
met_covariates_all <- cardis_meta_met %>% dplyr::select(ID, Proton.pump.inhibitor) %>% as.data.frame()
met_covariates_all$PPI_measured <- factor(ppi_measured, levels = c("No", "Yes"))
met_covariates_all$Proton.pump.inhibitor <- factor(met_covariates_all$Proton.pump.inhibitor, levels = c("No",
    "Yes"))
keep_met_samples <- !(met_covariates_all$Proton.pump.inhibitor == "No" & met_covariates_all$PPI_measured ==
    "Yes")
met_covariates_clean <- met_covariates_all[keep_met_samples, , drop = FALSE]
clean_ids <- met_covariates_clean$ID
reaction_abundance <- read.csv(project_input(file.path(meta_dir, "ReactionAbundance.csv")), check.names = FALSE)
colnames(reaction_abundance)[1] <- "rxn_name"
reaction_abundance$rxn_name <- trimws(as.character(reaction_abundance$rxn_name))
rxn_sample_ids <- colnames(reaction_abundance)[-1]
rxn_sample_ids <- sub("^X", "", rxn_sample_ids)
rxn_sample_ids <- trimws(rxn_sample_ids)
colnames(reaction_abundance) <- c("rxn_name", rxn_sample_ids)
common_ids <- clean_ids[clean_ids %in% colnames(reaction_abundance)]
cardis_meta_clean <- cardis_meta %>% dplyr::filter(ID %in% common_ids) %>% dplyr::arrange(match(ID, common_ids))
reaction_abundance_clean <- reaction_abundance[, c("rxn_name", common_ids), drop = FALSE]
stopifnot(all(colnames(reaction_abundance_clean)[-1] == cardis_meta_clean$ID))
med_vars <- c("Antidiabetic.drugs", "Insulin", "Metformin", "Sulfonurea", "DPP.IV.inhibitors", "GLP1.agonists",
    "Anti-hypertensive.drugs", "Beta-blockers", "ACE", "AT2.inhibitor", "Calcium.blockers", "Anse.diuretics",
    "Potassium-.sparing.diuretic", "Thiazides", "Lipid.lowering.drugs", "Statin", "Ezetimib", "Fibrate",
    "Nitrate", "Amiodarone", "Digitalic", "Aspirine", "Clopidogrel", "Coumarin")
cardis_meta_clean[med_vars] <- lapply(cardis_meta_clean[med_vars], function(x) {
    x <- as.character(x)
    ifelse(x %in% c("Yes", "1"), 1, ifelse(x %in% c("No", "0"), 0, NA))
})
med_broad <- c("Antidiabetic.drugs", "Anti-hypertensive.drugs", "Lipid.lowering.drugs", "Nitrate", "Amiodarone",
    "Digitalic", "Aspirine", "Clopidogrel", "Coumarin")
cardis_meta_clean$medication_count <- rowSums(cardis_meta_clean[med_broad], na.rm = TRUE)
rxn_abundance_dat <- t(reaction_abundance_clean[, -1, drop = FALSE])
colnames(rxn_abundance_dat) <- reaction_abundance_clean$rxn_name
rxn_abundance_dat <- as.matrix(rxn_abundance_dat)
storage.mode(rxn_abundance_dat) <- "numeric"
rxn_abundance_dat[is.na(rxn_abundance_dat)] <- 0
rxn_abundance_dat_filt <- rxn_abundance_dat[, MatrixGenerics::colVars(rxn_abundance_dat) > 0, drop = FALSE]
clr_pseudocount <- 1e-06
rxn_abundance_dat_clr <- log(rxn_abundance_dat_filt + clr_pseudocount)
rxn_abundance_dat_clr <- as.matrix(rxn_abundance_dat_clr)
cardis_meta_clean <- as.data.frame(cardis_meta_clean)
cardis_meta_clean$Proton.pump.inhibitor <- factor(cardis_meta_clean$Proton.pump.inhibitor, levels = c("No",
    "Yes"))
cardis_meta_clean$`Age.(years)` <- suppressWarnings(as.numeric(cardis_meta_clean$`Age.(years)`))
cardis_meta_clean$`BMI.(kg/m²)` <- suppressWarnings(as.numeric(cardis_meta_clean$`BMI.(kg/m²)`))
cardis_meta_clean$`Systolic.blood.pressure.(mmHg)` <- suppressWarnings(as.numeric(cardis_meta_clean$`Systolic.blood.pressure.(mmHg)`))
cardis_meta_clean$creatinine <- suppressWarnings(as.numeric(cardis_meta_clean$creatinine))
cardis_meta_clean$medication_count <- suppressWarnings(as.numeric(cardis_meta_clean$medication_count))
cardis_meta_clean$Gender <- as.factor(cardis_meta_clean$Gender)
cardis_meta_clean$Nationality <- as.factor(cardis_meta_clean$Nationality)
cardis_meta_clean$Diabetic.status <- as.factor(cardis_meta_clean$Diabetic.status)
cardis_meta_clean$Metformin <- as.factor(cardis_meta_clean$Metformin)
cardis_meta_clean$Statin <- as.factor(cardis_meta_clean$Statin)
cardis_meta_clean$Aspirine <- as.factor(cardis_meta_clean$Aspirine)
cardis_meta_clean$Proton.pump.inhibitor <- as.factor(cardis_meta_clean$Proton.pump.inhibitor)
rxn_total_abundance <- rowSums(rxn_abundance_dat_filt, na.rm = TRUE)
cardis_meta_clean$rxn_total_abundance <- rxn_total_abundance
cardis_meta_clean$log_rxn_total_abundance <- log1p(rxn_total_abundance)
ppi_term <- "Proton.pump.inhibitorYes"
fit_all_features_hc3 <- function(data_in, test_df, fml, transform_y = function(x) x) {
    out <- vector("list", ncol(data_in))
    for (i in seq_len(ncol(data_in))) {
        test_df$met <- transform_y(data_in[, i])
        fit <- tryCatch(lm(fml, data = test_df), error = function(e) NULL)
        if (is.null(fit)) {
            out[[i]] <- NULL
        }
        else {
            ct <- tryCatch(lmtest::coeftest(fit, vcov. = sandwich::vcovHC(fit, type = "HC3")), error = function(e) NULL)
            out[[i]] <- list(fit = fit, ct = ct)
        }
    }
    names(out) <- colnames(data_in)
    out
}
fml_full <- met ~ Proton.pump.inhibitor + Gender + `Age.(years)` + `BMI.(kg/m²)` + Nationality + Diabetic.status +
    medication_count + Metformin + Statin + Aspirine + creatinine + log_rxn_total_abundance
y_transform <- function(x) x
vars_full <- all.vars(fml_full)
test_df_full <- cardis_meta_clean %>% dplyr::select(all_of(vars_full[-1])) %>% as.data.frame()
keep_full <- complete.cases(test_df_full)
data_in_full <- rxn_abundance_dat_clr[keep_full, , drop = FALSE]
test_df_full <- cardis_meta_clean[keep_full, , drop = FALSE]
ppi_term <- "Proton.pump.inhibitorYes"
react_ab_glm_full <- fit_all_features_hc3(data_in = data_in_full, test_df = test_df_full, fml = fml_full,
    transform_y = y_transform)
metacardis_react_ab_results_preboot_full <- tibble(rxn = names(react_ab_glm_full), hc3_coef = sapply(react_ab_glm_full,
    function(x) {
        if (is.null(x) || is.null(x$ct))
            return(NA_real_)
        ct <- x$ct
        if (ppi_term %in% rownames(ct))
            ct[ppi_term, "Estimate"]
        else NA_real_
    }), hc3_se = sapply(react_ab_glm_full, function(x) {
    if (is.null(x) || is.null(x$ct))
        return(NA_real_)
    ct <- x$ct
    if (ppi_term %in% rownames(ct))
        ct[ppi_term, "Std. Error"]
    else NA_real_
}), hc3_p = sapply(react_ab_glm_full, function(x) {
    if (is.null(x) || is.null(x$ct))
        return(NA_real_)
    ct <- x$ct
    if (ppi_term %in% rownames(ct))
        ct[ppi_term, "Pr(>|t|)"]
    else NA_real_
}), hc3_df = sapply(react_ab_glm_full, function(x) {
    if (is.null(x) || is.null(x$fit))
        return(NA_real_)
    x$fit$df.residual
})) %>% mutate(hc3_crit = qt(0.975, df = hc3_df), hc3_ci_low = hc3_coef - hc3_crit * hc3_se, hc3_ci_high = hc3_coef +
    hc3_crit * hc3_se, hc3_fdr = p.adjust(hc3_p, method = "fdr")) %>% arrange(hc3_fdr, hc3_p)
write.csv(metacardis_react_ab_results_preboot_full, project_output(file.path(results_dir, "MetaCardis_mgpipe_reaction_abundance_cleanPPI_LOG_HC3_FULL_clean_results.csv"),
    name = "metacardis_reaction_abundance_hc3.csv"), row.names = FALSE)
}))
