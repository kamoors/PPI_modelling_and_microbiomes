# SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0
# Kiel final secretion regression
# Run through run.R so input and output paths resolve within the package.
# Inputs: see data/input_manifest.csv and docs/RUN_ORDER.md.
# Human participant-level input matrices and covariates are not distributed.
# Supply authorised data separately; see docs/RESTRICTED_DATA.md.

library(lmtest)
library(sandwich)
require(gplots)
require(ggfortify)
require(ggplot2)
require(readxl)
require(writexl)
library(reshape2)
library(ggsignif)
library(dplyr)
library(stringr)
library(RColorBrewer)
library(car)
library(boot.pval)
library(MatrixGenerics)
library(tibble)
ppi_source(ppi_code("scripts/00_shared_helpers/00_01_load_kiel_metadata.R"))
ppi_setwd(ppi_path("human_cohorts/kiel/"))
mg_sec <- read.csv(ppi_input(ppi_path("human_cohorts/kiel/Secretion_data_mgpipe_focus.csv")))
colnames(mg_sec) <- c("Net_secretion", stringr::str_sub(colnames(mg_sec[, -1]), 2))
stopifnot("Net_secretion" %in% colnames(mg_sec))
missing_in_mg <- setdiff(data$ID, colnames(mg_sec)[-1])
extra_in_mg <- setdiff(colnames(mg_sec)[-1], data$ID)
if (length(missing_in_mg) > 0) stop("Missing samples in mg_sec: ", paste(missing_in_mg, collapse = ", "))
if (length(extra_in_mg) > 0) message("Extra samples in mg_sec (will be dropped): ", paste(extra_in_mg, 
    collapse = ", "))
mg_sec <- mg_sec[, c("Net_secretion", data$ID)]
stopifnot(all(colnames(mg_sec[, -1]) == data$ID))
boot_dat <- t(mg_sec[, -1])
colnames(boot_dat) <- mg_sec[, 1]
boot_dat[is.na(boot_dat)] <- 0
boot_dat_filt <- boot_dat[, colVars(boot_dat) > 0, drop = FALSE]
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
fit_all_features_hc3 <- function(data_in, test_df, fml, transform_y = function(x) x) {
    out <- vector("list", ncol(data_in))
    for (i in seq_len(ncol(data_in))) {
        test_df$met <- transform_y(data_in[, i])
        fit <- lm(fml, data = test_df)
        out[[i]] <- lmtest::coeftest(fit, vcov. = sandwich::vcovHC(fit, type = "HC3"))
    }
    names(out) <- colnames(data_in)
    out
}
ppi_source(ppi_code("scripts/00_shared_helpers/00_07_secretion_bootstrap_functions.R"))
fml_full <- met ~ PPI + Age + Gender + BMI + diabetes_mb + medication_count + Metformin + STATIN + ASPIRIN
fml_base <- met ~ PPI + Age + Gender + BMI
fml_boot_full <- fml_full
fml_boot_base <- fml_base
y_transform <- function(x) log1p(x)
test_df <- focus_data
data_in <- boot_dat_filt
mg_sec_glm_full <- fit_all_features_hc3(data_in, test_df, fml_full, transform_y = y_transform)
focus_mgsec_results_preboot_full <- tibble(metabolite = names(mg_sec_glm_full), hc3_coef = sapply(mg_sec_glm_full, 
    function(ct) {
        if (is.null(ct)) 
            return(NA_real_)
        if (ppi_term %in% rownames(ct)) 
            ct[ppi_term, "Estimate"]
        else NA_real_
    }), hc3_se = sapply(mg_sec_glm_full, function(ct) {
    if (is.null(ct)) 
        return(NA_real_)
    if (ppi_term %in% rownames(ct)) 
        ct[ppi_term, "Std. Error"]
    else NA_real_
}), hc3_p = sapply(mg_sec_glm_full, function(ct) {
    if (is.null(ct)) 
        return(NA_real_)
    if (ppi_term %in% rownames(ct)) 
        ct[ppi_term, "Pr(>|t|)"]
    else NA_real_
})) %>% mutate(hc3_ci_low = hc3_coef - qnorm(0.975) * hc3_se, hc3_ci_high = hc3_coef + qnorm(0.975) * 
    hc3_se, hc3_fdr = p.adjust(hc3_p, method = "fdr")) %>% arrange(hc3_fdr, hc3_p)
load_metabolite_dict <- function(dict_path) {
    md <- readr::read_tsv(ppi_input(dict_path), show_col_types = FALSE)
    stopifnot(all(c("abbreviation", "fullName") %in% colnames(md)))
    md %>% dplyr::select(abbreviation, fullName) %>% dplyr::mutate(abbreviation = as.character(abbreviation), 
        fullName = as.character(fullName)) %>% dplyr::distinct(abbreviation, .keep_all = TRUE)
}
extract_met_abbr <- function(x) {
    x <- as.character(x)
    out <- str_match(x, "^EX_([^\\[]+)\\[")[, 2]
    idx_na <- is.na(out) | out == ""
    if (any(idx_na)) {
        out[idx_na] <- str_match(x[idx_na], "^[^_]+_([^\\[]+)\\[")[, 2]
    }
    idx_na <- is.na(out) | out == ""
    if (any(idx_na)) {
        tmp <- x[idx_na]
        tmp <- str_replace(tmp, "\\[.*\\]$", "")
        tmp <- str_replace(tmp, "^(EX|DM|SK|SINK|sink)_", "")
        tmp <- str_replace(tmp, "^[^_]+_", "")
        out[idx_na] <- tmp
    }
    out
}
dict_path <- ppi_path("human_cohorts/kiel/recon-store-metabolites-1 (1).tsv")
md <- load_metabolite_dict(dict_path)
rxn_col = "metabolite"
coef_col = "hc3_coef"
focus_mgsec_results_preboot_full <- focus_mgsec_results_preboot_full %>% dplyr::mutate(rxn = .data[[rxn_col]], 
    coef = as.numeric(.data[[coef_col]]), met_abbr = extract_met_abbr(.data[[rxn_col]])) %>% dplyr::left_join(md, 
    by = c(met_abbr = "abbreviation")) %>% dplyr::mutate(label = dplyr::coalesce(fullName, met_abbr, 
    rxn), label = as.character(label))
write.csv(focus_mgsec_results_preboot_full, file = ppi_output(ppi_path("human_cohorts/results/Focus_mgpipe_net_secretion_HC3_FULL_NEWPPI_clean_results.csv"), 
    name = "kiel_metabolite_secretion_hc3_full.csv"), row.names = FALSE)
ppi_source(ppi_code("scripts/00_shared_helpers/00_08_secretion_additional_analyses.R"))
