# SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0
# MetaCardis predicted secretion: manuscript-matching primary HC3 model.
# Run through run.R; inputs and exports use package-local routing.
# Recovered secretion model with metabolomics-based exclusion restored; systolic BP retained.
# All 269 metabolites match the filtered MetaCardis/Results table and its 253 shared downstream metabolites.
# Private participant inputs are required; no bootstrap or historical BASE analyses are run.

library(lmtest)
library(sandwich)
library(readxl)
library(dplyr)
library(stringr)
library(MatrixGenerics)
library(tibble)
ppi_setwd(ppi_path("Metacardis/data/MetaCardis"))
cardis_meta <- readxl::read_excel(ppi_input("Merged_Metacardis-Study.xlsx"))
cardis_mg <- read.csv(ppi_input("inputDiet_net_secretion_fluxes.csv"), check.names = FALSE)
rownames(cardis_mg) <- cardis_mg$`Net secretion`
mg_sec_cardis <- cardis_mg %>% dplyr::select(-`Net secretion`)
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
colnames(mg_sec_cardis) <- trimws(colnames(mg_sec_cardis))
common_ids <- clean_ids[clean_ids %in% colnames(mg_sec_cardis)]
cardis_meta2 <- cardis_meta %>% dplyr::filter(ID %in% common_ids)
mg_sec_cardis2 <- mg_sec_cardis[, common_ids, drop = FALSE]
stopifnot(all(colnames(mg_sec_cardis2) == cardis_meta2$ID))
missing_in_mg <- setdiff(cardis_meta$ID, colnames(mg_sec_cardis))
extra_in_mg <- setdiff(colnames(mg_sec_cardis), cardis_meta$ID)
boot_dat <- t(mg_sec_cardis2)
colnames(boot_dat) <- rownames(cardis_mg)
boot_dat <- as.matrix(boot_dat)
colnames(boot_dat) <- rownames(cardis_mg)
boot_dat[is.na(boot_dat)] <- 0
boot_dat_filt <- boot_dat[, MatrixGenerics::colVars(boot_dat) > 0, drop = FALSE]
focus_data <- as.data.frame(cardis_meta2)
med_vars <- c("Antidiabetic.drugs", "Insulin", "Metformin", "Sulfonurea", "DPP.IV.inhibitors", "GLP1.agonists", 
    "Anti-hypertensive.drugs", "Beta-blockers", "ACE", "AT2.inhibitor", "Calcium.blockers", "Anse.diuretics", 
    "Potassium-.sparing.diuretic", "Thiazides", "Lipid.lowering.drugs", "Statin", "Ezetimib", "Fibrate", 
    "Nitrate", "Amiodarone", "Digitalic", "Aspirine", "Clopidogrel", "Coumarin")
focus_data[med_vars] <- lapply(focus_data[med_vars], function(x) {
    x <- as.character(x)
    ifelse(x %in% c("Yes", "1"), 1, ifelse(x %in% c("No", "0"), 0, NA))
})
med_broad <- c("Antidiabetic.drugs", "Anti-hypertensive.drugs", "Lipid.lowering.drugs", "Nitrate", "Amiodarone", 
    "Digitalic", "Aspirine", "Clopidogrel", "Coumarin")
focus_data$medication_count <- rowSums(focus_data[med_broad], na.rm = TRUE)
focus_data$Proton.pump.inhibitor <- as.factor(focus_data$Proton.pump.inhibitor)
if (all(c("No", "Yes") %in% unique(na.omit(focus_data$Proton.pump.inhibitor)))) {
    focus_data$Proton.pump.inhibitor <- factor(focus_data$Proton.pump.inhibitor, levels = c("No", "Yes"))
}
focus_data$`Age.(years)` <- suppressWarnings(as.numeric(focus_data$`Age.(years)`))
focus_data$`BMI.(kg/m²)` <- suppressWarnings(as.numeric(focus_data$`BMI.(kg/m²)`))
focus_data$`Systolic.blood.pressure.(mmHg)` <- suppressWarnings(as.numeric(focus_data$`Systolic.blood.pressure.(mmHg)`))
focus_data$creatinine <- suppressWarnings(as.numeric(focus_data$creatinine))
focus_data$medication_count <- suppressWarnings(as.numeric(focus_data$medication_count))
focus_data$Gender <- as.factor(focus_data$Gender)
focus_data$Nationality <- as.factor(focus_data$Nationality)
focus_data$Diabetic.status <- as.factor(focus_data$Diabetic.status)
focus_data$Metformin <- as.factor(focus_data$Metformin)
focus_data$Statin <- as.factor(focus_data$Statin)
focus_data$Aspirine <- as.factor(focus_data$Aspirine)
if ("Anti-hypertensive.drugs" %in% colnames(focus_data)) {
    focus_data$`Anti-hypertensive.drugs` <- as.factor(focus_data$`Anti-hypertensive.drugs`)
}
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
fml_full <- met ~ Proton.pump.inhibitor + Gender + `Age.(years)` + `BMI.(kg/m²)` + Nationality + `Systolic.blood.pressure.(mmHg)` + 
    Diabetic.status + medication_count + Metformin + Statin + Aspirine + creatinine
y_transform <- function(x) log1p(x)
vars_full <- all.vars(fml_full)
test_df_full <- focus_data %>% dplyr::select(all_of(vars_full[-1])) %>% as.data.frame()
keep_full <- complete.cases(test_df_full)
data_in_full <- boot_dat_filt[keep_full, , drop = FALSE]
test_df_full <- focus_data[keep_full, , drop = FALSE]
mg_sec_glm_full <- fit_all_features_hc3(data_in = data_in_full, test_df = test_df_full, fml = fml_full, 
    transform_y = y_transform)
metacardis_secretion_results <- tibble(rxn = names(mg_sec_glm_full), hc3_coef = sapply(mg_sec_glm_full, 
    function(x) {
        if (is.null(x) || is.null(x$ct)) 
            return(NA_real_)
        ct <- x$ct
        if (ppi_term %in% rownames(ct)) 
            ct[ppi_term, "Estimate"]
        else NA_real_
    }), hc3_se = sapply(mg_sec_glm_full, function(x) {
    if (is.null(x) || is.null(x$ct)) 
        return(NA_real_)
    ct <- x$ct
    if (ppi_term %in% rownames(ct)) 
        ct[ppi_term, "Std. Error"]
    else NA_real_
}), hc3_p = sapply(mg_sec_glm_full, function(x) {
    if (is.null(x) || is.null(x$ct)) 
        return(NA_real_)
    ct <- x$ct
    if (ppi_term %in% rownames(ct)) 
        ct[ppi_term, "Pr(>|t|)"]
    else NA_real_
}), hc3_df = sapply(mg_sec_glm_full, function(x) {
    if (is.null(x) || is.null(x$fit)) 
        return(NA_real_)
    x$fit$df.residual
})) %>% mutate(hc3_crit = qt(0.975, df = hc3_df), hc3_ci_low = hc3_coef - hc3_crit * hc3_se, hc3_ci_high = hc3_coef + 
    hc3_crit * hc3_se, hc3_fdr = p.adjust(hc3_p, method = "fdr")) %>% arrange(hc3_fdr, hc3_p)
secretion_annotations <- read.csv(ppi_input(ppi_path("annotations/metacardis_secretion_annotations.csv")), 
    check.names = FALSE)
metacardis_secretion_results <- merge(metacardis_secretion_results, secretion_annotations, by = "rxn", 
    sort = FALSE)
metacardis_secretion_results <- metacardis_secretion_results[order(metacardis_secretion_results$X), c("X", 
    "rxn", "metabolite_name", "hc3_coef", "hc3_se", "hc3_p", "hc3_df", "hc3_crit", "hc3_ci_low", "hc3_ci_high", 
    "hc3_fdr", "metabr")]
write.csv(metacardis_secretion_results, ppi_output(ppi_path("human_cohorts/metacardis_results/Metacardis_mgsec_fullFormula_totalMedicationAdjusted_noPPIExcludedByMetabolomics_preBoot.csv"), 
    name = "metacardis_predicted_secretion_hc3.csv"), row.names = FALSE)
