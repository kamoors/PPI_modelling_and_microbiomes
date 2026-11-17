# SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0
# Kiel final reaction abundance regression
# Human participant-level input matrices and covariates are not distributed.

source("scripts/helpers/setup.R")
invisible(run_analysis(
  section = "01_human_microbiome",
  directory = "human_cohorts/kiel",
  script = "01_04_kiel_reaction_abundance.R",
  restricted = TRUE,
  code = {
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
project_source(project_code("scripts/00_shared_helpers/00_01_load_kiel_metadata.R"))
project_setwd(project_path("human_cohorts/kiel/"))
library(lmtest)
library(sandwich)
library(MatrixGenerics)
library(compositions)
library(tibble)
reaction_abundance <- read.csv(project_input(project_path("human_cohorts/kiel/Reaction_abundance_mgpipe_focus.csv")))
colnames(reaction_abundance) <- c("rxn_name", stringr::str_sub(colnames(reaction_abundance[, -1]), 2))
stopifnot("rxn_name" %in% colnames(reaction_abundance))
missing_in_rxn <- setdiff(data$ID, colnames(reaction_abundance)[-1])
extra_in_rxn <- setdiff(colnames(reaction_abundance)[-1], data$ID)
if (length(missing_in_rxn) > 0) stop("Missing samples in reaction_abundance: ", paste(missing_in_rxn,
    collapse = ", "))
if (length(extra_in_rxn) > 0) message("Extra samples in reaction_abundance (will be dropped): ", paste(extra_in_rxn,
    collapse = ", "))
reaction_abundance <- reaction_abundance[, c("rxn_name", data$ID)]
stopifnot(all(colnames(reaction_abundance[, -1]) == data$ID))
rxn_abundance_dat <- t(reaction_abundance[, -1])
colnames(rxn_abundance_dat) <- reaction_abundance[, 1]
rxn_abundance_dat[is.na(rxn_abundance_dat)] <- 0
rxn_abundance_dat_filt <- rxn_abundance_dat[, colVars(rxn_abundance_dat) > 0]
rxn_total_abundance <- rowSums(rxn_abundance_dat_filt, na.rm = TRUE)
clr_pseudocount <- 1e-06
rxn_abundance_dat_filt <- log(rxn_abundance_dat_filt + clr_pseudocount)
rxn_abundance_dat_filt <- as.matrix(rxn_abundance_dat_filt)
focus_data <- cbind(data, phenotype_all_drugs_filt[, 4:16], phenoData_filt[, c(2, 8, 9, 18, 68)])
focus_data$glucose_mg.dl_mb <- as.numeric(focus_data$glucose_mg.dl_mb)
focus_data$diabetes_mb <- as.factor(focus_data$diabetes_mb)
focus_data$BMI.g25.30 <- as.factor(focus_data$BMI.g25.30)
focus_data$Analgesic <- as.factor(focus_data$Analgesic)
focus_data$Analgesic_NSAID <- as.factor(focus_data$Analgesic_NSAID)
focus_data$Antibiotics <- as.factor(focus_data$Antibiotics)
focus_data$Antihypertensives <- as.factor(focus_data$Antihypertensives)
focus_data$value2_systolic_mb <- as.numeric(focus_data$value2_systolic_mb)
focus_data$rxn_total_abundance <- rxn_total_abundance
focus_data$log_rxn_total_abundance <- log1p(rxn_total_abundance)
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
test_df <- focus_data
ppi_term <- "PPI1"
fit_all_reactions_hc3 <- function(data_in, test_df, fml) {
    out <- vector("list", ncol(data_in))
    for (i in seq_len(ncol(data_in))) {
        test_df$met <- data_in[, i]
        fit <- lm(fml, data = test_df)
        out[[i]] <- lmtest::coeftest(fit, vcov. = sandwich::vcovHC(fit, type = "HC3"))
    }
    names(out) <- colnames(data_in)
    out
}
project_source(project_code("scripts/00_shared_helpers/00_09_reaction_bootstrap_functions.R"))
fml_full <- met ~ PPI + Age + Gender + BMI + diabetes_mb + medication_count + Metformin + STATIN + ASPIRIN +
    log_rxn_total_abundance
fml_base <- met ~ PPI + Age + Gender + BMI
fml_boot_full <- fml_full
fml_boot_base <- fml_base
data_in <- rxn_abundance_dat_filt
test_df <- focus_data
react_ab_glm_full <- fit_all_reactions_hc3(data_in, test_df, fml_full)
focus_react_ab_results_preboot_full <- tibble(metabolite = names(react_ab_glm_full), hc3_coef = sapply(react_ab_glm_full,
    function(ct) {
        if (is.null(ct))
            return(NA_real_)
        if (ppi_term %in% rownames(ct))
            ct[ppi_term, "Estimate"]
        else NA_real_
    }), hc3_se = sapply(react_ab_glm_full, function(ct) {
    if (is.null(ct))
        return(NA_real_)
    if (ppi_term %in% rownames(ct))
        ct[ppi_term, "Std. Error"]
    else NA_real_
}), hc3_p = sapply(react_ab_glm_full, function(ct) {
    if (is.null(ct))
        return(NA_real_)
    if (ppi_term %in% rownames(ct))
        ct[ppi_term, "Pr(>|t|)"]
    else NA_real_
})) %>% mutate(hc3_ci_low = hc3_coef - qnorm(0.975) * hc3_se, hc3_ci_high = hc3_coef + qnorm(0.975) *
    hc3_se, hc3_fdr = p.adjust(hc3_p, method = "fdr")) %>% arrange(hc3_fdr, hc3_p)
write.csv(focus_react_ab_results_preboot_full, file = project_output(project_path("human_cohorts/results/Focus_mgpipe_reaction_abundance_LOG_HC3_FULL_withNewPPI_results.csv"),
    name = "kiel_reaction_abundance_hc3_full.csv"), row.names = FALSE)
project_source(project_code("scripts/00_shared_helpers/00_10_reaction_additional_analyses.R"))
}))
