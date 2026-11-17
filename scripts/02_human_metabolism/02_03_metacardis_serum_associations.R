# SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0
# MetaCardis serum-metabolomics HC3 analysis
# Human participant-level input matrices and covariates are not distributed.

source("scripts/helpers/setup.R")
invisible(run_analysis(
  section = "02_human_metabolism",
  directory = "human_cohorts/metacardis",
  script = "02_03_metacardis_serum_associations.R",
  restricted = TRUE,
  code = {
library(lmtest)
library(sandwich)
library(gplots)
library(ggfortify)
library(ggplot2)
library(readxl)
library(writexl)
library(reshape2)
library(ggsignif)
library(dplyr)
library(stringr)
library(RColorBrewer)
library(car)
library(boot.pval)
library(MatrixGenerics)
library(tibble)
project_setwd(project_path("human_cohorts/metacardis"))
cardis_meta <- readxl::read_excel(project_input("Merged_Metacardis-Study.xlsx"))
med_vars <- c("Antidiabetic.drugs", "Insulin", "Metformin", "Sulfonurea", "DPP.IV.inhibitors", "GLP1.agonists",
    "Anti-hypertensive.drugs", "Beta-blockers", "ACE", "AT2.inhibitor", "Calcium.blockers", "Anse.diuretics",
    "Potassium-.sparing.diuretic", "Thiazides", "Lipid.lowering.drugs", "Statin", "Ezetimib", "Fibrate",
    "Nitrate", "Amiodarone", "Digitalic", "Aspirine", "Clopidogrel", "Coumarin")
cardis_meta[med_vars] <- lapply(cardis_meta[med_vars], function(x) {
    x <- as.character(x)
    ifelse(x %in% c("Yes", "1"), 1, ifelse(x %in% c("No", "0"), 0, NA))
})
med_broad <- c("Antidiabetic.drugs", "Anti-hypertensive.drugs", "Lipid.lowering.drugs", "Nitrate", "Amiodarone",
    "Digitalic", "Aspirine", "Clopidogrel", "Coumarin")
cardis_meta$medication_count <- rowSums(cardis_meta[med_broad], na.rm = TRUE)
met_data_raw <- cardis_meta[, 791:(ncol(cardis_meta) - 3)]
met_data_num <- met_data_raw %>% mutate(across(everything(), ~suppressWarnings(as.numeric(.))))
met_matrix_all <- as.matrix(met_data_num)
met_all_na_rows <- rowSums(is.na(met_matrix_all)) == ncol(met_matrix_all)
message("Samples with all metabolomics values missing: ", sum(met_all_na_rows))
if (sum(met_all_na_rows) > 0) {
    message("First all-NA sample IDs: ", paste(head(cardis_meta$ID[met_all_na_rows], 20), collapse = ", "))
}
cardis_meta_met <- cardis_meta[!met_all_na_rows, , drop = FALSE]
met_matrix_all <- met_matrix_all[!met_all_na_rows, , drop = FALSE]
stopifnot(nrow(cardis_meta_met) == nrow(met_matrix_all))
met_var_keep <- MatrixGenerics::colVars(met_matrix_all, na.rm = TRUE) > 0
met_matrix_all <- met_matrix_all[, met_var_keep, drop = FALSE]
message("Samples retained after removing all-NA metabolomics rows: ", nrow(met_matrix_all))
message("Number of metabolomics features after variance filter: ", ncol(met_matrix_all))
ppi_cols <- c("omeprazole", "pantoprazole", "dexlansoprazole")
ppi_signal <- met_matrix_all[, ppi_cols, drop = FALSE]
ppi_baseline <- colMins(ppi_signal, na.rm = TRUE)
ppi_baseline
ppi_detect_matrix <- sweep(ppi_signal, 2, ppi_baseline, FUN = function(x, base) x > base + 0.1)
ppi_detected <- rowSums(ppi_detect_matrix, na.rm = TRUE) > 0
ppi_metabolomics <- ifelse(ppi_detected, "Yes", "No")
met_covariates_all <- cardis_meta_met %>% dplyr::select(ID, Gender, Proton.pump.inhibitor, `Age.(years)`,
    `BMI.(kg/m²)`, Nationality, `Systolic.blood.pressure.(mmHg)`, Diabetic.status, `Anti-hypertensive.drugs`,
    Statin, Aspirine, Metformin, creatinine, medication_count, SMOKE) %>% as.data.frame()
met_covariates_all$PPI_measured <- factor(ppi_metabolomics)
met_covariates_all$Gender <- as.factor(met_covariates_all$Gender)
met_covariates_all$Proton.pump.inhibitor <- as.factor(met_covariates_all$Proton.pump.inhibitor)
met_covariates_all$`Age.(years)` <- suppressWarnings(as.numeric(met_covariates_all$`Age.(years)`))
met_covariates_all$`BMI.(kg/m²)` <- suppressWarnings(as.numeric(met_covariates_all$`BMI.(kg/m²)`))
met_covariates_all$Nationality <- as.factor(met_covariates_all$Nationality)
met_covariates_all$`Systolic.blood.pressure.(mmHg)` <- suppressWarnings(as.numeric(met_covariates_all$`Systolic.blood.pressure.(mmHg)`))
met_covariates_all$Diabetic.status <- as.factor(met_covariates_all$Diabetic.status)
met_covariates_all$`Anti-hypertensive.drugs` <- as.factor(met_covariates_all$`Anti-hypertensive.drugs`)
met_covariates_all$Statin <- as.factor(met_covariates_all$Statin)
met_covariates_all$Aspirine <- as.factor(met_covariates_all$Aspirine)
met_covariates_all$Metformin <- as.factor(met_covariates_all$Metformin)
met_covariates_all$creatinine <- suppressWarnings(as.numeric(met_covariates_all$creatinine))
met_covariates_all$medication_count <- suppressWarnings(as.numeric(met_covariates_all$medication_count))
met_covariates_all$SMOKE <- as.factor(met_covariates_all$SMOKE)
print(table(met_covariates_all$PPI_measured, useNA = "ifany"))
if (all(c("No", "Yes") %in% unique(na.omit(met_covariates_all$Proton.pump.inhibitor)))) {
    met_covariates_all$Proton.pump.inhibitor <- factor(met_covariates_all$Proton.pump.inhibitor, levels = c("No",
        "Yes"))
}
ppi_meta <- met_covariates_all$Proton.pump.inhibitor
keep_samples <- !(ppi_meta == "No" & ppi_detected)
met_covariates_clean <- met_covariates_all[keep_samples, , drop = FALSE]
met_matrix_clean <- met_matrix_all[keep_samples, , drop = FALSE]
met_covariates_clean$Proton.pump.inhibitor <- factor(met_covariates_clean$Proton.pump.inhibitor, levels = c("No",
    "Yes"))
met_fml_base <- met ~ Proton.pump.inhibitor + Gender + `Age.(years)` + `BMI.(kg/m²)` + Nationality
met_fml_full <- met ~ Proton.pump.inhibitor + Gender + `Age.(years)` + `BMI.(kg/m²)` + Nationality +
    Diabetic.status + medication_count + Metformin + Statin + Aspirine + creatinine
met_vars_base <- all.vars(met_fml_base)
met_vars_full <- all.vars(met_fml_full)
met_keep_base <- complete.cases(met_covariates_clean[, met_vars_base[-1], drop = FALSE])
met_keep_full <- complete.cases(met_covariates_clean[, met_vars_full[-1], drop = FALSE])
met_test_df_base <- met_covariates_clean[met_keep_base, , drop = FALSE]
met_test_df_full <- met_covariates_clean[met_keep_full, , drop = FALSE]
met_matrix_base <- met_matrix_clean[met_keep_base, , drop = FALSE]
met_matrix_full <- met_matrix_clean[met_keep_full, , drop = FALSE]
stopifnot(nrow(met_test_df_base) == nrow(met_matrix_base))
stopifnot(nrow(met_test_df_full) == nrow(met_matrix_full))
message("BASE complete cases: ", nrow(met_test_df_base))
message("FULL complete cases: ", nrow(met_test_df_full))
met_tmp_df <- met_test_df_base
met_tmp_df$met <- met_matrix_base[, 1]
met_tmp_fit <- lm(met_fml_base, data = met_tmp_df)
met_tmp_ct <- lmtest::coeftest(met_tmp_fit, vcov. = sandwich::vcovHC(met_tmp_fit, type = "HC3"))
print(rownames(met_tmp_ct))
met_ppi_term <- "Proton.pump.inhibitorYes"
if (!(met_ppi_term %in% rownames(met_tmp_ct))) {
    stop("met_ppi_term not found in coefficient names. Update met_ppi_term to one of: ", paste(rownames(met_tmp_ct),
        collapse = ", "))
}
met_fit_all_features_hc3 <- function(data_in, test_df, fml) {
    out <- vector("list", ncol(data_in))
    for (i in seq_len(ncol(data_in))) {
        test_df$met <- data_in[, i]
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
met_hc3_boot_one <- function(test_df, formula, ppi_term, R = 1000) {
    n <- nrow(test_df)
    beta_star <- rep(NA_real_, R)
    for (b in seq_len(R)) {
        idx <- sample.int(n, replace = TRUE)
        dboot <- test_df[idx, , drop = FALSE]
        fit <- tryCatch(lm(formula, data = dboot), error = function(e) NULL)
        if (is.null(fit))
            next
        ct <- tryCatch(lmtest::coeftest(fit, vcov. = sandwich::vcovHC(fit, type = "HC3")), error = function(e) NULL)
        if (is.null(ct))
            next
        if (!(ppi_term %in% rownames(ct)))
            next
        beta_star[b] <- ct[ppi_term, "Estimate"]
    }
    beta_star <- beta_star[is.finite(beta_star)]
    res <- data.frame(Estimate = NA_real_, p.value = NA_real_)
    rownames(res) <- ppi_term
    if (length(beta_star) == 0) {
        attr(res, "boot_beta") <- numeric(0)
        return(res)
    }
    res[ppi_term, "Estimate"] <- median(beta_star, na.rm = TRUE)
    res[ppi_term, "p.value"] <- 2 * min(mean(beta_star <= 0, na.rm = TRUE), mean(beta_star >= 0, na.rm = TRUE))
    attr(res, "boot_beta") <- beta_star
    res
}
met_bootstrap_subset <- function(data_in, test_df, fml_boot, ppi_term, R = 1000) {
    out <- vector("list", ncol(data_in))
    for (i in seq_len(ncol(data_in))) {
        test_df$met <- data_in[, i]
        out[[i]] <- met_hc3_boot_one(test_df = test_df, formula = fml_boot, ppi_term = ppi_term, R = R)
    }
    names(out) <- colnames(data_in)
    out
}
met_build_tables <- function(hc3_list, boot_list, ppi_term, hc3_fdr_cut = 0.2, boot_p_cut = 0.05) {
    met_hc3_tbl <- tibble(metabolite = names(hc3_list), hc3_coef = sapply(hc3_list, function(ct) {
        if (is.null(ct))
            return(NA_real_)
        if (ppi_term %in% rownames(ct))
            ct[ppi_term, "Estimate"]
        else NA_real_
    }), hc3_p = sapply(hc3_list, function(ct) {
        if (is.null(ct))
            return(NA_real_)
        if (ppi_term %in% rownames(ct))
            ct[ppi_term, "Pr(>|t|)"]
        else NA_real_
    })) %>% mutate(hc3_fdr = p.adjust(hc3_p, method = "fdr"))
    met_boot_tbl <- tibble(metabolite = names(boot_list), boot_coef_median = sapply(boot_list, function(x) x[ppi_term,
        "Estimate"]), boot_p = sapply(boot_list, function(x) x[ppi_term, "p.value"]), boot_ci_low = sapply(boot_list,
        function(x) {
            bb <- attr(x, "boot_beta")
            if (is.null(bb) || length(bb) == 0)
                NA_real_
            else as.numeric(quantile(bb, 0.025, na.rm = TRUE))
        }), boot_ci_high = sapply(boot_list, function(x) {
        bb <- attr(x, "boot_beta")
        if (is.null(bb) || length(bb) == 0)
            NA_real_
        else as.numeric(quantile(bb, 0.975, na.rm = TRUE))
    }))
    met_res_clean <- met_hc3_tbl %>% left_join(met_boot_tbl, by = "metabolite") %>% mutate(hc3_sig_fdr = hc3_fdr <
        hc3_fdr_cut, boot_sig_p = !is.na(boot_p) & boot_p < boot_p_cut, sig_both = hc3_sig_fdr & boot_sig_p) %>%
        arrange(desc(sig_both), hc3_fdr, hc3_p)
    met_res_clean_report <- met_res_clean %>% select(metabolite, hc3_coef, hc3_p, hc3_fdr, boot_coef_median,
        boot_ci_low, boot_ci_high, boot_p, sig_both)
    list(met_hc3_tbl = met_hc3_tbl, met_boot_tbl = met_boot_tbl, met_res_clean_report = met_res_clean_report)
}
met_glm_full <- met_fit_all_features_hc3(data_in = met_matrix_full, test_df = met_test_df_full, fml = met_fml_full)
met_results_preboot_full <- tibble(metabolite = names(met_glm_full), hc3_coef = sapply(met_glm_full,
    function(x) {
        if (is.null(x) || is.null(x$ct))
            return(NA_real_)
        ct <- x$ct
        if (met_ppi_term %in% rownames(ct))
            ct[met_ppi_term, "Estimate"]
        else NA_real_
    }), hc3_se = sapply(met_glm_full, function(x) {
    if (is.null(x) || is.null(x$ct))
        return(NA_real_)
    ct <- x$ct
    if (met_ppi_term %in% rownames(ct))
        ct[met_ppi_term, "Std. Error"]
    else NA_real_
}), hc3_p = sapply(met_glm_full, function(x) {
    if (is.null(x) || is.null(x$ct))
        return(NA_real_)
    ct <- x$ct
    if (met_ppi_term %in% rownames(ct))
        ct[met_ppi_term, "Pr(>|t|)"]
    else NA_real_
}), hc3_df = sapply(met_glm_full, function(x) {
    if (is.null(x) || is.null(x$fit))
        return(NA_real_)
    x$fit$df.residual
})) %>% mutate(hc3_crit = qt(0.975, df = hc3_df), hc3_ci_low = hc3_coef - hc3_crit * hc3_se, hc3_ci_high = hc3_coef +
    hc3_crit * hc3_se, hc3_fdr = p.adjust(hc3_p, method = "fdr")) %>% arrange(hc3_fdr, hc3_p)
write.csv(met_results_preboot_full, file = project_output(project_path("human_cohorts/results/Metacardis_metabolomics_fullFormula_totalMedicationAdjusted_noPPIExcludedByMetabolomics_preBoot.csv"),
    name = "metacardis_metabolomics_full_formula_total_medication_adjusted_no_ppi_excluded_by_metabolomics_pre_boot.csv"),
    row.names = FALSE)
met_hc3_p_full <- sapply(met_glm_full, function(ct) {
    if (is.null(ct))
        return(NA_real_)
    if (met_ppi_term %in% rownames(ct))
        ct[met_ppi_term, "Pr(>|t|)"]
    else NA_real_
})
met_hc3_fdr_full <- p.adjust(met_hc3_p_full, method = "fdr")
met_features_to_boot_full <- names(met_hc3_fdr_full)[is.finite(met_hc3_fdr_full) & met_hc3_fdr_full <
    0.2]
met_matrix_full_sig <- met_matrix_full[, met_features_to_boot_full, drop = FALSE]
met_sig_boot_full <- met_bootstrap_subset(data_in = met_matrix_full_sig, test_df = met_test_df_full,
    fml_boot = met_fml_full, ppi_term = met_ppi_term, R = 1000)
met_tables_full <- met_build_tables(hc3_list = met_glm_full, boot_list = met_sig_boot_full, ppi_term = met_ppi_term,
    hc3_fdr_cut = 0.2, boot_p_cut = 0.05)
met_res_list_full <- list(Metabolomics_HC3_full = met_glm_full, Metabolomics_preboot_full = met_results_preboot_full,
    Sig_metabolomics_boot_full_output = met_sig_boot_full, tables_full = met_tables_full)
saveRDS(met_res_list_full, project_output("MetaCardis_metabolomics_HC3_FULL_results.rds", name = "meta_cardis_metabolomics_hc3_full_results.rds"))
write.csv(met_results_preboot_full, file = project_output("MetaCardis_metabolomics_HC3_FULL_preboot_full_results.csv",
    name = "meta_cardis_metabolomics_hc3_full_preboot_full_results.csv"), row.names = FALSE)
write.csv(met_tables_full$met_res_clean_report, file = project_output("MetaCardis_metabolomics_HC3_FULL_clean_results.csv",
    name = "meta_cardis_metabolomics_hc3_full_clean_results.csv"), row.names = FALSE)
write.csv(met_tables_full$met_res_clean_report %>% filter(sig_both), file = project_output("MetaCardis_metabolomics_HC3_FULL_sig_both.csv",
    name = "meta_cardis_metabolomics_hc3_full_sig_both.csv"), row.names = FALSE)
met_glm_base <- met_fit_all_features_hc3(data_in = met_matrix_base, test_df = met_test_df_base, fml = met_fml_base)
met_results_preboot_base <- tibble(metabolite = names(met_glm_base), hc3_coef = sapply(met_glm_base,
    function(ct) {
        if (is.null(ct))
            return(NA_real_)
        if (met_ppi_term %in% rownames(ct))
            ct[met_ppi_term, "Estimate"]
        else NA_real_
    }), hc3_p = sapply(met_glm_base, function(ct) {
    if (is.null(ct))
        return(NA_real_)
    if (met_ppi_term %in% rownames(ct))
        ct[met_ppi_term, "Pr(>|t|)"]
    else NA_real_
})) %>% mutate(hc3_fdr = p.adjust(hc3_p, method = "fdr")) %>% arrange(hc3_fdr, hc3_p)
met_hc3_p_base <- sapply(met_glm_base, function(ct) {
    if (is.null(ct))
        return(NA_real_)
    if (met_ppi_term %in% rownames(ct))
        ct[met_ppi_term, "Pr(>|t|)"]
    else NA_real_
})
met_hc3_fdr_base <- p.adjust(met_hc3_p_base, method = "fdr")
met_features_to_boot_base <- names(met_hc3_fdr_base)[is.finite(met_hc3_fdr_base) & met_hc3_fdr_base <
    0.1]
met_matrix_base_sig <- met_matrix_base[, met_features_to_boot_base, drop = FALSE]
met_sig_boot_base <- met_bootstrap_subset(data_in = met_matrix_base_sig, test_df = met_test_df_base,
    fml_boot = met_fml_base, ppi_term = met_ppi_term, R = 1000)
met_tables_base <- met_build_tables(hc3_list = met_glm_base, boot_list = met_sig_boot_base, ppi_term = met_ppi_term,
    hc3_fdr_cut = 0.2, boot_p_cut = 0.05)
met_res_list_base <- list(Metabolomics_HC3_base = met_glm_base, Metabolomics_preboot_base = met_results_preboot_base,
    Sig_metabolomics_boot_base_output = met_sig_boot_base, tables_base = met_tables_base)
saveRDS(met_res_list_base, project_output("MetaCardis_metabolomics_HC3_BASE_results.rds", name = "meta_cardis_metabolomics_hc3_base_results.rds"))
write.csv(met_results_preboot_base, file = project_output("MetaCardis_metabolomics_HC3_BASE_preboot_full_results.csv",
    name = "meta_cardis_metabolomics_hc3_base_preboot_full_results.csv"), row.names = FALSE)
write.csv(met_tables_base$met_res_clean_report, file = project_output("MetaCardis_metabolomics_HC3_BASE_clean_results.csv",
    name = "meta_cardis_metabolomics_hc3_base_clean_results.csv"), row.names = FALSE)
write.csv(met_tables_base$met_res_clean_report %>% filter(sig_both), file = project_output("MetaCardis_metabolomics_HC3_BASE_sig_both.csv",
    name = "meta_cardis_metabolomics_hc3_base_sig_both.csv"), row.names = FALSE)
message("BASE model: bootstrapped metabolites = ", length(met_features_to_boot_base))
message("FULL model: bootstrapped metabolites = ", length(met_features_to_boot_full))
message("BASE model: sig_both = ", sum(met_tables_base$met_res_clean_report$sig_both, na.rm = TRUE))
message("FULL model: sig_both = ", sum(met_tables_full$met_res_clean_report$sig_both, na.rm = TRUE))
top20_met <- met_results_preboot_full %>% filter(!is.na(hc3_fdr)) %>% arrange(hc3_fdr) %>% slice(1:20) %>%
    pull(metabolite)
met_top20_matrix <- met_matrix_full[, top20_met, drop = FALSE]
met_top20_matrix <- met_top20_matrix[, colSums(is.na(met_top20_matrix)) == 0, drop = FALSE]
met_top20_scaled <- scale(met_top20_matrix)
all_metabolites <- scale(met_matrix_full)
annotation_df <- data.frame(PPI = met_test_df_full$Proton.pump.inhibitor)
rownames(annotation_df) <- met_test_df_full$ID
rownames(met_top20_scaled) <- met_test_df_full$ID
rownames(all_metabolites) <- met_test_df_full$ID
library(pheatmap)
pheatmap(t(all_metabolites), annotation_col = annotation_df, show_colnames = FALSE, fontsize_row = 8,
    clustering_method = "complete", main = "Top 20 PPI-associated metabolites (FULL model)")
library(dplyr)
library(pheatmap)
met_top20_res <- met_results_preboot_full %>% filter(!is.na(hc3_fdr), !is.na(hc3_coef)) %>% filter(hc3_fdr <
    0.1) %>% arrange(hc3_fdr, hc3_p) %>% arrange(hc3_coef)
met_top20_hc3_mat <- as.matrix(met_top20_res$hc3_coef)
rownames(met_top20_hc3_mat) <- met_top20_res$metabolite
colnames(met_top20_hc3_mat) <- "Metacardis Serum"
pheatmap(met_top20_hc3_mat, cluster_rows = FALSE, cluster_cols = FALSE, color = colorRampPalette(c("blue",
    "white", "red"))(100), breaks = seq(-max(abs(met_top20_hc3_mat), na.rm = TRUE), max(abs(met_top20_hc3_mat),
    na.rm = TRUE), length.out = 101), border_color = "grey80", fontsize_row = 9, main = "Metacardis Serum: Top 20 metabolites")
pheatmap(t(met_top20_hc3_mat), cluster_rows = FALSE, cluster_cols = FALSE, show_rownames = FALSE, fontsize_col = 14,
    angle_col = 45, color = colorRampPalette(c("blue", "white", "red"))(100), breaks = seq(-max(abs(met_top20_hc3_mat),
        na.rm = TRUE), max(abs(met_top20_hc3_mat), na.rm = TRUE), length.out = 101), border_color = "grey80",
    main = "Metacardis Serum: Top metabolites", cellheight = 15, cellwidth = 17)
}))
