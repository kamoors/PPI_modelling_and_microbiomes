# SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0
# MetaCardis MaAsLin2
# Run through run.R so input and output paths resolve within the package.
# Inputs: see data/input_manifest.csv and docs/RUN_ORDER.md.
# Human participant-level input matrices and covariates are not distributed.
# Supply authorised data separately; see docs/RESTRICTED_DATA.md.

library(readxl)
library(dplyr)
library(tibble)
library(Maaslin2)
library(matrixStats)
library(vegan)
library(ggplot2)
meta_dir <- ppi_path("human_cohorts/metacardis")
results_dir <- ppi_path("human_cohorts/results")
ppi_setwd(meta_dir)
cardis_meta <- readxl::read_excel(ppi_input("Merged_Metacardis-Study.xlsx"))
cardis_meta$ID <- trimws(cardis_meta$ID)
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
ppi_cols <- c("omeprazole", "pantoprazole", "dexlansoprazole")
ppi_cols <- ppi_cols[ppi_cols %in% colnames(met_matrix_all)]
if (length(ppi_cols) == 0) {
    stop("None of the expected PPI metabolite columns were found in metabolomics data.")
}
ppi_signal <- met_matrix_all[, ppi_cols, drop = FALSE]
ppi_baseline <- matrixStats::colMins(ppi_signal, na.rm = TRUE)
print(ppi_baseline)
ppi_detect_matrix <- sweep(ppi_signal, 2, ppi_baseline, FUN = function(x, base) x > base + 0.1)
ppi_detected <- rowSums(ppi_detect_matrix, na.rm = TRUE) > 0
ppi_measured <- ifelse(ppi_detected, "Yes", "No")
met_covariates_ppi <- cardis_meta_met %>% dplyr::select(ID, Proton.pump.inhibitor) %>% as.data.frame()
met_covariates_ppi$Proton.pump.inhibitor <- factor(met_covariates_ppi$Proton.pump.inhibitor, levels = c("No", 
    "Yes"))
met_covariates_ppi$PPI_measured <- factor(ppi_measured, levels = c("No", "Yes"))
keep_met_samples <- !(met_covariates_ppi$Proton.pump.inhibitor == "No" & met_covariates_ppi$PPI_measured == 
    "Yes")
met_covariates_clean <- met_covariates_ppi[keep_met_samples, , drop = FALSE]
message("Samples removed due to hidden PPI exposure: ", sum(!keep_met_samples))
print(table(Metadata = met_covariates_ppi$Proton.pump.inhibitor, Measured = met_covariates_ppi$PPI_measured, 
    useNA = "ifany"))
print(table(Metadata_clean = met_covariates_clean$Proton.pump.inhibitor, Measured_clean = met_covariates_clean$PPI_measured, 
    useNA = "ifany"))
clean_ids <- met_covariates_clean$ID
bact_ab <- cardis_meta[, 33:761]
bact_ab_num <- bact_ab %>% mutate(across(everything(), ~suppressWarnings(as.numeric(.))))
mg_matrix_all <- as.data.frame(bact_ab_num)
cardis_meta_mg <- cardis_meta %>% dplyr::filter(ID %in% clean_ids) %>% dplyr::arrange(match(ID, clean_ids))
mg_matrix_all <- mg_matrix_all[match(cardis_meta_mg$ID, cardis_meta$ID), , drop = FALSE]
stopifnot(nrow(cardis_meta_mg) == nrow(mg_matrix_all))
stopifnot(all(cardis_meta_mg$ID == clean_ids))
mg_feat_keep <- sapply(mg_matrix_all, function(x) {
    x <- as.numeric(x)
    v <- var(x, na.rm = TRUE)
    is.finite(v) && v > 0 && sum(x, na.rm = TRUE) > 0
})
mg_matrix_all <- mg_matrix_all[, mg_feat_keep, drop = FALSE]
message("Number of metagenomics features after initial filtering: ", ncol(mg_matrix_all))
med_vars <- c("Antidiabetic.drugs", "Insulin", "Metformin", "Sulfonurea", "DPP.IV.inhibitors", "GLP1.agonists", 
    "Anti-hypertensive.drugs", "Beta-blockers", "ACE", "AT2.inhibitor", "Calcium.blockers", "Anse.diuretics", 
    "Potassium-.sparing.diuretic", "Thiazides", "Lipid.lowering.drugs", "Statin", "Ezetimib", "Fibrate", 
    "Nitrate", "Amiodarone", "Digitalic", "Aspirine", "Clopidogrel", "Coumarin")
cardis_meta_mg[med_vars] <- lapply(cardis_meta_mg[med_vars], function(x) {
    x <- as.character(x)
    ifelse(x %in% c("Yes", "1"), 1, ifelse(x %in% c("No", "0"), 0, NA))
})
med_broad <- c("Antidiabetic.drugs", "Anti-hypertensive.drugs", "Lipid.lowering.drugs", "Nitrate", "Amiodarone", 
    "Digitalic", "Aspirine", "Clopidogrel", "Coumarin")
cardis_meta_mg$medication_count <- rowSums(cardis_meta_mg[med_broad], na.rm = TRUE)
mg_covariates_all <- cardis_meta_mg %>% dplyr::select(ID, Gender, Proton.pump.inhibitor, `Age.(years)`, 
    `BMI.(kg/m²)`, Nationality, Diabetic.status, Statin, Aspirine, Metformin, creatinine, medication_count) %>% 
    as.data.frame()
colnames(mg_covariates_all) <- make.names(colnames(mg_covariates_all), unique = TRUE)
mg_covariates_all$ID <- as.character(mg_covariates_all$ID)
mg_covariates_all$Gender <- as.factor(mg_covariates_all$Gender)
mg_covariates_all$Proton.pump.inhibitor <- factor(mg_covariates_all$Proton.pump.inhibitor, levels = c("No", 
    "Yes"))
mg_covariates_all$Age..years. <- suppressWarnings(as.numeric(mg_covariates_all$Age..years.))
mg_covariates_all$BMI..kg.m.. <- suppressWarnings(as.numeric(mg_covariates_all$BMI..kg.m..))
mg_covariates_all$Nationality <- as.factor(mg_covariates_all$Nationality)
mg_covariates_all$Systolic.blood.pressure..mmHg. <- suppressWarnings(as.numeric(mg_covariates_all$Systolic.blood.pressure..mmHg.))
mg_covariates_all$Diabetic.status <- as.factor(mg_covariates_all$Diabetic.status)
mg_covariates_all$Anti.hypertensive.drugs <- as.factor(mg_covariates_all$Anti.hypertensive.drugs)
mg_covariates_all$Statin <- as.factor(mg_covariates_all$Statin)
mg_covariates_all$Aspirine <- as.factor(mg_covariates_all$Aspirine)
mg_covariates_all$Metformin <- as.factor(mg_covariates_all$Metformin)
mg_covariates_all$creatinine <- suppressWarnings(as.numeric(mg_covariates_all$creatinine))
mg_covariates_all$medication_count <- suppressWarnings(as.numeric(mg_covariates_all$medication_count))
if (anyDuplicated(mg_covariates_all$ID) == 0) {
    rownames(mg_covariates_all) <- mg_covariates_all$ID
    rownames(mg_matrix_all) <- mg_covariates_all$ID
} else {
    stop("Duplicate sample IDs found in mg_covariates_all.")
}
mg_fixed_base <- c("Proton.pump.inhibitor", "Gender", "Age..years.", "BMI..kg.m..", "Nationality")
mg_fixed_full <- c("Proton.pump.inhibitor", "Gender", "Age..years.", "BMI..kg.m..", "Nationality", "Diabetic.status", 
    "medication_count", "Metformin", "Statin", "Aspirine", "creatinine")
filter_by_prevalence <- function(df, min_prev = 0.1) {
    keep <- sapply(df, function(x) mean(x > 0, na.rm = TRUE) >= min_prev)
    df[, keep, drop = FALSE]
}
filter_by_variance <- function(df) {
    keep <- sapply(df, function(x) {
        x <- as.numeric(x)
        v <- var(x, na.rm = TRUE)
        is.finite(v) && v > 0
    })
    df[, keep, drop = FALSE]
}
remove_na_samples_from_features <- function(data_df, meta_df) {
    na_samples <- apply(data_df, 1, function(x) any(is.na(x)))
    message("Samples with any feature NA: ", sum(na_samples))
    list(data = data_df[!na_samples, , drop = FALSE], meta = meta_df[!na_samples, , drop = FALSE], na_flag = na_samples)
}
mg_keep_base <- complete.cases(mg_covariates_all[, mg_fixed_base, drop = FALSE])
mg_meta_base <- mg_covariates_all[mg_keep_base, , drop = FALSE]
mg_data_base <- mg_matrix_all[mg_keep_base, , drop = FALSE]
message("BASE complete cases before feature-NA removal: ", nrow(mg_meta_base))
message("BASE PPI counts before feature-NA removal:")
print(table(mg_meta_base$Proton.pump.inhibitor, useNA = "ifany"))
mg_base_clean <- remove_na_samples_from_features(mg_data_base, mg_meta_base)
mg_data_base_clean <- mg_base_clean$data
mg_meta_base_clean <- mg_base_clean$meta
message("BASE samples after feature-NA removal: ", nrow(mg_meta_base_clean))
message("BASE PPI counts after feature-NA removal:")
print(table(mg_meta_base_clean$Proton.pump.inhibitor, useNA = "ifany"))
mg_data_base_final <- mg_data_base_clean %>% filter_by_prevalence(min_prev = 0.1) %>% filter_by_variance()
message("BASE features after prevalence/variance filtering: ", ncol(mg_data_base_final))
mg_keep_full <- complete.cases(mg_covariates_all[, mg_fixed_full, drop = FALSE])
mg_meta_full <- mg_covariates_all[mg_keep_full, , drop = FALSE]
mg_data_full <- mg_matrix_all[mg_keep_full, , drop = FALSE]
message("FULL complete cases before feature-NA removal: ", nrow(mg_meta_full))
message("FULL PPI counts before feature-NA removal:")
print(table(mg_meta_full$Proton.pump.inhibitor, useNA = "ifany"))
mg_full_clean <- remove_na_samples_from_features(mg_data_full, mg_meta_full)
mg_data_full_clean <- mg_full_clean$data
mg_meta_full_clean <- mg_full_clean$meta
message("FULL samples after feature-NA removal: ", nrow(mg_meta_full_clean))
message("FULL PPI counts after feature-NA removal:")
print(table(mg_meta_full_clean$Proton.pump.inhibitor, useNA = "ifany"))
mg_data_full_final <- mg_data_full_clean %>% filter_by_prevalence(min_prev = 0.1) %>% filter_by_variance()
message("FULL features after prevalence/variance filtering: ", ncol(mg_data_full_final))
stopifnot(identical(rownames(mg_data_base_final), rownames(mg_meta_base_clean)))
stopifnot(identical(rownames(mg_data_full_final), rownames(mg_meta_full_clean)))
if (ncol(mg_data_base_final) == 0) stop("No BASE features left after filtering.")
if (ncol(mg_data_full_final) == 0) stop("No FULL features left after filtering.")
dir.create(ppi_output(file.path(results_dir, "Maaslin2_BASE"), directory = TRUE), showWarnings = FALSE, 
    recursive = TRUE)
fit_base <- Maaslin2(input_data = mg_data_base_final, input_metadata = mg_meta_base_clean, output = file.path(results_dir, 
    "Maaslin2_BASE"), fixed_effects = mg_fixed_base, random_effects = NULL, reference = c("Proton.pump.inhibitor,No"), 
    normalization = "NONE", transform = "NONE", analysis_method = "LM", correction = "BH", standardize = FALSE, 
    min_abundance = 0, min_prevalence = 0, plot_heatmap = TRUE, plot_scatter = FALSE)
dir.create(ppi_output(file.path(results_dir, "Maaslin2_FULL"), directory = TRUE), showWarnings = FALSE, 
    recursive = TRUE)
fit_full <- Maaslin2(input_data = mg_data_full_final, input_metadata = mg_meta_full_clean, output = file.path(results_dir, 
    "Maaslin2_FULL"), fixed_effects = mg_fixed_full, random_effects = NULL, reference = c("Proton.pump.inhibitor,No"), 
    normalization = "NONE", transform = "NONE", analysis_method = "LM", correction = "BH", standardize = FALSE, 
    min_abundance = 0, min_prevalence = 0, plot_heatmap = TRUE, plot_scatter = FALSE)
mg_res_base <- read.delim(ppi_input(file.path(results_dir, "Maaslin2_BASE", "all_results.tsv")), check.names = FALSE)
mg_res_base_clean <- mg_res_base %>% filter(metadata == "Proton.pump.inhibitor", value == "Yes") %>% 
    transmute(feature, metadata, level = value, coef = coef, stderr = stderr, pval = pval, qval = qval) %>% 
    arrange(qval, pval)
write.csv(mg_res_base_clean, ppi_output(file.path(results_dir, "MetaCardis_metagenomics_MaAsLin2_BASE_clean_results.csv"), 
    name = "meta_cardis_metagenomics_ma_as_lin2_base_clean_results.csv"), row.names = FALSE)
write.csv(mg_res_base_clean %>% filter(qval < 0.1), ppi_output(file.path(results_dir, "MetaCardis_metagenomics_MaAsLin2_BASE_q_lt_0.1.csv"), 
    name = "meta_cardis_metagenomics_ma_as_lin2_base_q_lt_0_1.csv"), row.names = FALSE)
mg_res_full <- read.delim(ppi_input(file.path(results_dir, "Maaslin2_FULL", "all_results.tsv")), check.names = FALSE)
mg_res_full_clean <- mg_res_full %>% filter(metadata == "Proton.pump.inhibitor", value == "Yes") %>% 
    transmute(feature, metadata, level = value, coef = coef, stderr = stderr, pval = pval, qval = qval) %>% 
    arrange(qval, pval)
write.csv(mg_res_full_clean, ppi_output(file.path(results_dir, "MetaCardis_metagenomics_MaAsLin2_FULL_clean_results.csv"), 
    name = "meta_cardis_metagenomics_ma_as_lin2_full_clean_results.csv"), row.names = FALSE)
write.csv(mg_res_full_clean %>% filter(qval < 0.1), ppi_output(file.path(results_dir, "MetaCardis_metagenomics_MaAsLin2_FULL_q_lt_0.1.csv"), 
    name = "meta_cardis_metagenomics_ma_as_lin2_full_q_lt_0_1.csv"), row.names = FALSE)
message("BASE model: q < 0.2 = ", sum(mg_res_base_clean$qval < 0.2, na.rm = TRUE))
message("FULL model: q < 0.2 = ", sum(mg_res_full_clean$qval < 0.2, na.rm = TRUE))
message("BASE retained samples: ", nrow(mg_meta_base_clean))
message("FULL retained samples: ", nrow(mg_meta_full_clean))
message("BASE retained features: ", ncol(mg_data_base_final))
message("FULL retained features: ", ncol(mg_data_full_final))
saveRDS(list(fit_base = fit_base, mg_res_base = mg_res_base, mg_res_base_clean = mg_res_base_clean, mg_meta_base_clean = mg_meta_base_clean, 
    mg_data_base_final = mg_data_base_final), ppi_output(file.path(results_dir, "MetaCardis_metagenomics_MaAsLin2_BASE_results.rds"), 
    name = "meta_cardis_metagenomics_ma_as_lin2_base_results.rds"))
saveRDS(list(fit_full = fit_full, mg_res_full = mg_res_full, mg_res_full_clean = mg_res_full_clean, mg_meta_full_clean = mg_meta_full_clean, 
    mg_data_full_final = mg_data_full_final), ppi_output(file.path(results_dir, "MetaCardis_metagenomics_MaAsLin2_FULL_results.rds"), 
    name = "meta_cardis_metagenomics_ma_as_lin2_full_results.rds"))
extract_ppi_term <- function(model, metric_name, term = "Proton.pump.inhibitorYes") {
    sm <- summary(model)$coefficients
    if (!(term %in% rownames(sm))) {
        return(tibble(metric = metric_name, term = term, estimate = NA_real_, std_error = NA_real_, t_value = NA_real_, 
            p_value = NA_real_))
    }
    tibble(metric = metric_name, term = term, estimate = sm[term, "Estimate"], std_error = sm[term, "Std. Error"], 
        t_value = sm[term, "t value"], p_value = sm[term, "Pr(>|t|)"])
}
stopifnot(identical(rownames(mg_data_base_final), rownames(mg_meta_base_clean)))
alpha_df_base <- as.data.frame(mg_data_base_final)
alpha_df_base[] <- lapply(alpha_df_base, function(x) as.numeric(as.character(x)))
Observed_base <- rowSums(alpha_df_base > 0, na.rm = TRUE)
Shannon_base <- vegan::diversity(alpha_df_base, index = "shannon")
alpha_meta_base <- mg_meta_base_clean %>% mutate(Observed = Observed_base, Shannon = Shannon_base)
write.csv(alpha_meta_base, ppi_output(file.path(results_dir, "MetaCardis_alpha_diversity_BASE_table.csv"), 
    name = "meta_cardis_alpha_diversity_base_table.csv"), row.names = FALSE)
alpha_summary_base <- alpha_meta_base %>% group_by(Proton.pump.inhibitor) %>% summarise(n = n(), Observed_mean = mean(Observed, 
    na.rm = TRUE), Observed_sd = sd(Observed, na.rm = TRUE), Observed_median = median(Observed, na.rm = TRUE), 
    Shannon_mean = mean(Shannon, na.rm = TRUE), Shannon_sd = sd(Shannon, na.rm = TRUE), Shannon_median = median(Shannon, 
        na.rm = TRUE), .groups = "drop")
write.csv(alpha_summary_base, ppi_output(file.path(results_dir, "MetaCardis_alpha_diversity_BASE_summary.csv"), 
    name = "meta_cardis_alpha_diversity_base_summary.csv"), row.names = FALSE)
wilcox_observed_base <- wilcox.test(Observed ~ Proton.pump.inhibitor, data = alpha_meta_base)
wilcox_shannon_base <- wilcox.test(Shannon ~ Proton.pump.inhibitor, data = alpha_meta_base)
alpha_unadjusted_tests_base <- tibble(metric = c("Observed", "Shannon"), test = "Wilcoxon rank-sum", 
    p_value = c(wilcox_observed_base$p.value, wilcox_shannon_base$p.value))
write.csv(alpha_unadjusted_tests_base, ppi_output(file.path(results_dir, "MetaCardis_alpha_diversity_BASE_unadjusted_tests.csv"), 
    name = "meta_cardis_alpha_diversity_base_unadjusted_tests.csv"), row.names = FALSE)
lm_observed_base <- lm(Observed ~ Proton.pump.inhibitor + Gender + Age..years. + BMI..kg.m.. + Nationality, 
    data = alpha_meta_base)
lm_shannon_base <- lm(Shannon ~ Proton.pump.inhibitor + Gender + Age..years. + BMI..kg.m.. + Nationality, 
    data = alpha_meta_base)
alpha_adjusted_results_base <- bind_rows(extract_ppi_term(lm_observed_base, "Observed"), extract_ppi_term(lm_shannon_base, 
    "Shannon")) %>% mutate(fdr = p.adjust(p_value, method = "fdr"))
write.csv(alpha_adjusted_results_base, ppi_output(file.path(results_dir, "MetaCardis_alpha_diversity_BASE_adjusted_results.csv"), 
    name = "meta_cardis_alpha_diversity_base_adjusted_results.csv"), row.names = FALSE)
stopifnot(identical(rownames(mg_data_full_final), rownames(mg_meta_full_clean)))
alpha_df_full <- as.data.frame(mg_data_full_final)
alpha_df_full[] <- lapply(alpha_df_full, function(x) as.numeric(as.character(x)))
Observed_full <- rowSums(alpha_df_full > 0, na.rm = TRUE)
Shannon_full <- vegan::diversity(alpha_df_full, index = "shannon")
alpha_meta_full <- mg_meta_full_clean %>% mutate(Observed = Observed_full, Shannon = Shannon_full)
write.csv(alpha_meta_full, ppi_output(file.path(results_dir, "MetaCardis_alpha_diversity_FULL_table.csv"), 
    name = "meta_cardis_alpha_diversity_full_table.csv"), row.names = FALSE)
alpha_summary_full <- alpha_meta_full %>% group_by(Proton.pump.inhibitor) %>% summarise(n = n(), Observed_mean = mean(Observed, 
    na.rm = TRUE), Observed_sd = sd(Observed, na.rm = TRUE), Observed_median = median(Observed, na.rm = TRUE), 
    Shannon_mean = mean(Shannon, na.rm = TRUE), Shannon_sd = sd(Shannon, na.rm = TRUE), Shannon_median = median(Shannon, 
        na.rm = TRUE), .groups = "drop")
write.csv(alpha_summary_full, ppi_output(file.path(results_dir, "MetaCardis_alpha_diversity_FULL_summary.csv"), 
    name = "meta_cardis_alpha_diversity_full_summary.csv"), row.names = FALSE)
wilcox_observed_full <- wilcox.test(Observed ~ Proton.pump.inhibitor, data = alpha_meta_full)
wilcox_shannon_full <- wilcox.test(Shannon ~ Proton.pump.inhibitor, data = alpha_meta_full)
alpha_unadjusted_tests_full <- tibble(metric = c("Observed_FULL", "Shannon_FULL"), test = "Wilcoxon rank-sum", 
    p_value = c(wilcox_observed_full$p.value, wilcox_shannon_full$p.value))
write.csv(alpha_unadjusted_tests_full, ppi_output(file.path(results_dir, "MetaCardis_alpha_diversity_FULL_unadjusted_tests.csv"), 
    name = "meta_cardis_alpha_diversity_full_unadjusted_tests.csv"), row.names = FALSE)
lm_observed_full <- lm(Observed ~ Proton.pump.inhibitor + Gender + Age..years. + BMI..kg.m.. + Nationality + 
    Systolic.blood.pressure..mmHg. + Diabetic.status + medication_count + Metformin + Statin + Aspirine + 
    creatinine, data = alpha_meta_full)
lm_shannon_full <- lm(Shannon ~ Proton.pump.inhibitor + Gender + Age..years. + BMI..kg.m.. + Nationality + 
    Systolic.blood.pressure..mmHg. + Diabetic.status + medication_count + Metformin + Statin + Aspirine + 
    creatinine, data = alpha_meta_full)
alpha_adjusted_results_full <- bind_rows(extract_ppi_term(lm_observed_full, "Observed_FULL"), extract_ppi_term(lm_shannon_full, 
    "Shannon_FULL")) %>% mutate(fdr = p.adjust(p_value, method = "fdr"))
write.csv(alpha_adjusted_results_full, ppi_output(file.path(results_dir, "MetaCardis_alpha_diversity_FULL_adjusted_results.csv"), 
    name = "meta_cardis_alpha_diversity_full_adjusted_results.csv"), row.names = FALSE)
p_observed <- ggplot(alpha_meta_full, aes(x = Proton.pump.inhibitor, y = Observed)) + geom_boxplot(outlier.shape = NA) + 
    geom_jitter(width = 0.15, alpha = 0.35) + theme_classic() + labs(title = "Observed richness by PPI", 
    x = "Proton pump inhibitor", y = "Observed richness")
ggsave(filename = ppi_output(file.path(results_dir, "MetaCardis_alpha_diversity_FULL_Observed_boxplot.png"), 
    name = "meta_cardis_alpha_diversity_full_observed_boxplot.png"), plot = p_observed, width = 5, height = 4, 
    dpi = 300)
p_shannon <- ggplot(alpha_meta_full, aes(x = Proton.pump.inhibitor, y = Shannon)) + geom_boxplot(outlier.shape = NA) + 
    geom_jitter(width = 0.15, alpha = 0.35) + theme_classic() + labs(title = "Shannon diversity by PPI", 
    x = "Proton pump inhibitor", y = "Shannon diversity")
ggsave(filename = ppi_output(file.path(results_dir, "MetaCardis_alpha_diversity_FULL_Shannon_boxplot.png"), 
    name = "meta_cardis_alpha_diversity_full_shannon_boxplot.png"), plot = p_shannon, width = 5, height = 4, 
    dpi = 300)
message("Alpha diversity analysis complete.")
print(alpha_summary_full)
print(alpha_adjusted_results_full)
