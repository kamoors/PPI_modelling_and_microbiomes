# SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0
# Kiel ALDEx2
# Human participant-level input matrices and covariates are not distributed.

source("scripts/helpers/setup.R")
invisible(run_analysis(
  section = "01_human_microbiome",
  directory = "",
  script = "01_01_kiel_taxonomy.R",
  restricted = TRUE,
  code = {
library(dplyr)
library(readr)
library(stringr)
library(ALDEx2)
library(tibble)
project_source(project_code("scripts/00_shared_helpers/00_01_load_kiel_metadata.R"))
results_dir <- project_path("Results/Aldex2")
dir.create(project_output(results_dir, directory = TRUE), showWarnings = FALSE, recursive = TRUE)
focus_data <- cbind(data, phenotype_all_drugs_filt[, 4:16], phenoData_filt[, c(2, 8, 9, 18, 68)]) %>%
    as.data.frame()
med_broad <- c("Analgesic", "Antibiotics", "Antidepressant", "Antihypertensives", "Antiphlogistic", "Antidiabetic")
focus_data[med_broad] <- phenotype_all_drugs_filt[, med_broad]
focus_data[med_broad] <- lapply(focus_data[med_broad], function(x) {
    as.numeric(as.character(x))
})
focus_data$medication_count <- rowSums(focus_data[, med_broad, drop = FALSE], na.rm = TRUE)
focus_data$PPI <- as.factor(focus_data$PPI)
focus_data$Gender <- as.factor(focus_data$Gender)
focus_data$diabetes_mb <- as.factor(focus_data$diabetes_mb)
focus_data$Metformin <- as.factor(focus_data$Metformin)
focus_data$STATIN <- as.factor(focus_data$STATIN)
focus_data$ASPIRIN <- as.factor(focus_data$ASPIRIN)
focus_data$Age <- suppressWarnings(as.numeric(focus_data$Age))
focus_data$BMI <- suppressWarnings(as.numeric(focus_data$BMI))
focus_data$medication_count <- suppressWarnings(as.numeric(focus_data$medication_count))
full_formula <- ~PPI + Age + Gender + BMI + diabetes_mb + medication_count + Metformin + STATIN + ASPIRIN
model_vars <- all.vars(full_formula)
keep_samples <- complete.cases(focus_data[, model_vars, drop = FALSE])
focus_data_full <- focus_data[keep_samples, , drop = FALSE]
message("Kiel complete cases retained: ", nrow(focus_data_full))
print(table(focus_data_full$PPI, useNA = "ifany"))
mm_full <- model.matrix(full_formula, focus_data_full)
run_aldex_full <- function(counts, covariates, model_matrix, feature_level, mc_samples = 256, ppi_term = "PPI1") {
    counts <- as.data.frame(counts)
    counts[] <- lapply(counts, function(x) {
        suppressWarnings(as.numeric(as.character(x)))
    })
    counts[is.na(counts)] <- 0
    counts <- counts[, rownames(covariates), drop = FALSE]
    stopifnot(all(colnames(counts) == rownames(covariates)))
    counts <- counts[rowSums(counts, na.rm = TRUE) > 0, , drop = FALSE]
    message(feature_level, " features retained after non-zero filtering: ", nrow(counts))
    x <- aldex.clr(counts, model_matrix, mc.samples = mc_samples, denom = "all", verbose = FALSE)
    glm_test <- aldex.glm(x, model_matrix)
    glm_effect <- aldex.glm.effect(x)
    glm_test_out <- glm_test %>% as.data.frame() %>% rownames_to_column("feature")
    if (ppi_term %in% names(glm_effect)) {
        effect_out <- glm_effect[[ppi_term]] %>% as.data.frame() %>% rownames_to_column("feature")
    }
    else {
        effect_out <- data.frame(feature = rownames(glm_test))
    }
    clean_out <- glm_test_out %>% select(feature, any_of(c(paste0(ppi_term, ":Est"), paste0(ppi_term,
        ":SE"), paste0(ppi_term, ":t.val"), paste0(ppi_term, ":pval"), paste0(ppi_term, ":pval.holm"),
        paste0(ppi_term, ":pval.padj")))) %>% rename_with(~c("feature", "coef", "se", "t_value", "pval",
        "pval_holm", "qval")[seq_along(.x)])
    list(clr = x, glm_test = glm_test, glm_effect = glm_effect, glm_test_out = glm_test_out, effect_out = effect_out,
        clean_out = clean_out, model_matrix = model_matrix, formula = full_formula, feature_level = feature_level)
}
species_counts <- counts_filt
species_counts <- species_counts[, keep_samples, drop = FALSE]
colnames(species_counts) <- focus_data_full$ID
rownames(focus_data_full) <- focus_data_full$ID
species_results <- run_aldex_full(counts = species_counts, covariates = focus_data_full, model_matrix = mm_full,
    feature_level = "species", mc_samples = 256, ppi_term = "PPI1")
write.csv(species_results$glm_test_out, project_output(file.path(results_dir, "Kiel_species_ALDEx2_FULL_all_results.csv"),
    name = "kiel_species_ald_ex2_full_all_results.csv"), row.names = FALSE)
write.csv(species_results$effect_out, project_output(file.path(results_dir, "Kiel_species_ALDEx2_FULL_effects.csv"),
    name = "kiel_species_ald_ex2_full_effects.csv"), row.names = FALSE)
write.csv(species_results$clean_out, project_output(file.path(results_dir, "Kiel_species_ALDEx2_FULL_PPI_clean_results.csv"),
    name = "kiel_species_ald_ex2_full_ppi_clean_results.csv"), row.names = FALSE)
write.csv(species_results$clean_out %>% filter(qval < 0.1), project_output(file.path(results_dir, "Kiel_species_ALDEx2_FULL_PPI_q_lt_0.1.csv"),
    name = "kiel_species_ald_ex2_full_ppi_q_lt_0_1.csv"), row.names = FALSE)
saveRDS(species_results, project_output(file.path(results_dir, "Kiel_species_ALDEx2_FULL_results.rds"), name = "kiel_species_ald_ex2_full_results.rds"))
genus_anon <- read.csv(project_input(project_path("data/Community_sim/FOCUS/FOCUS/genus_counts_filt.csv")))
genus_counts <- genus_anon[, -1]
rownames(genus_counts) <- genus_anon$X
colnames(genus_counts) <- stringr::str_sub(colnames(genus_counts), start = 2)
genus_counts <- genus_counts[, focus_data_full$ID, drop = FALSE]
genus_results <- run_aldex_full(counts = genus_counts, covariates = focus_data_full, model_matrix = mm_full,
    feature_level = "genus", mc_samples = 256, ppi_term = "PPI1")
write.csv(genus_results$glm_test_out, project_output(file.path(results_dir, "Kiel_genus_ALDEx2_FULL_all_results.csv"),
    name = "kiel_genus_ald_ex2_full_all_results.csv"), row.names = FALSE)
write.csv(genus_results$effect_out, project_output(file.path(results_dir, "Kiel_genus_ALDEx2_FULL_effects.csv"),
    name = "kiel_genus_ald_ex2_full_effects.csv"), row.names = FALSE)
write.csv(genus_results$clean_out, project_output(file.path(results_dir, "Kiel_genus_ALDEx2_FULL_PPI_clean_results.csv"),
    name = "kiel_genus_ald_ex2_full_ppi_clean_results.csv"), row.names = FALSE)
write.csv(genus_results$clean_out %>% filter(qval < 0.1), project_output(file.path(results_dir, "Kiel_genus_ALDEx2_FULL_PPI_q_lt_0.1.csv"),
    name = "kiel_genus_ald_ex2_full_ppi_q_lt_0_1.csv"), row.names = FALSE)
saveRDS(genus_results, project_output(file.path(results_dir, "Kiel_genus_ALDEx2_FULL_results.rds"), name = "kiel_genus_ald_ex2_full_results.rds"))
message("Species significant q < 0.1: ", sum(species_results$clean_out$qval < 0.1, na.rm = TRUE))
message("Genus significant q < 0.1: ", sum(genus_results$clean_out$qval < 0.1, na.rm = TRUE))
}))
