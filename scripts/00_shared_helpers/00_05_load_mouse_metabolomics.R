# SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0
# Metabolomics loading for RNA pathway scores

library(readxl)
library(dplyr)
library(stringr)
gut_dir <- project_path("metabolomics/20250714_MassSpecData/20250714_MassSpecData")
serum_dir <- project_path("metabolomics/Serum_metabolomics/20260121_data")
gut_process_script <- project_code("scripts/00_shared_helpers/00_04_import_metabolite_sheet.R")
serum_process_script <- project_code("scripts/00_shared_helpers/00_02_import_serum_metabolite_sheet.R")
name_map_file <- project_path("metabolomics/20250714_MassSpecData/LCMS_HILIC_lmer_metabolite_all_results_2025_07_aw.xlsx")
read_xlsx_folder <- function(path) {
    files <- project_list_files(path, pattern = "\\.xlsx$", full.names = FALSE)
    out <- lapply(file.path(path, files), project_reader(readxl::read_excel))
    names(out) <- files
    out
}
clean_adduct_names <- function(x) {
    sub("(_M\\+Cl|_M\\+AceticAcid-H|_2M-H|_M\\+H|_M\\+NH4|_NH4|_2M\\+H)$", "", x)
}
sample_key <- function(x) {
    x %>% stringr::str_remove("_01_[0-9]+\\.d$") %>% tolower()
}
add_donor_from_sample <- function(meta, sample_col = "Sample_Name") {
    meta %>% mutate(donor = stringr::str_extract(tolower(.data[[sample_col]]), "(ctrl\\d+|ppi\\d+)"))
}
apply_name_map <- function(mat, map_tbl, old_col = "metabolite_name", new_col = "final_name") {
    old_names <- colnames(mat)
    idx <- match(old_names, map_tbl[[old_col]])
    new_names <- old_names
    ok <- !is.na(idx) & !is.na(map_tbl[[new_col]][idx]) & map_tbl[[new_col]][idx] != ""
    new_names[ok] <- map_tbl[[new_col]][idx[ok]]
    colnames(mat) <- new_names
    mat
}
load_gut_metabolomics <- function(gut_dir, process_script, name_map_file) {
    project_source(process_script)
    f <- read_xlsx_folder(gut_dir)
    pos_mice <- process_metabolite_sheet(data = f$`2025-07-11_20250526_AW_PPI_KIEL_STUDYI_MS_C8_RP_BA_pos_data.xlsx`,
        metid = f$`2025-07-11_20250526_AW_PPI_KIEL_STUDYI_MS_C8_RP_BA_pos_metid.xlsx`)
    neg_mice <- process_metabolite_sheet(data = f$`2025-07-11_20250523_AW_PPI_KIEL_STUDYI_MS_C8_RP_BA_neg_data.xlsx`,
        metid = f$`2025-07-11_20250523_AW_PPI_KIEL_STUDYI_MS_C8_RP_BA_neg_metid.xlsx`)
    new_names_pos <- readxl::read_excel(project_input(name_map_file), sheet = "organs_pos")
    new_names_neg <- readxl::read_excel(project_input(name_map_file), sheet = "organs_neg")
    pos_mice$high_conf <- apply_name_map(pos_mice$high_conf, new_names_pos)
    neg_mice$high_conf <- apply_name_map(neg_mice$high_conf, new_names_neg)
    if (!identical(sample_key(pos_mice$metadata$Sample_Name), sample_key(neg_mice$metadata$Sample_Name))) {
        warning("Gut pos/neg sample names are not identical after key cleaning. Check row order.")
    }
    mice_combined <- cbind(pos_mice$high_conf, neg_mice$high_conf)
    colnames(mice_combined) <- make.unique(colnames(mice_combined), sep = "_dup")
    meta_gut <- pos_mice$metadata %>% add_donor_from_sample("Sample_Name")
    list(pos_mice = pos_mice, neg_mice = neg_mice, mice_combined = mice_combined, metadata = meta_gut,
        meta_cecum = meta_gut %>% filter(Organ == "cecum"), feat_cecum = mice_combined[meta_gut$Organ ==
            "cecum", , drop = FALSE], meta_feces = meta_gut %>% filter(Organ == "feces"), feat_feces = mice_combined[meta_gut$Organ ==
            "feces", , drop = FALSE])
}
load_serum_metabolomics <- function(serum_dir, process_script) {
    project_source(process_script)
    f <- read_xlsx_folder(serum_dir)
    pos_serum <- process_metabolite_sheet(data = f$`2026-01-21_20251205_AW_PPI_KIEL_STUDYII_MS_C8_RP_BA_pos_data.xlsx`,
        metid = f$`2026-01-21_20251205_AW_PPI_KIEL_STUDYII_MS_C8_RP_BA_pos_metid.xlsx`, meta_cols = 7)
    neg_serum <- process_metabolite_sheet(data = f$`2026-01-21_20251202_AW_PPI_KIEL_STUDYII_MS_C8_RP_BA_neg_data.xlsx`,
        metid = f$`2026-01-21_20251202_AW_PPI_KIEL_STUDYII_MS_C8_RP_BA_neg_metid.xlsx`, meta_cols = 7)
    serum_combined <- cbind(pos_serum$high_conf, neg_serum$high_conf)
    colnames(serum_combined) <- clean_adduct_names(colnames(serum_combined))
    colnames(serum_combined) <- make.unique(colnames(serum_combined), sep = "_dup")
    meta_serum <- pos_serum$metadata %>% add_donor_from_sample("Sample_Name")
    list(pos_serum = pos_serum, neg_serum = neg_serum, serum_combined = serum_combined, metadata = meta_serum)
}
gut_metab <- load_gut_metabolomics(gut_dir = gut_dir, process_script = gut_process_script, name_map_file = name_map_file)
serum_metab <- load_serum_metabolomics(serum_dir = serum_dir, process_script = serum_process_script)
pos_mice <- gut_metab$pos_mice
neg_mice <- gut_metab$neg_mice
mice_combined <- gut_metab$mice_combined
meta_cecum_comb <- gut_metab$meta_cecum
feat_cecum_comb <- gut_metab$feat_cecum
meta_feces_comb <- gut_metab$meta_feces
feat_feces_comb <- gut_metab$feat_feces
pos_serum <- serum_metab$pos_serum
neg_serum <- serum_metab$neg_serum
serum_combined <- serum_metab$serum_combined
cat("\nGut metabolomics:\n")
cat("Total samples:", nrow(mice_combined), "\n")
cat("Total metabolites:", ncol(mice_combined), "\n")
print(table(gut_metab$metadata$Organ, gut_metab$metadata$Separation, useNA = "ifany"))
cat("\nCecum:\n")
cat("Samples:", nrow(feat_cecum_comb), "\n")
cat("Metabolites:", ncol(feat_cecum_comb), "\n")
cat("\nFeces:\n")
cat("Samples:", nrow(feat_feces_comb), "\n")
cat("Metabolites:", ncol(feat_feces_comb), "\n")
cat("\nSerum metabolomics:\n")
cat("Samples:", nrow(serum_combined), "\n")
cat("Metabolites:", ncol(serum_combined), "\n")
print(table(serum_metab$metadata$Separation, useNA = "ifany"))
cat("\nDuplicate gut metabolite names after make.unique():", anyDuplicated(colnames(mice_combined)),
    "\n")
cat("Duplicate serum metabolite names after make.unique():", anyDuplicated(colnames(serum_combined)),
    "\n")
