# SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0
# Metabolomics/RNA sample matching
# Run through run.R so input and output paths resolve within the package.
# Inputs: see data/input_manifest.csv and docs/RUN_ORDER.md.
# Writes generated files under results/; bundled inputs are read-only.

library(dplyr)
library(stringr)
library(writexl)
bridge_dir <- ppi_path("rnaseq/Results/metabolomics_rnaseq_bridges")
dir.create(ppi_output(bridge_dir, directory = TRUE), recursive = TRUE, showWarnings = FALSE)
.norm <- function(x) {
    tolower(gsub("\\s+", "", x))
}
extract_metab_mouse_info <- function(sample_name) {
    snorm <- tolower(gsub("[^[:alnum:]_\\-]+", "_", sample_name))
    m <- stringr::str_match(snorm, "(ctrl|ppi)\\s*[-_]*\\s*(\\d+)")
    grp_raw <- m[, 2]
    donor_num <- m[, 3]
    token_pos <- stringr::str_locate(snorm, "(ctrl|ppi)\\s*[-_]*\\s*(\\d+)")
    after <- ifelse(is.na(token_pos[, 2]), NA_character_, substr(snorm, token_pos[, 2] + 1, nchar(snorm)))
    mouse_num <- stringr::str_match(after, "(\\d+)")[, 2]
    group <- dplyr::case_when(grp_raw == "ctrl" ~ "Ctrl", grp_raw == "ppi" ~ "PPI", TRUE ~ NA_character_)
    tibble::tibble(Sample_Name = sample_name, group_metab = group, donor_num = donor_num, mouse_num = mouse_num)
}
build_metab_rnaseq_bridge <- function(metab_meta, rnaseq_meta, organ_label, sample_col = "Sample_Name") {
    parsed <- extract_metab_mouse_info(metab_meta[[sample_col]])
    parsed <- parsed %>% dplyr::mutate(external_name_key = ifelse(is.na(group_metab) | is.na(donor_num) | 
        is.na(mouse_num), NA_character_, paste0(group_metab, "_", donor_num, "_m_", mouse_num, "_", organ_label)))
    rnaseq_small <- rnaseq_meta %>% dplyr::select(rna_sample = name, external_name, donor, group_rnaseq = group, 
        organ)
    bridge <- parsed %>% dplyr::left_join(rnaseq_small, by = c(external_name_key = "external_name")) %>% 
        dplyr::mutate(matched = !is.na(rna_sample))
    bridge
}
summarise_bridge <- function(bridge_tbl, dataset, compartment, organ_label) {
    bridge_tbl %>% dplyr::summarise(dataset = dataset, compartment = compartment, organ = organ_label, 
        n_metab_samples = dplyr::n(), n_matched = sum(matched), n_unmatched = sum(!matched), matched_fraction = n_matched/n_metab_samples, 
        .groups = "drop")
}
rnaseq_meta_by_org <- list(Brain = meta_brain, Colon = meta_colon, Liver = meta_liver)
gut_meta <- pos_mice$metadata
gut_meta_cecum <- gut_meta %>% dplyr::filter(.norm(Organ) == "cecum")
gut_meta_feces <- gut_meta %>% dplyr::filter(.norm(Organ) %in% c("feces", "faeces", "stool"))
bridges_gut <- list()
for (org in names(rnaseq_meta_by_org)) {
    bridges_gut[[org]] <- list(cecum = build_metab_rnaseq_bridge(metab_meta = gut_meta_cecum, rnaseq_meta = rnaseq_meta_by_org[[org]], 
        organ_label = org), feces = build_metab_rnaseq_bridge(metab_meta = gut_meta_feces, rnaseq_meta = rnaseq_meta_by_org[[org]], 
        organ_label = org))
}
serum_meta <- pos_serum$metadata
bridges_serum <- list()
for (org in names(rnaseq_meta_by_org)) {
    bridges_serum[[org]] <- build_metab_rnaseq_bridge(metab_meta = serum_meta, rnaseq_meta = rnaseq_meta_by_org[[org]], 
        organ_label = org)
}
bridge_diag_gut <- dplyr::bind_rows(lapply(names(bridges_gut), function(org) {
    dplyr::bind_rows(summarise_bridge(bridges_gut[[org]]$cecum, dataset = "gut", compartment = "cecum", 
        organ_label = org), summarise_bridge(bridges_gut[[org]]$feces, dataset = "gut", compartment = "feces", 
        organ_label = org))
}))
bridge_diag_serum <- dplyr::bind_rows(lapply(names(bridges_serum), function(org) {
    summarise_bridge(bridges_serum[[org]], dataset = "serum", compartment = "serum", organ_label = org)
}))
bridge_diagnostics <- dplyr::bind_rows(bridge_diag_gut, bridge_diag_serum)
print(bridge_diagnostics)
unmatched_gut <- dplyr::bind_rows(lapply(names(bridges_gut), function(org) {
    dplyr::bind_rows(bridges_gut[[org]]$cecum %>% dplyr::filter(!matched) %>% dplyr::mutate(dataset = "gut", 
        compartment = "cecum", target_organ = org), bridges_gut[[org]]$feces %>% dplyr::filter(!matched) %>% 
        dplyr::mutate(dataset = "gut", compartment = "feces", target_organ = org))
}))
unmatched_serum <- dplyr::bind_rows(lapply(names(bridges_serum), function(org) {
    bridges_serum[[org]] %>% dplyr::filter(!matched) %>% dplyr::mutate(dataset = "serum", compartment = "serum", 
        target_organ = org)
}))
unmatched_all <- dplyr::bind_rows(unmatched_gut, unmatched_serum)
bridges <- bridges_gut
serum_bridges <- bridges_serum
cat("\nBridge objects created:\n")
cat("bridges        = gut cecum/feces bridges\n")
cat("serum_bridges  = serum bridges\n")
