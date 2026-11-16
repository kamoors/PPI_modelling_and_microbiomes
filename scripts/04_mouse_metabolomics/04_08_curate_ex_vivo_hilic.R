# SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0
# Ex vivo HILIC curation
# Run through run.R so input and output paths resolve within the package.
# Inputs: see data/input_manifest.csv and docs/RUN_ORDER.md.
# Writes generated files under results/; bundled inputs are read-only.

required_packages <- c("readxl", "dplyr", "tidyr", "purrr", "stringr", "tibble", "readr", "writexl")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]
if (length(missing_packages) > 0) {
    stop("Install these packages first: ", paste(missing_packages, collapse = ", "))
}
pos_data_path <- ppi_path("metabolomics/20250714_MassSpecData/20250714_MassSpecData/2025-07-11_20250516_AW_PPI_Supernatants_HILIC_StudyI_pos_data.xlsx")
pos_metid_path <- ppi_path("metabolomics/20250714_MassSpecData/20250714_MassSpecData/2025-07-11_20250516_AW_PPI_Supernatants_HILIC_StudyI_pos_metid.xlsx")
neg_data_path <- ppi_path("metabolomics/20250714_MassSpecData/20250714_MassSpecData/2025-07-11_20250519_AW_PPI_Supernatants_HILIC_StudyI_neg_data.xlsx")
neg_metid_path <- ppi_path("metabolomics/20250714_MassSpecData/20250714_MassSpecData/2025-07-11_20250519_AW_PPI_Supernatants_HILIC_StudyI_neg_metid.xlsx")
out_dir <- ppi_path("metabolomics/20250714_MassSpecData/Ex_vivo_cleaned_one_per_metabolite")
dir.create(ppi_output(out_dir, directory = TRUE), recursive = TRUE, showWarnings = FALSE)
meta_columns <- c("Sample_Name", "Injektionsorder", "Group", "Label", "Day", "Order", "Separation", "Label_Day", 
    "Color")
representative_min_positive_per_group <- 4L
representative_min_donors_per_group <- 2L
manual_feature_name_overrides <- c()
manual_exclude_feature_keys <- c()
clean_text_na <- function(x) {
    x <- as.character(x)
    x[is.na(x)] <- ""
    x <- stringr::str_squish(x)
    x[x == ""] <- NA_character_
    x
}
strip_method_suffix <- function(x) {
    x <- clean_text_na(x)
    x <- stringr::str_remove(x, stringr::regex(";\\s*(LC-|AIF|CE[0-9]*|MS2|CorrDec).*$", ignore_case = TRUE))
    stringr::str_squish(x)
}
canonicalize_metabolite_name <- function(x) {
    x <- strip_method_suffix(x)
    xl <- stringr::str_to_lower(x)
    dplyr::case_when(is.na(x) ~ NA_character_, xl == "isobutyric acid_butyric acid" ~ "Isobutyrate / Butyrate (unresolved)", 
        xl == "hexanoic acid_4-methylvaleric acid" ~ "Hexanoate / 4-Methylvalerate (unresolved)", xl == 
            "isovaleric acid_valeric acid" ~ "Isovalerate / Valerate (unresolved)", xl == "l-leucine_isoleucine" ~ 
            "Leucine / Isoleucine (unresolved)", xl == "glucose_galactose" ~ "Glucose / Galactose (unresolved)", 
        xl == "3-(4-hydroxyphenyl)-propionic acid_dl-3-phenyllactic acid" ~ "3-(4-Hydroxyphenyl)propionate / 3-Phenyllactate (unresolved)", 
        xl %in% c("proline", "l-proline") ~ "Proline", xl %in% c("histidine", "l-histidine", "histidin") ~ 
            "Histidine", xl %in% c("aspartic acid", "l-aspartic acid", "aspartate") ~ "Aspartate", xl == 
            "l-glutamic acid" ~ "Glutamate", xl == "l-glutamine" ~ "Glutamine", xl == "l-alanine" ~ "Alanine", 
        xl == "l-arginine" ~ "Arginine", xl == "l-asparagine" ~ "Asparagine", xl == "l-citrulline" ~ 
            "Citrulline", xl == "l-leucine" ~ "Leucine", xl == "l-lysine" ~ "Lysine", xl %in% c("l-methionine", 
            "methionin") ~ "Methionine", xl %in% c("l-phenylalanine", "phenylalanine") ~ "Phenylalanine", 
        xl == "l-tryptophan" ~ "Tryptophan", xl == "l-tyrosine" ~ "Tyrosine", xl == "l-valine" ~ "Valine", 
        xl == "lactic acid" ~ "Lactate", xl == "succinic acid" ~ "Succinate", xl == "fumaric acid" ~ 
            "Fumarate", xl == "citric acid" ~ "Citrate", xl == "malate" ~ "Malate", xl == "d-glucose" ~ 
            "Glucose", xl == "gluconic acid" ~ "Gluconate", xl == "orotic acid" ~ "Orotate", xl == "meso-tartrate" ~ 
            "meso-Tartrate", xl == "tricarballylic acid" ~ "Tricarballylate", xl %in% c("pantothenic acid", 
            "vitamin b5") ~ "Pantothenate", xl == "sn-glycero-3-phosphocholine" ~ "Glycerophosphocholine", 
        xl == "pyroglutamic acid" ~ "Pyroglutamate", xl == "pipecolate" ~ "Pipecolate", xl == "5-aminopentanoic acid" ~ 
            "5-Aminopentanoate", xl == "4-guanidinobutanoic acid" ~ "4-Guanidinobutanoate", xl == "2-hydroxy-3-methylpentanoic acid" ~ 
            "2-Hydroxy-3-methylpentanoate", xl == "3-hydroxyglutaric acid" ~ "3-Hydroxyglutarate", xl == 
            "adenine" ~ "Adenine", xl == "guanine" ~ "Guanine", xl == "hypoxanthine" ~ "Hypoxanthine", 
        xl == "uracil" ~ "Uracil", xl == "xanthine" ~ "Xanthine", xl == "pyridoxine" ~ "Pyridoxine", 
        xl == "vanillin" ~ "Vanillin", TRUE ~ x)
}
annotation_flag <- function(name, rank) {
    dplyr::case_when(stringr::str_detect(name, "\\(unresolved\\)$") ~ "unresolved_isomer_or_compound_pair", 
        stringr::str_detect(name, stringr::regex("^NCGC|^Tocris|B-Factor", ignore_case = TRUE)) ~ "database_or_nonstandard_label", 
        rank == "Rank3" ~ "rank3_tentative", TRUE ~ "standard")
}
prepare_annotation <- function(metid_path, mode) {
    x <- readxl::read_xlsx(ppi_input(metid_path), .name_repair = "minimal")
    required <- c("id", "Metabolite_Rank1", "Metabolite_Rank2", "Metabolite_Rank3", "Rank")
    missing <- setdiff(required, names(x))
    if (length(missing) > 0) {
        stop(mode, " annotation file is missing: ", paste(missing, collapse = ", "))
    }
    x %>% dplyr::transmute(feature_id = as.character(id), mode = mode, feature_key = paste0(mode, "::", 
        feature_id), rank = as.character(Rank), rank1 = clean_text_na(Metabolite_Rank1), rank2 = clean_text_na(Metabolite_Rank2), 
        rank3 = clean_text_na(Metabolite_Rank3), selected_annotation_raw = dplyr::case_when(rank == "Rank1" ~ 
            rank1, rank == "Rank2" ~ rank2, rank == "Rank3" ~ rank3, TRUE ~ NA_character_), rank_priority = dplyr::case_when(rank == 
            "Rank1" ~ 1L, rank == "Rank2" ~ 2L, rank == "Rank3" ~ 3L, TRUE ~ 99L)) %>% dplyr::mutate(selected_annotation_clean = strip_method_suffix(selected_annotation_raw), 
        canonical_name = canonicalize_metabolite_name(selected_annotation_clean), canonical_key = stringr::str_to_lower(stringr::str_squish(canonical_name)), 
        annotation_flag = annotation_flag(canonical_name, rank))
}
read_data_file <- function(data_path, mode) {
    x <- readxl::read_xlsx(ppi_input(data_path), .name_repair = "minimal")
    missing_meta <- setdiff(meta_columns, names(x))
    if (length(missing_meta) > 0) {
        stop(mode, " data file is missing metadata columns: ", paste(missing_meta, collapse = ", "))
    }
    feature_columns <- names(x)[stringr::str_detect(names(x), "^ID_[0-9]+$")]
    if (length(feature_columns) == 0) {
        stop("No feature-ID columns found in the ", mode, " data file.")
    }
    metadata <- x %>% dplyr::select(dplyr::all_of(meta_columns)) %>% dplyr::mutate(Separation = dplyr::recode(as.character(Separation), 
        ctlr = "ctrl"), Donor = as.character(Label), Day = as.character(Day), SampleKey = paste(Donor, 
        Day, sep = "__"))
    feature_matrix <- x %>% dplyr::select(dplyr::all_of(feature_columns)) %>% dplyr::mutate(dplyr::across(dplyr::everything(), 
        as.numeric)) %>% as.matrix()
    colnames(feature_matrix) <- paste0(mode, "::", feature_columns)
    rownames(feature_matrix) <- metadata$SampleKey
    list(metadata = metadata, feature_matrix = feature_matrix, feature_ids = feature_columns)
}
safe_median <- function(x) {
    x <- x[is.finite(x)]
    if (length(x) == 0) 
        NA_real_
    else stats::median(x)
}
safe_min_positive <- function(x) {
    x <- x[is.finite(x) & x > 0]
    if (length(x) == 0) 
        NA_real_
    else min(x)
}
pos <- read_data_file(pos_data_path, "Positive")
neg <- read_data_file(neg_data_path, "Negative")
if (anyDuplicated(pos$metadata$SampleKey) || anyDuplicated(neg$metadata$SampleKey)) {
    stop("Duplicate donor-day sample keys exist in one of the data files.")
}
if (!setequal(pos$metadata$SampleKey, neg$metadata$SampleKey)) {
    stop("Positive and negative files do not contain the same donor-day samples.")
}
neg_order <- match(pos$metadata$SampleKey, neg$metadata$SampleKey)
neg_metadata_aligned <- neg$metadata[neg_order, , drop = FALSE]
neg_matrix_aligned <- neg$feature_matrix[neg_order, , drop = FALSE]
alignment_columns <- c("Donor", "Day", "Separation", "Label_Day")
for (nm in alignment_columns) {
    if (!identical(as.character(pos$metadata[[nm]]), as.character(neg_metadata_aligned[[nm]]))) {
        stop("Positive/negative sample mismatch in column: ", nm)
    }
}
metadata_all <- pos$metadata %>% dplyr::mutate(Sample_Name_Positive = Sample_Name, Sample_Name_Negative = neg_metadata_aligned$Sample_Name, 
    Sample_Name = SampleKey, is_medium = Separation == "medium", is_day0 = Day == "d0" & Separation %in% 
        c("ctrl", "ppi"), is_analysis_sample = Separation %in% c("ctrl", "ppi") & Day %in% c("d1", "d2", 
        "d3", "d4"))
feature_matrix_all <- cbind(pos$feature_matrix, neg_matrix_aligned)
storage.mode(feature_matrix_all) <- "numeric"
if (any(feature_matrix_all < 0, na.rm = TRUE)) {
    warning("Negative values were present. They are retained in the audit, but the limma ", "script will truncate negative values to zero before log transformation.")
}
annotation_all <- dplyr::bind_rows(prepare_annotation(pos_metid_path, "Positive"), prepare_annotation(neg_metid_path, 
    "Negative"))
matrix_feature_keys <- colnames(feature_matrix_all)
if (!setequal(annotation_all$feature_key, matrix_feature_keys)) {
    missing_annotation <- setdiff(matrix_feature_keys, annotation_all$feature_key)
    missing_data <- setdiff(annotation_all$feature_key, matrix_feature_keys)
    stop("Annotation/data feature mismatch. Missing annotations: ", paste(missing_annotation, collapse = ", "), 
        "; annotation rows without data: ", paste(missing_data, collapse = ", "))
}
if (length(manual_feature_name_overrides) > 0) {
    unknown_override_keys <- setdiff(names(manual_feature_name_overrides), annotation_all$feature_key)
    if (length(unknown_override_keys) > 0) {
        stop("Unknown manual override feature keys: ", paste(unknown_override_keys, collapse = ", "))
    }
    annotation_all <- annotation_all %>% dplyr::mutate(manual_name_override = unname(manual_feature_name_overrides[feature_key]), 
        canonical_name = dplyr::if_else(!is.na(manual_name_override), as.character(manual_name_override), 
            canonical_name), canonical_key = stringr::str_to_lower(stringr::str_squish(canonical_name)), 
        annotation_flag = annotation_flag(canonical_name, rank))
} else {
    annotation_all$manual_name_override <- NA_character_
}
annotation_all <- annotation_all %>% dplyr::mutate(manually_excluded = feature_key %in% manual_exclude_feature_keys)
analysis_index <- which(metadata_all$is_analysis_sample)
all_biological_index <- which(metadata_all$Separation %in% c("ctrl", "ppi") & metadata_all$Day %in% c("d0", 
    "d1", "d2", "d3", "d4"))
medium_index <- which(metadata_all$is_medium)
day0_index <- which(metadata_all$is_day0)
feature_stats <- purrr::map_dfr(matrix_feature_keys, function(feature_key) {
    values <- feature_matrix_all[, feature_key]
    group_summary <- purrr::map_dfr(c("ctrl", "ppi"), function(group_value) {
        idx <- which(metadata_all$is_analysis_sample & metadata_all$Separation == group_value)
        group_donors <- unique(metadata_all$Donor[idx])
        donor_detected <- vapply(group_donors, function(donor_value) {
            donor_idx <- idx[metadata_all$Donor[idx] == donor_value]
            any(values[donor_idx] > 0, na.rm = TRUE)
        }, FUN.VALUE = logical(1))
        tibble::tibble(group = group_value, n_positive = sum(values[idx] > 0, na.rm = TRUE), n_nonmissing = sum(is.finite(values[idx])), 
            n_donors_detected = sum(donor_detected), median_group = safe_median(values[idx]))
    })
    ctrl_row <- group_summary[group_summary$group == "ctrl", ]
    ppi_row <- group_summary[group_summary$group == "ppi", ]
    medium_value <- if (length(medium_index) == 1) {
        values[medium_index]
    }
    else {
        NA_real_
    }
    median_analysis <- safe_median(values[analysis_index])
    tibble::tibble(feature_key = feature_key, n_positive_analysis = sum(values[analysis_index] > 0, na.rm = TRUE), 
        n_positive_ctrl = ctrl_row$n_positive, n_positive_ppi = ppi_row$n_positive, n_donors_detected_ctrl = ctrl_row$n_donors_detected, 
        n_donors_detected_ppi = ppi_row$n_donors_detected, n_donors_detected_total = ctrl_row$n_donors_detected + 
            ppi_row$n_donors_detected, median_ctrl = ctrl_row$median_group, median_ppi = ppi_row$median_group, 
        median_analysis = median_analysis, min_positive_all_biological = safe_min_positive(values[all_biological_index]), 
        median_day0 = safe_median(values[day0_index]), medium_value = medium_value, biological_to_medium_ratio = dplyr::case_when(!is.finite(medium_value) ~ 
            NA_real_, medium_value == 0 & median_analysis > 0 ~ Inf, medium_value == 0 & median_analysis == 
            0 ~ NA_real_, TRUE ~ median_analysis/medium_value), zero_variance_analysis = stats::var(values[analysis_index], 
            na.rm = TRUE) == 0, passes_representative_prevalence = ctrl_row$n_positive >= representative_min_positive_per_group & 
            ppi_row$n_positive >= representative_min_positive_per_group & ctrl_row$n_donors_detected >= 
            representative_min_donors_per_group & ppi_row$n_donors_detected >= representative_min_donors_per_group)
})
annotation_with_stats <- annotation_all %>% dplyr::left_join(feature_stats, by = "feature_key")
calculate_duplicate_concordance <- function(group_table) {
    member_keys <- group_table$feature_key
    if (length(member_keys) == 1) {
        return(group_table %>% dplyr::transmute(feature_key, median_spearman_to_other_features = NA_real_))
    }
    x <- feature_matrix_all[all_biological_index, member_keys, drop = FALSE]
    rho <- suppressWarnings(stats::cor(x, method = "spearman", use = "pairwise.complete.obs"))
    diag(rho) <- NA_real_
    tibble::tibble(feature_key = member_keys, median_spearman_to_other_features = apply(rho, 1, function(z) {
        z <- z[is.finite(z)]
        if (length(z) == 0) 
            NA_real_
        else stats::median(z)
    }))
}
feature_concordance <- annotation_with_stats %>% dplyr::filter(!is.na(canonical_key), !manually_excluded) %>% 
    dplyr::group_by(canonical_key) %>% dplyr::group_modify(~calculate_duplicate_concordance(.x)) %>% 
    dplyr::ungroup()
annotation_with_stats <- annotation_with_stats %>% dplyr::left_join(
    feature_concordance %>% dplyr::select(feature_key, median_spearman_to_other_features), by = "feature_key")
duplicate_groups <- annotation_with_stats %>% dplyr::filter(!is.na(canonical_key), !manually_excluded) %>% 
    dplyr::group_by(canonical_key, canonical_name) %>% dplyr::summarise(n_features = dplyr::n(), feature_keys = paste(feature_key, 
    collapse = "; "), modes = paste(sort(unique(mode)), collapse = "; "), ranks = paste(sort(unique(rank)), 
    collapse = "; "), median_pairwise_spearman = {
    x <- median_spearman_to_other_features[is.finite(median_spearman_to_other_features)]
    if (length(x) == 0) 
        NA_real_
    else stats::median(x)
}, duplicate_concordance_class = dplyr::case_when(n_features == 1 ~ "single_feature", is.na(median_pairwise_spearman) ~ 
    "not_estimable", median_pairwise_spearman >= 0.7 ~ "high", median_pairwise_spearman >= 0.4 ~ "moderate", 
    TRUE ~ "low_or_discordant"), .groups = "drop")
select_representatives <- function(annotation_table, allowed_ranks, set_name) {
    decision <- annotation_table %>% dplyr::filter(rank %in% allowed_ranks, !is.na(canonical_key), !manually_excluded) %>% 
        dplyr::left_join(duplicate_groups %>% dplyr::select(canonical_key, n_features, median_pairwise_spearman, 
            duplicate_concordance_class), by = "canonical_key") %>% dplyr::group_by(canonical_key) %>% 
        dplyr::arrange(dplyr::desc(passes_representative_prevalence), rank_priority, dplyr::desc(n_positive_analysis), 
            dplyr::desc(n_donors_detected_total), dplyr::desc(dplyr::coalesce(median_spearman_to_other_features, 
                -Inf)), feature_key, .by_group = TRUE) %>% dplyr::mutate(representative_rank_within_metabolite = dplyr::row_number(), 
        selected_representative = representative_rank_within_metabolite == 1L, analysis_set = set_name, 
        selection_reason = dplyr::case_when(selected_representative & passes_representative_prevalence ~ 
            paste("Selected without treatment statistics: prevalence, annotation rank,", "detection and duplicate concordance."), 
            selected_representative & !passes_representative_prevalence ~ paste("No feature in this metabolite group necessarily passed the representative", 
                "prevalence rule; selected by annotation rank and detection for audit."), TRUE ~ "Alternative feature for the same canonical metabolite.")) %>% 
        dplyr::ungroup()
    selected <- decision %>% dplyr::filter(selected_representative) %>% dplyr::arrange(canonical_name)
    if (anyDuplicated(selected$canonical_key)) {
        stop("Representative selection failed: duplicated canonical keys remain.")
    }
    list(decision = decision, selected = selected)
}
primary_selection <- select_representatives(annotation_with_stats, allowed_ranks = c("Rank1", "Rank2"), 
    set_name = "primary_rank1_rank2")
exploratory_selection <- select_representatives(annotation_with_stats, allowed_ranks = c("Rank1", "Rank2", 
    "Rank3"), set_name = "exploratory_all_ranks")
build_selected_matrix <- function(selected_table) {
    x <- feature_matrix_all[, selected_table$feature_key, drop = FALSE]
    colnames(x) <- selected_table$canonical_name
    x
}
primary_matrix_all_samples <- build_selected_matrix(primary_selection$selected)
exploratory_matrix_all_samples <- build_selected_matrix(exploratory_selection$selected)
metadata_export <- metadata_all %>% dplyr::select(Sample_Name, Sample_Name_Positive, Sample_Name_Negative, 
    Donor, Day, Separation, Label_Day, is_medium, is_day0, is_analysis_sample)
primary_export <- dplyr::bind_cols(metadata_export, tibble::as_tibble(primary_matrix_all_samples, .name_repair = "minimal"))
exploratory_export <- dplyr::bind_cols(metadata_export, tibble::as_tibble(exploratory_matrix_all_samples, 
    .name_repair = "minimal"))
readr::write_csv(primary_export, ppi_output(file.path(out_dir, "HILIC_exvivo_primary_one_per_metabolite_all_samples.csv"), 
    name = "hilic_exvivo_primary_one_per_metabolite_all_samples.csv"))
readr::write_csv(exploratory_export, ppi_output(file.path(out_dir, "HILIC_exvivo_exploratory_all_ranks_one_per_metabolite_all_samples.csv"), 
    name = "hilic_exvivo_exploratory_all_ranks_one_per_metabolite_all_samples.csv"))
analysis_object <- list(primary_matrix_all_samples = primary_matrix_all_samples, exploratory_matrix_all_samples = exploratory_matrix_all_samples, 
    metadata_all = metadata_export, annotation_all = annotation_with_stats, duplicate_groups = duplicate_groups, 
    primary_decisions = primary_selection$decision, exploratory_decisions = exploratory_selection$decision, 
    primary_selected_map = primary_selection$selected, exploratory_selected_map = exploratory_selection$selected, 
    settings = list(primary_annotation_ranks = c("Rank1", "Rank2"), exploratory_annotation_ranks = c("Rank1", 
        "Rank2", "Rank3"), representative_min_positive_per_group = representative_min_positive_per_group, 
        representative_min_donors_per_group = representative_min_donors_per_group, day0_retained_for_baseline_sensitivity = TRUE, 
        medium_retained_for_audit_only = TRUE))
saveRDS(analysis_object, ppi_output(file.path(out_dir, "HILIC_exvivo_cleaned_analysis_objects.rds"), 
    name = "ex_vivo_hilic_analysis_inputs.rds"))
name_harmonization <- annotation_all %>% dplyr::filter(!is.na(selected_annotation_raw), selected_annotation_raw != 
    canonical_name) %>% dplyr::select(feature_key, mode, rank, selected_annotation_raw, selected_annotation_clean, 
    canonical_name, annotation_flag, manual_name_override) %>% dplyr::arrange(canonical_name, feature_key)
medium_day0_audit <- annotation_with_stats %>% dplyr::select(feature_key, mode, rank, canonical_name, 
    median_analysis, median_day0, medium_value, biological_to_medium_ratio, n_positive_ctrl, n_positive_ppi, 
    n_donors_detected_ctrl, n_donors_detected_ppi) %>% dplyr::arrange(biological_to_medium_ratio)
settings_table <- tibble::tribble(~setting, ~value, "Primary annotation ranks", "Rank1 and Rank2", "Exploratory annotation ranks", 
    "Rank1, Rank2 and Rank3", "Representative minimum positive samples per group", as.character(representative_min_positive_per_group), 
    "Representative minimum detected donors per group", as.character(representative_min_donors_per_group), 
    "Primary analysis days", "d1-d4", "Day d0", "Retained for baseline-adjusted sensitivity analysis", 
    "Medium sample", "Audit only; not subtracted and not included in limma", "Cross-mode merge limitation", 
    "No RT or m/z supplied; one representative is selected and all alternatives are audited", "Treatment statistics used for selection", 
    "No")
writexl::write_xlsx(list(Sample_metadata = metadata_export, All_annotations = annotation_with_stats, 
    Name_harmonization = name_harmonization, Feature_prevalence = feature_stats, Duplicate_groups = duplicate_groups, 
    Primary_decisions = primary_selection$decision, Primary_selected = primary_selection$selected, Exploratory_decisions = exploratory_selection$decision, 
    Exploratory_selected = exploratory_selection$selected, Medium_day0_audit = medium_day0_audit, Settings = settings_table), 
    path = ppi_output(file.path(out_dir, "HILIC_exvivo_annotation_audit.xlsx"), name = "hilic_exvivo_annotation_audit.xlsx"))
message("Completed HILIC cleaning. Primary Rank1/Rank2 metabolites: ", ncol(primary_matrix_all_samples), 
    "; exploratory all-rank metabolites: ", ncol(exploratory_matrix_all_samples), ".")
