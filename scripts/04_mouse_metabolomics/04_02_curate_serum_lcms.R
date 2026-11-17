# SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0
# Serum feature curation

source("scripts/helpers/setup.R")
invisible(run_analysis(
  section = "04_mouse_metabolomics",
  directory = "metabolomics",
  script = "04_02_curate_serum_lcms.R",
  restricted = FALSE,
  code = {
required_packages <- c("dplyr", "tidyr", "purrr", "stringr", "tibble", "readr", "writexl")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]
if (length(missing_packages) > 0) {
    stop("Install these packages first: ", paste(missing_packages, collapse = ", "))
}
review_path <- project_path("metabolomics/Serum_metabolomics/20260121_data/20260630_duplicate_metabolite_original_features_id_version01.txt")
out_dir <- project_path("metabolomics/Serum_metabolomics/20260121_data/Serum_cleaned_after_collaborator_review")
dir.create(project_output(out_dir, directory = TRUE), recursive = TRUE, showWarnings = FALSE)
meta_cols <- 7
run_mixed_model <- TRUE
sample_column <- "Sample_Name"
group_column <- "Separation"
group_levels <- c("ctrl", "ppi")
donor_column <- "Donor"
rt_tolerance_minutes <- 0.1
exclude_unresolved_rank2_conflicts <- FALSE
manual_exclude_metabolites <- c("Tributyl phosphate")
pos_data_name <- "2026-01-21_20251205_AW_PPI_KIEL_STUDYII_MS_C8_RP_BA_pos_data.xlsx"
pos_metid_name <- "2026-01-21_20251205_AW_PPI_KIEL_STUDYII_MS_C8_RP_BA_pos_metid.xlsx"
neg_data_name <- "2026-01-21_20251202_AW_PPI_KIEL_STUDYII_MS_C8_RP_BA_neg_data.xlsx"
neg_metid_name <- "2026-01-21_20251202_AW_PPI_KIEL_STUDYII_MS_C8_RP_BA_neg_metid.xlsx"
if (!exists("f")) {
    stop("Object `f` does not exist. Run your original serum-loading script first.")
}
if (!exists("pos_serum") || !exists("neg_serum")) {
    stop("pos_serum and/or neg_serum do not exist. Run the original serum-loading script first.")
}
if (run_mixed_model && !exists("run_mixedmeta_cross_sectional_single_compartment")) {
    stop("run_mixed_model is TRUE, but ", "run_mixedmeta_cross_sectional_single_compartment() is not loaded. ",
        "Source 00_03_fit_serum_mixed_models.R first or set ", "run_mixed_model <- FALSE.")
}
needed_f_objects <- c(pos_data_name, pos_metid_name, neg_data_name, neg_metid_name)
missing_f_objects <- setdiff(needed_f_objects, names(f))
if (length(missing_f_objects) > 0) {
    stop("These expected items are missing from `f`: ", paste(missing_f_objects, collapse = ", "))
}
if (!project_file_exists(review_path)) {
    stop("Collaborator review file not found: ", review_path)
}
first_matching_column <- function(data, candidates = character(), regex = NULL, required = TRUE, label = "column") {
    direct <- candidates[candidates %in% names(data)]
    if (length(direct) > 0) {
        return(direct[[1]])
    }
    if (!is.null(regex)) {
        regex_hit <- names(data)[stringr::str_detect(tolower(names(data)), regex)]
        if (length(regex_hit) > 0) {
            return(regex_hit[[1]])
        }
    }
    if (required) {
        stop("Could not identify ", label, ". Available columns: ", paste(names(data), collapse = ", "))
    }
    NA_character_
}
clean_text_na <- function(x) {
    x <- as.character(x)
    x[is.na(x)] <- ""
    x <- stringr::str_trim(x)
    x[x == ""] <- NA_character_
    x
}
strip_adduct_suffix <- function(x) {
    x <- as.character(x)
    old <- NULL
    while (!identical(old, x)) {
        old <- x
        x <- stringr::str_remove(x, stringr::regex(paste0("(", "_M\\+Cl", "|_M\\+AceticAcid-H", "|_2M-H", "|_M\\+H",
            "|_M\\+NH4", "|_NH4", "|_2M\\+H", ")$"), ignore_case = TRUE))
    }
    stringr::str_trim(x)
}
extract_renamed_annotation <- function(comment) {
    comment <- clean_text_na(comment)
    out <- rep(NA_character_, length(comment))
    pattern_as <- stringr::regex("rename\\s+as\\s+(.+)$", ignore_case = TRUE)
    pattern_underscore <- stringr::regex("rename_([^;]+)$", ignore_case = TRUE)
    has_as <- !is.na(comment) & stringr::str_detect(comment, pattern_as)
    out[has_as] <- stringr::str_match(comment[has_as], pattern_as)[, 2]
    has_under <- is.na(out) & !is.na(comment) & stringr::str_detect(comment, pattern_underscore)
    out[has_under] <- stringr::str_match(comment[has_under], pattern_underscore)[, 2]
    out <- stringr::str_remove(out, stringr::regex("_(isomer|remove.*|exclude.*|keep.*)$", ignore_case = TRUE))
    clean_text_na(out)
}
make_metid_key <- function(metid, mode) {
    id_col <- first_matching_column(metid, candidates = c("id", "ID", "Feature ID", "Feature_ID", "feature_id"),
        regex = "(^|_)id$|feature.*id", label = "feature-ID column")
    rank_col <- first_matching_column(metid, candidates = c("Rank", "rank"), regex = "^rank$", label = "rank column")
    rank1_col <- first_matching_column(metid, candidates = c("Metabolite_Rank1", "metabolite_rank1"),
        regex = "^metabolite.*rank.?1$", required = FALSE, label = "Rank-1 annotation column")
    rank2_cols <- names(metid)[stringr::str_detect(tolower(names(metid)), "^metabolite_rank2($|_)")]
    rank2_cols <- c(intersect(c("Metabolite_Rank2", "metabolite_rank2"), rank2_cols), setdiff(rank2_cols,
        c("Metabolite_Rank2", "metabolite_rank2")))
    general_annotation_col <- first_matching_column(metid, candidates = c("Metabolite", "metabolite",
        "Annotation", "annotation", "name"), regex = "annotation$|metabolite$|compound.*name", required = FALSE,
        label = "general annotation column")
    mz_col <- first_matching_column(metid, candidates = c("m/z", "mz", "MZ"), regex = "^m.?z$|mass.*charge",
        required = FALSE, label = "m/z column")
    rt_col <- first_matching_column(metid, candidates = c("RT", "rt", "RT (min)", "retention_time", "retention time"),
        regex = "^rt|retention", required = FALSE, label = "retention-time column")
    rank_values <- as.character(metid[[rank_col]])
    rank1_annotation <- rep(NA_character_, nrow(metid))
    if (!is.na(rank1_col)) {
        rank1_annotation <- clean_text_na(metid[[rank1_col]])
    }
    if (length(rank2_cols) > 0) {
        rank2_matrix <- as.data.frame(lapply(metid[rank2_cols], clean_text_na), stringsAsFactors = FALSE,
            check.names = FALSE)
        rank2_candidate_list <- lapply(seq_len(nrow(rank2_matrix)), function(i) {
            values <- unlist(rank2_matrix[i, , drop = FALSE], use.names = FALSE)
            values <- values[!is.na(values) & trimws(values) != ""]
            unique(as.character(values))
        })
        rank2_annotation <- vapply(rank2_candidate_list, function(values) {
            if (length(values) == 0) {
                NA_character_
            }
            else {
                values[[1]]
            }
        }, FUN.VALUE = character(1))
        rank2_candidate_count <- lengths(rank2_candidate_list)
        rank2_candidates <- vapply(rank2_candidate_list, function(values) {
            if (length(values) == 0) {
                NA_character_
            }
            else {
                paste(values, collapse = " | ")
            }
        }, FUN.VALUE = character(1))
    }
    else {
        rank2_annotation <- rep(NA_character_, nrow(metid))
        rank2_candidate_count <- integer(nrow(metid))
        rank2_candidates <- rep(NA_character_, nrow(metid))
    }
    annotation_general <- rep(NA_character_, nrow(metid))
    if (!is.na(general_annotation_col)) {
        annotation_general <- clean_text_na(metid[[general_annotation_col]])
    }
    annotation <- dplyr::case_when(rank_values == "Rank1" ~ rank1_annotation, rank_values == "Rank2" ~
        rank2_annotation, TRUE ~ annotation_general)
    annotation <- dplyr::coalesce(annotation, rank1_annotation, rank2_annotation, annotation_general)
    tibble::tibble(feature_id = as.character(metid[[id_col]]), mode = mode, rank = rank_values, original_annotation = annotation,
        rank2_candidates = rank2_candidates, rank2_candidate_count = rank2_candidate_count, rank2_annotation_conflict = rank2_candidate_count >
            1, mz_metid = if (!is.na(mz_col)) {
            suppressWarnings(as.numeric(metid[[mz_col]]))
        }
        else {
            NA_real_
        }, rt_metid = if (!is.na(rt_col)) {
            suppressWarnings(as.numeric(metid[[rt_col]]))
        }
        else {
            NA_real_
        }) %>% dplyr::mutate(feature_key = paste(mode, feature_id, sep = "::")) %>% dplyr::distinct(feature_key,
        .keep_all = TRUE)
}
review <- readr::read_tsv(project_input(review_path), show_col_types = FALSE, progress = FALSE) %>% dplyr::transmute(reviewed_metabolite = as.character(`Metabolite annotation`),
    feature_id = as.character(`Original feature ID`), mode = stringr::str_to_title(as.character(Mode)),
    collaborator_rank = as.character(Rank), collaborator_original_annotation = as.character(`Original assigned annotation`),
    mz_review = suppressWarnings(as.numeric(`m/z`)), rt_review = suppressWarnings(as.numeric(`RT (min)`)),
    collaborator_status = stringr::str_to_lower(stringr::str_trim(as.character(Statistics))), observation = as.character(Observation),
    comment = as.character(Comment), collaborator_p_value = suppressWarnings(as.numeric(p.value)), renamed_annotation = extract_renamed_annotation(Comment),
    feature_key = paste(mode, feature_id, sep = "::")) %>% dplyr::mutate(reviewed_include = collaborator_status ==
    "valid")
status_values <- unique(review$collaborator_status)
unexpected_status <- setdiff(status_values, c("valid", "x"))
if (length(unexpected_status) > 0) {
    stop("Unexpected values in Statistics: ", paste(unexpected_status, collapse = ", "))
}
if (anyDuplicated(review$feature_key)) {
    stop("The collaborator table contains duplicate mode/feature-ID keys.")
}
pos_data <- as.data.frame(f[[pos_data_name]], check.names = FALSE)
neg_data <- as.data.frame(f[[neg_data_name]], check.names = FALSE)
pos_metid <- as.data.frame(f[[pos_metid_name]], check.names = FALSE)
neg_metid <- as.data.frame(f[[neg_metid_name]], check.names = FALSE)
pos_key <- make_metid_key(pos_metid, "Positive")
neg_key <- make_metid_key(neg_metid, "Negative")
feature_key_all <- dplyr::bind_rows(pos_key, neg_key) %>% dplyr::filter(rank %in% c("Rank1", "Rank2"))
rank2_annotation_conflicts <- feature_key_all %>% dplyr::filter(rank == "Rank2", rank2_annotation_conflict) %>%
    dplyr::select(feature_key, feature_id, mode, rank, original_annotation, rank2_candidates, rank2_candidate_count)
if (nrow(rank2_annotation_conflicts) > 0) {
    warning(nrow(rank2_annotation_conflicts), " Rank-2 feature(s) contain more than one distinct annotation ",
        "across Metabolite_Rank2 columns. The first annotation was used; ", "inspect Serum_Rank2_annotation_conflicts.csv.")
    readr::write_csv(rank2_annotation_conflicts, project_output(file.path(out_dir, "Serum_Rank2_annotation_conflicts.csv"),
        name = "serum_rank2_annotation_conflicts.csv"))
}
pos_feature_ids <- intersect(pos_key$feature_id[pos_key$rank %in% c("Rank1", "Rank2")], names(pos_data))
neg_feature_ids <- intersect(neg_key$feature_id[neg_key$rank %in% c("Rank1", "Rank2")], names(neg_data))
if (length(pos_feature_ids) == 0) {
    stop("No positive-mode Rank1/Rank2 feature IDs were found in the raw data.")
}
if (length(neg_feature_ids) == 0) {
    stop("No negative-mode Rank1/Rank2 feature IDs were found in the raw data.")
}
meta_pos <- as.data.frame(pos_serum$metadata)
meta_neg <- as.data.frame(neg_serum$metadata)
if (!"Sample_Name" %in% names(meta_pos)) {
    stop("pos_serum$metadata does not contain Sample_Name.")
}
if (!"Sample_Name" %in% names(meta_neg)) {
    stop("neg_serum$metadata does not contain Sample_Name.")
}
if (anyDuplicated(meta_pos$Sample_Name)) {
    stop("Positive metadata contains duplicate Sample_Name values.")
}
if (anyDuplicated(meta_neg$Sample_Name)) {
    stop("Negative metadata contains duplicate Sample_Name values.")
}
make_sample_key <- function(x) {
    x <- as.character(x)
    x <- stringr::str_trim(x)
    x <- stringr::str_to_lower(x)
    x <- stringr::str_remove(x, "_[0-9]+\\.d$")
    x
}
meta_pos$sample_key <- make_sample_key(meta_pos$Sample_Name)
meta_neg$sample_key <- make_sample_key(meta_neg$Sample_Name)
if (anyDuplicated(meta_pos$sample_key)) {
    duplicated_keys <- unique(meta_pos$sample_key[duplicated(meta_pos$sample_key) | duplicated(meta_pos$sample_key,
        fromLast = TRUE)])
    stop("Positive-mode metadata contains duplicated sample keys after ", "removing acquisition numbers: ",
        paste(duplicated_keys, collapse = ", "))
}
if (anyDuplicated(meta_neg$sample_key)) {
    duplicated_keys <- unique(meta_neg$sample_key[duplicated(meta_neg$sample_key) | duplicated(meta_neg$sample_key,
        fromLast = TRUE)])
    stop("Negative-mode metadata contains duplicated sample keys after ", "removing acquisition numbers: ",
        paste(duplicated_keys, collapse = ", "))
}
common_sample_keys <- intersect(meta_pos$sample_key, meta_neg$sample_key)
if (length(common_sample_keys) == 0) {
    stop("No biological samples are shared between positive and negative ", "mode after removing terminal acquisition numbers.")
}
common_sample_keys <- meta_pos$sample_key[meta_pos$sample_key %in% common_sample_keys]
pos_row_index <- match(common_sample_keys, meta_pos$sample_key)
neg_row_index <- match(common_sample_keys, meta_neg$sample_key)
metadata_aligned <- meta_pos[pos_row_index, , drop = FALSE]
if (!identical(metadata_aligned$sample_key, meta_neg$sample_key[neg_row_index])) {
    stop("Positive and negative biological samples could not be aligned.")
}
common_samples <- metadata_aligned$Sample_Name
X_pos <- as.matrix(pos_data[pos_row_index, pos_feature_ids, drop = FALSE])
X_neg <- as.matrix(neg_data[neg_row_index, neg_feature_ids, drop = FALSE])
storage.mode(X_pos) <- "numeric"
storage.mode(X_neg) <- "numeric"
colnames(X_pos) <- paste("Positive", pos_feature_ids, sep = "::")
colnames(X_neg) <- paste("Negative", neg_feature_ids, sep = "::")
X_all <- cbind(X_pos, X_neg)
rownames(X_all) <- common_samples
feature_decisions <- feature_key_all %>% dplyr::filter(feature_key %in% colnames(X_all)) %>% dplyr::left_join(review,
    by = c("feature_key", "feature_id", "mode")) %>% dplyr::mutate(reviewed = !is.na(collaborator_status),
    include_after_review = dplyr::case_when(reviewed ~ reviewed_include, TRUE ~ TRUE), corrected_annotation_raw = dplyr::coalesce(renamed_annotation,
        reviewed_metabolite, original_annotation), corrected_metabolite = strip_adduct_suffix(corrected_annotation_raw),
    mz = dplyr::coalesce(mz_review, mz_metid), rt_min = dplyr::coalesce(rt_review, rt_metid), manual_contaminant_exclusion = corrected_metabolite %in%
        manual_exclude_metabolites, collaborator_resolved_annotation = !is.na(renamed_annotation) | !is.na(reviewed_metabolite),
    unresolved_rank2_conflict = rank == "Rank2" & rank2_annotation_conflict & !collaborator_resolved_annotation,
    include_final = include_after_review & !manual_contaminant_exclusion & !(exclude_unresolved_rank2_conflicts &
        unresolved_rank2_conflict), exclusion_reason = dplyr::case_when(reviewed & !reviewed_include ~
        "Collaborator marked x", manual_contaminant_exclusion ~ "Manual contaminant exclusion", exclude_unresolved_rank2_conflicts &
        unresolved_rank2_conflict ~ "Unresolved competing Rank-2 annotations", TRUE ~ NA_character_))
missing_review_ids <- setdiff(review$feature_key, feature_decisions$feature_key)
if (length(missing_review_ids) > 0) {
    warning("These reviewed features were not present among the reconstructed ", "Rank1/Rank2 raw features:\n",
        paste(missing_review_ids, collapse = "\n"))
}
unresolved_rank2_conflicts <- feature_decisions %>% dplyr::filter(unresolved_rank2_conflict) %>% dplyr::select(feature_key,
    feature_id, mode, rank, original_annotation, rank2_candidates, rank2_candidate_count, reviewed_metabolite,
    renamed_annotation, collaborator_status, observation, comment, include_final, exclusion_reason)
if (nrow(unresolved_rank2_conflicts) > 0) {
    warning(nrow(unresolved_rank2_conflicts), " feature(s) have competing Rank-2 annotations. They are retained ",
        "under Metabolite_Rank2 for the primary analysis and exported for ", "audit in Serum_competing_Rank2_annotations_audit.csv.")
    readr::write_csv(unresolved_rank2_conflicts, project_output(file.path(out_dir, "Serum_competing_Rank2_annotations_audit.csv"),
        name = "serum_competing_rank2_annotations_audit.csv"))
}
missing_annotation_features <- feature_decisions %>% dplyr::filter(include_final, is.na(corrected_metabolite) |
    corrected_metabolite == "") %>% dplyr::select(feature_key, feature_id, mode, rank, original_annotation,
    reviewed_metabolite, renamed_annotation, collaborator_status, rt_min, mz)
if (nrow(missing_annotation_features) > 0) {
    message("\nFeatures lacking a usable metabolite annotation: ", nrow(missing_annotation_features))
    print(missing_annotation_features, n = min(50, nrow(missing_annotation_features)))
    readr::write_csv(missing_annotation_features, project_output(file.path(out_dir, "Serum_features_missing_annotation_diagnostic.csv"),
        name = "serum_features_missing_annotation_diagnostic.csv"))
    missing_fraction <- nrow(missing_annotation_features)/sum(feature_decisions$include_final)
    if (nrow(missing_annotation_features) > 20 || missing_fraction > 0.05) {
        stop("Too many included features lack an annotation (", nrow(missing_annotation_features), "; ",
            round(100 * missing_fraction, 1), "%). This likely indicates that the annotation columns in the ",
            "metid sheets were not parsed correctly. Inspect ", "Serum_features_missing_annotation_diagnostic.csv and run ",
            "names(pos_metid) / names(neg_metid).")
    }
    warning(nrow(missing_annotation_features), " feature(s) had no assigned metabolite annotation and will be ",
        "excluded, matching the original annotated-feature analysis.")
    feature_decisions <- feature_decisions %>% dplyr::mutate(missing_annotation = is.na(corrected_metabolite) |
        corrected_metabolite == "", include_final = include_final & !missing_annotation, exclusion_reason = dplyr::case_when(missing_annotation ~
        "No assigned metabolite annotation", TRUE ~ exclusion_reason))
} else {
    feature_decisions <- feature_decisions %>% dplyr::mutate(missing_annotation = FALSE)
}
included_features <- feature_decisions %>% dplyr::filter(include_final)
sample_names_lower <- tolower(metadata_aligned$Sample_Name)
biological_index <- stringr::str_detect(sample_names_lower, "ctrl[0-9]+|ppi[0-9]+")
qc_index <- stringr::str_detect(sample_names_lower, "(^|[^a-z])qc([^a-z]|$)")
if (sum(biological_index) == 0) {
    stop("No biological samples were identified from Sample_Name. ", "Edit biological_index if your naming scheme differs.")
}
safe_cv <- function(x) {
    x <- x[is.finite(x) & x > 0]
    if (length(x) < 3) {
        return(NA_real_)
    }
    mean_x <- mean(x)
    if (!is.finite(mean_x) || mean_x == 0) {
        return(NA_real_)
    }
    stats::sd(x)/mean_x
}
quality_metrics <- lapply(colnames(X_all), function(feature_key) {
    x <- X_all[, feature_key]
    x_bio <- x[biological_index]
    x_qc <- x[qc_index]
    tibble::tibble(feature_key = feature_key, biological_detection_n = sum(is.finite(x_bio) & x_bio >
        0), biological_missing_n = sum(!is.finite(x_bio) | x_bio <= 0), biological_median_signal = if (any(is.finite(x_bio) &
        x_bio > 0)) {
        stats::median(x_bio[is.finite(x_bio) & x_bio > 0])
    }
    else {
        NA_real_
    }, qc_n = sum(is.finite(x_qc) & x_qc > 0), qc_cv = safe_cv(x_qc))
}) %>% dplyr::bind_rows()
included_features <- included_features %>% dplyr::left_join(quality_metrics, by = "feature_key") %>%
    dplyr::mutate(rank_priority = dplyr::case_when(rank == "Rank1" ~ 1L, rank == "Rank2" ~ 2L, TRUE ~
        9L), qc_cv_for_ranking = dplyr::if_else(is.na(qc_cv), Inf, qc_cv), median_for_ranking = dplyr::if_else(is.na(biological_median_signal),
        -Inf, biological_median_signal))
assign_cross_mode_clusters <- function(one_metabolite, rt_tolerance) {
    one_metabolite <- one_metabolite %>% dplyr::arrange(rt_min, mode, feature_id) %>% dplyr::mutate(local_cluster = NA_integer_)
    next_cluster <- 1L
    pos_rows <- which(one_metabolite$mode == "Positive" & is.finite(one_metabolite$rt_min))
    neg_rows <- which(one_metabolite$mode == "Negative" & is.finite(one_metabolite$rt_min))
    if (length(pos_rows) > 0 && length(neg_rows) > 0) {
        candidate_pairs <- tidyr::expand_grid(pos_row = pos_rows, neg_row = neg_rows) %>% dplyr::mutate(rt_difference = abs(one_metabolite$rt_min[pos_row] -
            one_metabolite$rt_min[neg_row])) %>% dplyr::filter(rt_difference <= rt_tolerance) %>% dplyr::arrange(rt_difference)
        used_pos <- integer()
        used_neg <- integer()
        if (nrow(candidate_pairs) > 0) {
            for (pair_index in seq_len(nrow(candidate_pairs))) {
                current_pos <- candidate_pairs$pos_row[pair_index]
                current_neg <- candidate_pairs$neg_row[pair_index]
                if (!current_pos %in% used_pos && !current_neg %in% used_neg) {
                  one_metabolite$local_cluster[c(current_pos, current_neg)] <- next_cluster
                  used_pos <- c(used_pos, current_pos)
                  used_neg <- c(used_neg, current_neg)
                  next_cluster <- next_cluster + 1L
                }
            }
        }
    }
    unmatched_rows <- which(is.na(one_metabolite$local_cluster))
    if (length(unmatched_rows) > 0) {
        one_metabolite$local_cluster[unmatched_rows] <- seq.int(from = next_cluster, length.out = length(unmatched_rows))
    }
    one_metabolite
}
clustered_features <- included_features %>% dplyr::group_by(corrected_metabolite) %>% dplyr::group_modify(~assign_cross_mode_clusters(.x,
    rt_tolerance = rt_tolerance_minutes)) %>% dplyr::ungroup() %>% dplyr::mutate(analyte_cluster = paste(corrected_metabolite,
    local_cluster, sep = "::cluster"))
representative_table <- clustered_features %>% dplyr::group_by(analyte_cluster) %>% dplyr::arrange(dplyr::desc(biological_detection_n),
    rank_priority, qc_cv_for_ranking, dplyr::desc(median_for_ranking), feature_key, .by_group = TRUE) %>%
    dplyr::mutate(representative_rank = dplyr::row_number(), selected_representative = representative_rank ==
        1L, number_features_in_cluster = dplyr::n(), number_modes_in_cluster = dplyr::n_distinct(mode),
        cluster_rt = if (any(is.finite(rt_min))) {
            stats::median(rt_min[is.finite(rt_min)])
        } else {
            NA_real_
        }, detection_tie = sum(biological_detection_n == max(biological_detection_n, na.rm = TRUE)) >
            1, selection_used_mean_tiebreaker = selected_representative & number_features_in_cluster >
            1 & detection_tie & !is.finite(min(qc_cv_for_ranking, na.rm = TRUE))) %>% dplyr::ungroup()
cluster_names <- representative_table %>% dplyr::distinct(corrected_metabolite, analyte_cluster, cluster_rt) %>%
    dplyr::group_by(corrected_metabolite) %>% dplyr::arrange(cluster_rt, analyte_cluster, .by_group = TRUE) %>%
    dplyr::mutate(number_analytes_with_same_name = dplyr::n(), analyte_index = dplyr::row_number(), final_analyte_name = dplyr::case_when(number_analytes_with_same_name ==
        1L ~ corrected_metabolite, is.finite(cluster_rt) ~ paste0(corrected_metabolite, " [RT ", formatC(cluster_rt,
        format = "f", digits = 2), "]"), TRUE ~ paste0(corrected_metabolite, " [isomer ", analyte_index,
        "]"))) %>% dplyr::ungroup()
representative_table <- representative_table %>% dplyr::left_join(cluster_names, by = c("corrected_metabolite",
    "analyte_cluster", "cluster_rt"))
cluster_candidates <- representative_table %>% dplyr::filter(selected_representative) %>% dplyr::mutate(detection_rate = biological_detection_n/sum(biological_index),
    qc_available = is.finite(qc_cv), standard_confirmed = FALSE, expected_ion = TRUE, peak_quality_pass = TRUE,
    blank_quality_pass = TRUE)
manual_cluster_priority <- tibble::tribble(~analyte_cluster, ~standard_confirmed, ~expected_ion, ~peak_quality_pass,
    ~blank_quality_pass, ~manual_comment)
cluster_candidates <- cluster_candidates %>% dplyr::left_join(manual_cluster_priority, by = "analyte_cluster",
    suffix = c("", "_manual")) %>% dplyr::mutate(standard_confirmed = dplyr::coalesce(standard_confirmed_manual,
    standard_confirmed), expected_ion = dplyr::coalesce(expected_ion_manual, expected_ion), peak_quality_pass = dplyr::coalesce(peak_quality_pass_manual,
    peak_quality_pass), blank_quality_pass = dplyr::coalesce(blank_quality_pass_manual, blank_quality_pass))
metabolite_representative_audit <- cluster_candidates %>% dplyr::group_by(corrected_metabolite) %>% dplyr::arrange(dplyr::desc(standard_confirmed),
    dplyr::desc(expected_ion), dplyr::desc(peak_quality_pass), dplyr::desc(blank_quality_pass), rank_priority,
    dplyr::desc(biological_detection_n), dplyr::desc(qc_available), qc_cv_for_ranking, feature_key, .by_group = TRUE) %>%
    dplyr::mutate(number_rt_clusters = dplyr::n(), metabolite_candidate_rank = dplyr::row_number(), selected_metabolite_representative = metabolite_candidate_rank ==
        1L) %>% dplyr::ungroup()
final_metabolite_representatives <- metabolite_representative_audit %>% dplyr::filter(selected_metabolite_representative)
selected_features <- final_metabolite_representatives %>% dplyr::arrange(corrected_metabolite)
duplicate_selected_metabolites <- selected_features %>% dplyr::count(corrected_metabolite, name = "n") %>%
    dplyr::filter(n != 1L)
if (nrow(duplicate_selected_metabolites) > 0) {
    print(duplicate_selected_metabolites, n = Inf)
    stop("There is not exactly one selected feature per ", "corrected metabolite name.")
}
if (anyNA(selected_features$corrected_metabolite) || any(selected_features$corrected_metabolite == "")) {
    stop("A selected feature has a missing metabolite name.")
}
if (anyDuplicated(selected_features$feature_key)) {
    stop("The same original feature was selected more than once.")
}
selected_keys <- selected_features$feature_key
serum_clean <- X_all[, selected_keys, drop = FALSE]
colnames(serum_clean) <- selected_features$corrected_metabolite
if (anyDuplicated(colnames(serum_clean))) {
    stop("Duplicate corrected metabolite names remain.")
}
if (ncol(serum_clean) != dplyr::n_distinct(selected_features$corrected_metabolite)) {
    stop("The cleaned matrix is not exactly one feature ", "per corrected metabolite name.")
}
if (!identical(rownames(serum_clean), metadata_aligned$Sample_Name)) {
    stop("Cleaned serum matrix and metadata are misaligned.")
}
serum_clean_bio <- serum_clean[biological_index, , drop = FALSE]
metadata_clean_bio <- metadata_aligned[biological_index, , drop = FALSE]
if (nrow(serum_clean_bio) != nrow(metadata_clean_bio)) {
    stop("Biological metadata and cleaned serum matrix differ in row count.")
}
if (!sample_column %in% names(metadata_clean_bio)) {
    stop("Biological metadata does not contain the sample column: ", sample_column)
}
if (!group_column %in% names(metadata_clean_bio)) {
    stop("Biological metadata does not contain the group column: ", group_column)
}
if (!identical(rownames(serum_clean_bio), as.character(metadata_clean_bio[[sample_column]]))) {
    stop("Biological feature matrix and metadata sample order are misaligned.")
}
if (anyDuplicated(colnames(serum_clean_bio))) {
    stop("Duplicate metabolite names remain in the limma input matrix.")
}
metadata_clean_bio[[donor_column]] <- stringr::str_extract(stringr::str_to_lower(as.character(metadata_clean_bio[[sample_column]])),
    "(?:ctrl|ppi)[0-9]+")
if (anyNA(metadata_clean_bio[[donor_column]])) {
    stop("At least one donor ID could not be extracted from ", sample_column, ". Expected strings such as ctrl1 or ppi2.")
}
metadata_clean_bio[[group_column]] <- factor(stringr::str_to_lower(as.character(metadata_clean_bio[[group_column]])),
    levels = group_levels)
if (anyNA(metadata_clean_bio[[group_column]])) {
    stop("Unexpected or missing group labels. Expected: ", paste(group_levels, collapse = ", "))
}
donor_group_check <- metadata_clean_bio %>% dplyr::transmute(Donor = .data[[donor_column]], Group = .data[[group_column]]) %>%
    dplyr::distinct() %>% dplyr::count(Donor, name = "n_groups") %>% dplyr::filter(n_groups != 1L)
if (nrow(donor_group_check) > 0) {
    stop("At least one donor maps to more than one treatment group.")
}
donor_sample_counts <- metadata_clean_bio %>% dplyr::transmute(Donor = .data[[donor_column]], Group = .data[[group_column]]) %>%
    dplyr::count(Donor, Group, name = "n_recipient_mice") %>% dplyr::arrange(Group, Donor)
donor_group_counts <- donor_sample_counts %>% dplyr::distinct(Donor, Group) %>% dplyr::count(Group, name = "n_donors")
if (nrow(donor_group_counts) != 2L || any(donor_group_counts$n_donors != 3L)) {
    warning("Expected three control donors and three PPI donors. Observed: ", paste(paste0(donor_group_counts$Group,
        "=", donor_group_counts$n_donors), collapse = ", "))
}
limma_feature_matrix <- serum_clean_bio
limma_metadata <- metadata_clean_bio
limma_input_summary <- tibble::tibble(item = c("Recipient mice", "Independent donors", "Control donors",
    "PPI donors", "Metabolites", "Duplicate metabolite column names"), value = c(nrow(limma_feature_matrix),
    dplyr::n_distinct(limma_metadata[[donor_column]]), sum(donor_group_counts$n_donors[donor_group_counts$Group ==
        group_levels[1]]), sum(donor_group_counts$n_donors[donor_group_counts$Group == group_levels[2]]),
    ncol(limma_feature_matrix), anyDuplicated(colnames(limma_feature_matrix))))
if (run_mixed_model) {
    res_serum_clean <- run_mixedmeta_cross_sectional_single_compartment(metadata = metadata_clean_bio,
        features = serum_clean_bio, allow_metabolite_duplicates = FALSE, log_transform = TRUE)
} else {
    res_serum_clean <- list(metabolite_level = tibble::tibble(), feature_level = tibble::tibble())
}
post_rename_duplicates <- representative_table %>% dplyr::group_by(corrected_metabolite) %>% dplyr::summarise(retained_features = dplyr::n(),
    retained_analytes = dplyr::n_distinct(analyte_cluster), positive_features = sum(mode == "Positive"),
    negative_features = sum(mode == "Negative"), selected_feature_ids = paste(feature_id[selected_representative],
        collapse = "; "), final_analyte_names = paste(unique(final_analyte_name), collapse = "; "), .groups = "drop") %>%
    dplyr::filter(retained_features > 1 | retained_analytes > 1) %>% dplyr::arrange(dplyr::desc(retained_features),
    corrected_metabolite)
manual_review_needed <- representative_table %>% dplyr::filter(number_features_in_cluster > 1 & selection_used_mean_tiebreaker) %>%
    dplyr::arrange(corrected_metabolite, analyte_cluster, representative_rank)
multiple_rt_cluster_choices <- metabolite_representative_audit %>% dplyr::filter(number_rt_clusters >
    1) %>% dplyr::select(corrected_metabolite, analyte_cluster, feature_key, feature_id, mode, rank,
    rt_min, mz, biological_detection_n, biological_missing_n, biological_median_signal, qc_n, qc_cv,
    standard_confirmed, expected_ion, peak_quality_pass, blank_quality_pass, manual_comment, metabolite_candidate_rank,
    selected_metabolite_representative) %>% dplyr::arrange(corrected_metabolite, metabolite_candidate_rank)
excluded_features <- feature_decisions %>% dplyr::filter(!include_final) %>% dplyr::arrange(corrected_metabolite,
    mode, rt_min)
cleaned_matrix_export <- data.frame(Sample_Name = rownames(serum_clean_bio), serum_clean_bio, check.names = FALSE)
limma_metadata_export <- metadata_clean_bio %>% dplyr::mutate(dplyr::across(dplyr::where(is.factor),
    base::as.character))
model_metabolite_level <- if ("metabolite_level" %in% names(res_serum_clean)) {
    res_serum_clean$metabolite_level
} else {
    tibble::tibble()
}
model_feature_level <- if ("feature_level" %in% names(res_serum_clean)) {
    res_serum_clean$feature_level
} else {
    tibble::tibble()
}
audit_workbook_path <- file.path(out_dir, "Serum_LCMS_collaborator_review_audit_and_results.xlsx")
writexl::write_xlsx(list(Feature_decisions = feature_decisions, Included_features = included_features,
    Analyte_clusters = representative_table, Metabolite_selection_audit = metabolite_representative_audit,
    Final_one_per_metabolite = final_metabolite_representatives, Multiple_RT_cluster_choices = multiple_rt_cluster_choices,
    Selected_representatives = selected_features, Post_rename_duplicates = post_rename_duplicates, Manual_review_needed = manual_review_needed,
    Competing_Rank2 = unresolved_rank2_conflicts, Excluded_features = excluded_features, Limma_input_summary = limma_input_summary,
    Donor_sample_counts = donor_sample_counts, Limma_metadata = limma_metadata_export, Model_metabolite_level = model_metabolite_level,
    Model_feature_level = model_feature_level), path = project_output(audit_workbook_path))
readr::write_csv(cleaned_matrix_export, project_output(file.path(out_dir, "Serum_LCMS_cleaned_feature_matrix.csv"),
    name = "serum_lcms_cleaned_feature_matrix.csv"))
readr::write_csv(limma_metadata_export, project_output(file.path(out_dir, "Serum_LCMS_limma_metadata.csv"),
    name = "serum_lcms_limma_metadata.csv"))
readr::write_csv(selected_features, project_output(file.path(out_dir, "Serum_LCMS_selected_representatives.csv"),
    name = "serum_lcms_selected_representatives.csv"))
readr::write_csv(multiple_rt_cluster_choices, project_output(file.path(out_dir, "Serum_LCMS_multiple_RT_cluster_choices.csv"),
    name = "serum_lcms_multiple_rt_cluster_choices.csv"))
if (nrow(model_metabolite_level) > 0) {
    readr::write_csv(model_metabolite_level, project_output(file.path(out_dir, "Serum_LCMS_cleaned_mixed_model_results.csv"),
        name = "serum_lcms_cleaned_mixed_model_results.csv"))
}
main_rds_path <- file.path(out_dir, "Serum_LCMS_cleaned_analysis_objects.rds")
saveRDS(list(review = review, feature_decisions = feature_decisions, included_features = included_features,
    unresolved_rank2_conflicts = unresolved_rank2_conflicts, clustered_features = clustered_features,
    representative_table = representative_table, cluster_candidates = cluster_candidates, metabolite_representative_audit = metabolite_representative_audit,
    final_metabolite_representatives = final_metabolite_representatives, multiple_rt_cluster_choices = multiple_rt_cluster_choices,
    selected_features = selected_features, X_all = X_all, serum_clean = serum_clean, serum_clean_bio = serum_clean_bio,
    metadata_clean_bio = metadata_clean_bio, limma_feature_matrix = limma_feature_matrix, limma_metadata = limma_metadata,
    donor_sample_counts = donor_sample_counts, donor_group_counts = donor_group_counts, limma_input_summary = limma_input_summary,
    run_mixed_model = run_mixed_model, res_serum_clean = res_serum_clean), project_output(main_rds_path))
limma_rds_path <- file.path(out_dir, "Serum_LCMS_limma_input_objects.rds")
saveRDS(list(serum_clean_bio = serum_clean_bio, metadata_clean_bio = metadata_clean_bio, selected_features = selected_features,
    final_metabolite_representatives = final_metabolite_representatives, metabolite_representative_audit = metabolite_representative_audit,
    multiple_rt_cluster_choices = multiple_rt_cluster_choices, donor_sample_counts = donor_sample_counts,
    input_notes = c("Rows are recipient mice and columns are exact corrected metabolite names.", "Exactly one original LC-MS feature is retained per exact corrected metabolite name.",
        "Values are not log transformed in this RDS.", "Use feature-specific log2(intensity + half-minimum-positive) transformation in limma.",
        "Donor is stored explicitly in metadata_clean_bio$Donor.", "Separation is a factor with levels ctrl then ppi.")),
    project_output(limma_rds_path))
message("\nSerum cleaning and limma-input preparation completed.\n", "Original Rank1/Rank2 features reconstructed: ",
    ncol(X_all), "\nFeatures excluded after review/manual contaminant rule: ", nrow(excluded_features),
    "\nFeatures retained before cross-mode representative selection: ", nrow(included_features), "\nFinal exact metabolite names retained: ",
    ncol(serum_clean_bio), "\nRecipient mice in limma input: ", nrow(serum_clean_bio), "\nIndependent donors in limma input: ",
    dplyr::n_distinct(metadata_clean_bio[[donor_column]]), "\nMetabolites with more than one RT-cluster candidate: ",
    dplyr::n_distinct(multiple_rt_cluster_choices$corrected_metabolite), "\nCross-mode clusters requiring a representative choice: ",
    sum(representative_table$number_features_in_cluster > 1 & representative_table$selected_representative),
    "\nClusters for which sample mean was used as the final tie-breaker: ", dplyr::n_distinct(manual_review_needed$analyte_cluster),
    "\nMixed model run: ", run_mixed_model, "\n\nAudit workbook:\n", audit_workbook_path, "\n\nMain analysis RDS:\n",
    main_rds_path, "\n\nLimma input RDS:\n", limma_rds_path)
}))
