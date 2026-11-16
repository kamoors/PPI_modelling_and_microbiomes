# SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0
# Caecum/faeces feature curation
# Run through run.R so input and output paths resolve within the package.
# Inputs: see data/input_manifest.csv and docs/RUN_ORDER.md.
# Writes generated files under results/; bundled inputs are read-only.

required_packages <- c("readxl", "dplyr", "tidyr", "purrr", "stringr", "tibble", "readr", "writexl")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]
if (length(missing_packages) > 0) {
    stop("Install these packages first: ", paste(missing_packages, collapse = ", "))
}
pos_data_path <- ppi_path("metabolomics/20250714_MassSpecData/20250714_MassSpecData/2025-07-11_20250526_AW_PPI_KIEL_STUDYI_MS_C8_RP_BA_pos_data.xlsx")
pos_metid_path <- ppi_path("metabolomics/20250714_MassSpecData/20250714_MassSpecData/2025-07-11_20250526_AW_PPI_KIEL_STUDYI_MS_C8_RP_BA_pos_metid.xlsx")
neg_data_path <- ppi_path("metabolomics/20250714_MassSpecData/20250714_MassSpecData/2025-07-11_20250523_AW_PPI_KIEL_STUDYI_MS_C8_RP_BA_neg_data.xlsx")
neg_metid_path <- ppi_path("metabolomics/20250714_MassSpecData/20250714_MassSpecData/2025-07-11_20250523_AW_PPI_KIEL_STUDYI_MS_C8_RP_BA_neg_metid.xlsx")
review_path <- ppi_path("metabolomics/20250714_MassSpecData/20250714_MassSpecData/20260702_duplicate_metabolite_original_features_id_version02.txt")
out_dir <- ppi_path("metabolomics/20250714_MassSpecData/Cecum_Feces_cleaned_after_collaborator_review")
dir.create(ppi_output(out_dir, directory = TRUE), recursive = TRUE, showWarnings = FALSE)
meta_cols <- 9L
rt_tolerance_minutes <- 0.1
exclude_unresolved_rank2_conflicts <- FALSE
manual_exclude_metabolites <- character(0)
manual_review_overrides <- tibble::tribble(~feature_key, ~override_status, ~override_reason, "Negative::ID_1051", 
    "valid", "PROVISIONAL: collaborator status was blank; paired positive ID_1786 was x in a merged-consensus pair. Confirm with collaborator.")
manual_analyte_priority <- tibble::tribble(~analyte_cluster, ~priority, ~manual_comment)
manual_feature_cluster_overrides <- tibble::tribble(~feature_key, ~manual_analyte_cluster, ~manual_cluster_reason, 
    "Positive::ID_1598", "Cholic acid::manual_cluster1", "Cholic-acid monomer/adduct pair; retain one representative", 
    "Positive::ID_1932", "Cholic acid::manual_cluster1", "Cholic-acid monomer/adduct pair; retain one representative")
input_paths <- c(pos_data_path, pos_metid_path, neg_data_path, neg_metid_path, review_path)
missing_paths <- input_paths[!ppi_file_exists(input_paths)]
if (length(missing_paths) > 0) {
    stop("These input files were not found:\n", paste(missing_paths, collapse = "\n"))
}
pos_data <- readxl::read_xlsx(ppi_input(pos_data_path), .name_repair = "minimal") %>% as.data.frame(check.names = FALSE)
neg_data <- readxl::read_xlsx(ppi_input(neg_data_path), .name_repair = "minimal") %>% as.data.frame(check.names = FALSE)
pos_metid <- readxl::read_xlsx(ppi_input(pos_metid_path), .name_repair = "minimal") %>% as.data.frame(check.names = FALSE)
neg_metid <- readxl::read_xlsx(ppi_input(neg_metid_path), .name_repair = "minimal") %>% as.data.frame(check.names = FALSE)
if (ncol(pos_data) <= meta_cols || ncol(neg_data) <= meta_cols) {
    stop("meta_cols is too large or the data files contain no feature columns.")
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
        x <- stringr::str_remove(x, stringr::regex(paste0("(", "_M\\+Cl", "|_M\\+AceticAcid-H", "|_2M-H", 
            "|_M\\+H", "|_M\\+NH4", "|_NH4", "|_2M\\+H", ")$"), ignore_case = TRUE))
    }
    stringr::str_trim(x)
}
clean_metabolite_name <- function(x) {
    x <- strip_adduct_suffix(x)
    x <- stringr::str_remove(x, "\\s*\\([0-9]*\\.?[0-9]+\\)\\s*$")
    stringr::str_squish(x)
}
strip_nonchemical_suffixes <- function(x) {
    x <- as.character(x)
    x <- stringr::str_remove(x, stringr::regex("\\s*;\\s*(?:LC-|AIF\\b|CE[0-9]*\\b|MS2(?:Dec)?\\b).*$", 
        ignore_case = TRUE))
    x <- stringr::str_remove(x, stringr::regex("\\s*;\\s*\\[[^\\]]+\\][+-]?\\s*C[0-9].*$", ignore_case = TRUE))
    stringr::str_squish(x)
}
canonicalize_metabolite_name <- function(x) {
    x <- clean_metabolite_name(x)
    x <- strip_nonchemical_suffixes(x)
    x <- stringr::str_squish(x)
    x <- stringr::str_replace(x, stringr::regex("^PC\\(0:0/([^\\)]+)\\)$", ignore_case = TRUE), "LPC \\1")
    x <- stringr::str_replace(x, stringr::regex("^PC\\(([^\\)]+)/0:0\\)$", ignore_case = TRUE), "LPC \\1")
    x <- stringr::str_replace(x, stringr::regex("^PE\\(0:0/([^\\)]+)\\)$", ignore_case = TRUE), "LPE \\1")
    x <- stringr::str_replace(x, stringr::regex("^PE\\(([^\\)]+)/0:0\\)$", ignore_case = TRUE), "LPE \\1")
    exact_map <- c(`cholic acid` = "Cholic acid", `1-stearoyl-sn-glycero-3-phosphocholine` = "LPC 18:0", 
        `sm(d18:1/15:0)` = "SM 33:1;2O", `22-hydoxy-2-hopen-1-one` = "22-Hydroxy-2-hopen-1-one", `3-(3-hydroxyphenyl)propionic acid` = "3-(3-Hydroxyphenyl)propionic acid", 
        `acetyl-carnitine` = "Acetylcarnitine", `azelaic acid` = "Azelaic acid", glucose = "Glucose", 
        proline = "Proline", succinate = "Succinate", sphinganine = "Sphinganine", `oleoyl ethanolamide` = "Oleoylethanolamide")
    lower_x <- stringr::str_to_lower(x)
    matched <- match(lower_x, names(exact_map))
    replace_index <- !is.na(matched)
    x[replace_index] <- unname(exact_map[matched[replace_index]])
    x <- stringr::str_replace(x, stringr::regex("^lpc\\s*", ignore_case = TRUE), "LPC ")
    x <- stringr::str_replace(x, stringr::regex("^lpe\\s*", ignore_case = TRUE), "LPE ")
    x <- stringr::str_replace(x, stringr::regex("^sm\\s*", ignore_case = TRUE), "SM ")
    stringr::str_squish(x)
}
make_canonical_metabolite_key <- function(x) {
    x <- canonicalize_metabolite_name(x)
    x <- stringr::str_to_lower(stringr::str_squish(x))
    clean_text_na(x)
}
canonical_rank2_candidate_count <- function(candidate_string) {
    if (is.na(candidate_string) || trimws(candidate_string) == "") {
        return(0L)
    }
    candidates <- stringr::str_split(candidate_string, stringr::fixed(" | "))[[1]]
    canonical <- make_canonical_metabolite_key(candidates)
    canonical <- canonical[!is.na(canonical)]
    length(unique(canonical))
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
make_sample_key <- function(x) {
    x <- stringr::str_to_lower(stringr::str_trim(as.character(x)))
    stringr::str_remove(x, "_[0-9]+\\.d$")
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
        rank2_annotation <- vapply(rank2_candidate_list, function(values) if (length(values) == 0) 
            NA_character_
        else values[[1]], FUN.VALUE = character(1))
        rank2_candidate_count <- lengths(rank2_candidate_list)
        rank2_candidates <- vapply(rank2_candidate_list, function(values) {
            if (length(values) == 0) 
                NA_character_
            else paste(values, collapse = " | ")
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
            1, mz_metid = if (!is.na(mz_col)) 
            suppressWarnings(as.numeric(metid[[mz_col]]))
        else NA_real_, rt_metid = if (!is.na(rt_col)) 
            suppressWarnings(as.numeric(metid[[rt_col]]))
        else NA_real_) %>% dplyr::mutate(feature_key = paste(mode, feature_id, sep = "::")) %>% dplyr::distinct(feature_key, 
        .keep_all = TRUE)
}
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
matrix_export <- function(x) {
    data.frame(Sample_Name = rownames(x), x, check.names = FALSE)
}
review_raw <- readr::read_tsv(ppi_input(review_path), show_col_types = FALSE, progress = FALSE)
required_review_columns <- c("Metabolite annotation", "Original feature ID", "Mode", "Rank", "Original assigned annotation", 
    "m/z", "RT (min)", "Consolidation outcome", "Statistics", "Observation", "Comment")
missing_review_columns <- setdiff(required_review_columns, names(review_raw))
if (length(missing_review_columns) > 0) {
    stop("Collaborator review is missing: ", paste(missing_review_columns, collapse = ", "))
}
review <- review_raw %>% dplyr::transmute(reviewed_metabolite = clean_text_na(`Metabolite annotation`), 
    feature_id = as.character(`Original feature ID`), mode = stringr::str_to_title(as.character(Mode)), 
    collaborator_rank = as.character(Rank), collaborator_original_annotation = as.character(`Original assigned annotation`), 
    mz_review = suppressWarnings(as.numeric(`m/z`)), rt_review = suppressWarnings(as.numeric(`RT (min)`)), 
    consolidation_outcome = as.character(`Consolidation outcome`), collaborator_status_original = stringr::str_to_lower(stringr::str_trim(as.character(Statistics))), 
    observation = as.character(Observation), comment = as.character(Comment), review_statistic_unused = if ("p_adj_meta" %in% 
        names(review_raw)) {
        suppressWarnings(as.numeric(p_adj_meta))
    } else {
        NA_real_
    }, rename_from_comment = extract_renamed_annotation(Comment), rename_from_observation = extract_renamed_annotation(Observation), 
    renamed_annotation = dplyr::coalesce(rename_from_comment, rename_from_observation), rename_source = dplyr::case_when(!is.na(rename_from_comment) ~ 
        "Comment", !is.na(rename_from_observation) ~ "Observation", TRUE ~ NA_character_), feature_key = paste(mode, 
        feature_id, sep = "::")) %>% dplyr::mutate(collaborator_status_original = dplyr::na_if(collaborator_status_original, 
    "")) %>% dplyr::left_join(manual_review_overrides, by = "feature_key") %>% dplyr::mutate(collaborator_status = dplyr::coalesce(stringr::str_to_lower(override_status), 
    collaborator_status_original), status_was_overridden = !is.na(override_status), reviewed_include = collaborator_status == 
    "valid")
if (anyDuplicated(review$feature_key)) {
    stop("The collaborator table contains duplicate mode/feature-ID keys.")
}
unresolved_review_status <- review %>% dplyr::filter(is.na(collaborator_status) | !collaborator_status %in% 
    c("valid", "x"))
if (nrow(unresolved_review_status) > 0) {
    print(unresolved_review_status, n = Inf)
    stop("At least one collaborator Statistics value is unresolved. Add an ", "explicit entry to manual_review_overrides before continuing.")
}
if (any(review$status_was_overridden)) {
    warning("Manual collaborator-status override(s) are active: ", paste(review$feature_key[review$status_was_overridden], 
        collapse = ", "))
}
pos_key <- make_metid_key(pos_metid, "Positive")
neg_key <- make_metid_key(neg_metid, "Negative")
feature_key_all <- dplyr::bind_rows(pos_key, neg_key) %>% dplyr::filter(rank %in% c("Rank1", "Rank2"))
rank2_annotation_conflicts <- feature_key_all %>% dplyr::filter(rank == "Rank2", rank2_annotation_conflict) %>% 
    dplyr::select(feature_key, feature_id, mode, rank, original_annotation, rank2_candidates, rank2_candidate_count)
pos_feature_ids <- intersect(pos_key$feature_id[pos_key$rank %in% c("Rank1", "Rank2")], names(pos_data))
neg_feature_ids <- intersect(neg_key$feature_id[neg_key$rank %in% c("Rank1", "Rank2")], names(neg_data))
if (length(pos_feature_ids) == 0 || length(neg_feature_ids) == 0) {
    stop("No Rank1/Rank2 feature IDs could be reconstructed in one or both modes.")
}
meta_pos <- pos_data[, seq_len(meta_cols), drop = FALSE]
meta_neg <- neg_data[, seq_len(meta_cols), drop = FALSE]
if (!"Sample_Name" %in% names(meta_pos) || !"Sample_Name" %in% names(meta_neg)) {
    stop("The first metadata columns must contain Sample_Name.")
}
meta_pos$sample_key <- make_sample_key(meta_pos$Sample_Name)
meta_neg$sample_key <- make_sample_key(meta_neg$Sample_Name)
if (anyDuplicated(meta_pos$sample_key) || anyDuplicated(meta_neg$sample_key)) {
    stop("Duplicated biological sample keys were produced within an ion mode.")
}
common_sample_keys <- meta_pos$sample_key[meta_pos$sample_key %in% meta_neg$sample_key]
if (length(common_sample_keys) == 0) {
    stop("No shared positive/negative biological samples were found.")
}
pos_row_index <- match(common_sample_keys, meta_pos$sample_key)
neg_row_index <- match(common_sample_keys, meta_neg$sample_key)
metadata_aligned <- meta_pos[pos_row_index, , drop = FALSE]
if (!identical(metadata_aligned$sample_key, meta_neg$sample_key[neg_row_index])) {
    stop("Positive and negative biological samples could not be aligned.")
}
metadata_fields_to_check <- intersect(c("Order", "Separation", "Organ", "Weight"), intersect(names(meta_pos), 
    names(meta_neg)))
for (field in metadata_fields_to_check) {
    pos_value <- as.character(meta_pos[[field]][pos_row_index])
    neg_value <- as.character(meta_neg[[field]][neg_row_index])
    if (!identical(pos_value, neg_value)) {
        stop("Positive/negative metadata disagree for field: ", field)
    }
}
common_samples <- as.character(metadata_aligned$Sample_Name)
X_pos <- as.matrix(pos_data[pos_row_index, pos_feature_ids, drop = FALSE])
X_neg <- as.matrix(neg_data[neg_row_index, neg_feature_ids, drop = FALSE])
storage.mode(X_pos) <- "numeric"
storage.mode(X_neg) <- "numeric"
colnames(X_pos) <- paste("Positive", pos_feature_ids, sep = "::")
colnames(X_neg) <- paste("Negative", neg_feature_ids, sep = "::")
X_all <- cbind(X_pos, X_neg)
rownames(X_all) <- common_samples
if (any(!is.finite(X_all) & !is.na(X_all))) {
    stop("The combined matrix contains non-finite values other than NA.")
}
if (any(X_all < 0, na.rm = TRUE)) {
    stop("Negative LC-MS intensities were found. This script does not silently ", "truncate LC-MS values; inspect the input normalization first.")
}
metadata_aligned$Separation <- stringr::str_to_lower(as.character(metadata_aligned$Separation))
metadata_aligned$Separation <- dplyr::recode(metadata_aligned$Separation, ctlr = "ctrl", control = "ctrl", 
    ppi = "ppi", .default = metadata_aligned$Separation)
metadata_aligned$Separation <- factor(metadata_aligned$Separation, levels = c("ctrl", "ppi"))
metadata_aligned$Organ <- stringr::str_to_lower(as.character(metadata_aligned$Organ))
metadata_aligned$Donor <- stringr::str_extract(metadata_aligned$sample_key, "(?:ctrl|ppi)[0-9]+")
mouse_match <- stringr::str_match(metadata_aligned$sample_key, "_(?:ctrl|ppi)[0-9]+_([0-9]+)_")
metadata_aligned$Mouse_Number <- mouse_match[, 2]
metadata_aligned$Mouse_ID <- paste0("mouse", metadata_aligned$Mouse_Number)
metadata_aligned$Donor_Mouse_ID <- paste(metadata_aligned$Donor, metadata_aligned$Mouse_ID, sep = "_")
if (anyNA(metadata_aligned$Separation) || anyNA(metadata_aligned$Donor) || anyNA(metadata_aligned$Mouse_Number)) {
    stop("Group, donor, or mouse ID could not be parsed for at least one sample.")
}
if (!all(metadata_aligned$Organ %in% c("cecum", "feces"))) {
    stop("Unexpected organ labels: ", paste(unique(metadata_aligned$Organ), collapse = ", "))
}
paired_mouse_check <- metadata_aligned %>% dplyr::count(Donor_Mouse_ID, Organ, name = "n") %>% tidyr::pivot_wider(names_from = Organ, 
    values_from = n, values_fill = 0)
if (!all(c("cecum", "feces") %in% names(paired_mouse_check)) || any(paired_mouse_check$cecum != 1L) || 
    any(paired_mouse_check$feces != 1L)) {
    warning("Not every parsed mouse has exactly one caecum and one faeces sample.")
}
feature_decisions <- feature_key_all %>% dplyr::filter(feature_key %in% colnames(X_all)) %>% dplyr::left_join(review, 
    by = c("feature_key", "feature_id", "mode")) %>% dplyr::mutate(reviewed = !is.na(collaborator_status), 
    include_after_review = dplyr::case_when(reviewed ~ reviewed_include, TRUE ~ TRUE), corrected_annotation_raw = dplyr::coalesce(renamed_annotation, 
        reviewed_metabolite, original_annotation), corrected_metabolite = canonicalize_metabolite_name(corrected_annotation_raw), 
    canonical_metabolite_key = make_canonical_metabolite_key(corrected_metabolite), mz = dplyr::coalesce(mz_review, 
        mz_metid), rt_min = dplyr::coalesce(rt_review, rt_metid), manual_contaminant_exclusion = canonical_metabolite_key %in% 
        make_canonical_metabolite_key(manual_exclude_metabolites), collaborator_resolved_annotation = !is.na(renamed_annotation) | 
        !is.na(reviewed_metabolite), rank2_canonical_candidate_count = vapply(rank2_candidates, canonical_rank2_candidate_count, 
        FUN.VALUE = integer(1)), rank2_conflict_resolved_by_canonicalization = rank == "Rank2" & rank2_annotation_conflict & 
        rank2_canonical_candidate_count == 1L, unresolved_rank2_conflict = rank == "Rank2" & rank2_annotation_conflict & 
        !collaborator_resolved_annotation & !rank2_conflict_resolved_by_canonicalization, missing_annotation = is.na(corrected_metabolite) | 
        corrected_metabolite == "", include_final = include_after_review & !manual_contaminant_exclusion & 
        !missing_annotation & !(exclude_unresolved_rank2_conflicts & unresolved_rank2_conflict), exclusion_reason = dplyr::case_when(reviewed & 
        !reviewed_include ~ "Collaborator marked x", manual_contaminant_exclusion ~ "Manual contaminant exclusion", 
        missing_annotation ~ "No assigned metabolite annotation", exclude_unresolved_rank2_conflicts & 
            unresolved_rank2_conflict ~ "Unresolved competing Rank-2 annotations", TRUE ~ NA_character_))
missing_review_ids <- setdiff(review$feature_key, feature_decisions$feature_key)
if (length(missing_review_ids) > 0) {
    warning("Reviewed features absent from reconstructed Rank1/Rank2 features: ", paste(missing_review_ids, 
        collapse = ", "))
}
included_features <- feature_decisions %>% dplyr::filter(include_final)
if (nrow(included_features) == 0) {
    stop("No features remain after collaborator curation.")
}
organ_vector <- metadata_aligned$Organ
cecum_index <- organ_vector == "cecum"
feces_index <- organ_vector == "feces"
quality_metrics <- lapply(colnames(X_all), function(feature_key) {
    x <- X_all[, feature_key]
    positive_all <- is.finite(x) & x > 0
    positive_cecum <- positive_all & cecum_index
    positive_feces <- positive_all & feces_index
    tibble::tibble(feature_key = feature_key, total_detection_n = sum(positive_all), cecum_detection_n = sum(positive_cecum), 
        feces_detection_n = sum(positive_feces), minimum_organ_detection_n = min(sum(positive_cecum), 
            sum(positive_feces)), total_missing_or_zero_n = sum(!positive_all), biological_median_signal = if (any(positive_all)) {
            stats::median(x[positive_all])
        }
        else {
            NA_real_
        })
}) %>% dplyr::bind_rows()
included_features <- included_features %>% dplyr::left_join(quality_metrics, by = "feature_key") %>% 
    dplyr::mutate(rank_priority = dplyr::case_when(rank == "Rank1" ~ 1L, rank == "Rank2" ~ 2L, TRUE ~ 
        9L), median_for_audit = dplyr::if_else(is.na(biological_median_signal), -Inf, biological_median_signal))
clustered_features <- included_features %>% dplyr::group_by(canonical_metabolite_key) %>% dplyr::group_modify(~assign_cross_mode_clusters(.x, 
    rt_tolerance = rt_tolerance_minutes)) %>% dplyr::ungroup() %>% dplyr::left_join(manual_feature_cluster_overrides, 
    by = "feature_key") %>% dplyr::mutate(automatic_analyte_cluster = paste(corrected_metabolite, local_cluster, 
    sep = "::cluster"), analyte_cluster = dplyr::coalesce(manual_analyte_cluster, automatic_analyte_cluster), 
    manual_cluster_applied = !is.na(manual_analyte_cluster))
representative_table <- clustered_features %>% dplyr::group_by(analyte_cluster) %>% dplyr::arrange(dplyr::desc(minimum_organ_detection_n), 
    dplyr::desc(total_detection_n), rank_priority, feature_key, .by_group = TRUE) %>% dplyr::mutate(representative_rank = dplyr::row_number(), 
    selected_analyte_representative = representative_rank == 1L, number_features_in_cluster = dplyr::n(), 
    number_modes_in_cluster = dplyr::n_distinct(mode), cluster_rt = if (any(is.finite(rt_min))) {
        stats::median(rt_min[is.finite(rt_min)])
    } else {
        NA_real_
    }) %>% dplyr::ungroup()
cluster_names <- representative_table %>% dplyr::distinct(canonical_metabolite_key, corrected_metabolite, 
    analyte_cluster, cluster_rt) %>% dplyr::group_by(canonical_metabolite_key) %>% dplyr::arrange(cluster_rt, 
    analyte_cluster, .by_group = TRUE) %>% dplyr::mutate(number_analytes_with_same_name = dplyr::n(), 
    analyte_index = dplyr::row_number(), final_analyte_name = dplyr::case_when(number_analytes_with_same_name == 
        1L ~ corrected_metabolite, is.finite(cluster_rt) ~ paste0(corrected_metabolite, " [RT ", formatC(cluster_rt, 
        format = "f", digits = 2), "]"), TRUE ~ paste0(corrected_metabolite, " [isomer ", analyte_index, 
        "]"))) %>% dplyr::ungroup()
representative_table <- representative_table %>% dplyr::left_join(cluster_names, by = c("canonical_metabolite_key", 
    "corrected_metabolite", "analyte_cluster", "cluster_rt"))
selected_analytes <- representative_table %>% dplyr::filter(selected_analyte_representative) %>% dplyr::arrange(final_analyte_name)
if (anyDuplicated(selected_analytes$final_analyte_name)) {
    stop("RT-resolved analyte names are not unique.")
}
rt_resolved_all <- X_all[, selected_analytes$feature_key, drop = FALSE]
colnames(rt_resolved_all) <- selected_analytes$final_analyte_name
metabolite_representative_audit <- selected_analytes %>% dplyr::left_join(manual_analyte_priority, by = "analyte_cluster") %>% 
    dplyr::mutate(priority = dplyr::coalesce(as.integer(priority), 999L)) %>% dplyr::group_by(canonical_metabolite_key) %>% 
    dplyr::arrange(priority, rank_priority, dplyr::desc(minimum_organ_detection_n), dplyr::desc(total_detection_n), 
        feature_key, .by_group = TRUE) %>% dplyr::mutate(number_rt_clusters = dplyr::n(), metabolite_candidate_rank = dplyr::row_number(), 
    selected_metabolite_representative = metabolite_candidate_rank == 1L) %>% dplyr::ungroup()
selected_metabolites <- metabolite_representative_audit %>% dplyr::filter(selected_metabolite_representative) %>% 
    dplyr::arrange(corrected_metabolite)
if (anyDuplicated(selected_metabolites$canonical_metabolite_key)) {
    stop("There is not exactly one selected feature per canonical metabolite.")
}
one_per_metabolite_all <- X_all[, selected_metabolites$feature_key, drop = FALSE]
colnames(one_per_metabolite_all) <- selected_metabolites$corrected_metabolite
split_matrix_and_metadata <- function(matrix_all, organ_name) {
    keep <- metadata_aligned$Organ == organ_name
    matrix_organ <- matrix_all[keep, , drop = FALSE]
    metadata_organ <- metadata_aligned[keep, , drop = FALSE]
    if (!identical(rownames(matrix_organ), as.character(metadata_organ$Sample_Name))) {
        stop("Matrix/metadata mismatch for organ: ", organ_name)
    }
    list(matrix = matrix_organ, metadata = metadata_organ)
}
rt_cecum <- split_matrix_and_metadata(rt_resolved_all, "cecum")
rt_feces <- split_matrix_and_metadata(rt_resolved_all, "feces")
met_cecum <- split_matrix_and_metadata(one_per_metabolite_all, "cecum")
met_feces <- split_matrix_and_metadata(one_per_metabolite_all, "feces")
if (!setequal(rt_cecum$metadata$Donor_Mouse_ID, rt_feces$metadata$Donor_Mouse_ID)) {
    warning("Caecum and faeces do not contain exactly the same parsed mice.")
}
post_rename_duplicates <- representative_table %>% dplyr::group_by(canonical_metabolite_key, corrected_metabolite) %>% 
    dplyr::summarise(retained_original_features = dplyr::n(), retained_analytes = dplyr::n_distinct(analyte_cluster), 
        positive_features = sum(mode == "Positive"), negative_features = sum(mode == "Negative"), selected_feature_keys = paste(feature_key[selected_analyte_representative], 
            collapse = "; "), final_analyte_names = paste(unique(final_analyte_name), collapse = "; "), 
        .groups = "drop") %>% dplyr::filter(retained_original_features > 1 | retained_analytes > 1) %>% 
    dplyr::arrange(dplyr::desc(retained_original_features), corrected_metabolite)
multiple_rt_cluster_choices <- metabolite_representative_audit %>% dplyr::filter(number_rt_clusters > 
    1) %>% dplyr::arrange(corrected_metabolite, metabolite_candidate_rank)
excluded_features <- feature_decisions %>% dplyr::filter(!include_final) %>% dplyr::arrange(corrected_metabolite, 
    mode, rt_min)
metadata_export <- metadata_aligned %>% dplyr::mutate(dplyr::across(dplyr::where(is.factor), as.character))
name_harmonization_audit <- feature_decisions %>% dplyr::filter(include_final, !is.na(corrected_annotation_raw), 
    corrected_annotation_raw != corrected_metabolite | !is.na(renamed_annotation) | rank2_conflict_resolved_by_canonicalization) %>% 
    dplyr::select(feature_key, mode, rank, original_annotation, reviewed_metabolite, rename_from_comment, 
        rename_from_observation, rename_source, corrected_annotation_raw, corrected_metabolite, canonical_metabolite_key, 
        rank2_candidates, rank2_canonical_candidate_count, rank2_conflict_resolved_by_canonicalization) %>% 
    dplyr::arrange(corrected_metabolite, feature_key)
manual_cluster_audit <- representative_table %>% dplyr::filter(manual_cluster_applied) %>% dplyr::select(feature_key, 
    corrected_metabolite, mz, rt_min, manual_analyte_cluster, manual_cluster_reason, selected_analyte_representative) %>% 
    dplyr::arrange(manual_analyte_cluster, feature_key)
summary_table <- tibble::tibble(item = c("Samples total", "Caecum samples", "Faeces samples", "Independent donors", 
    "Control donors", "PPI donors", "Rank1/Rank2 original features reconstructed", "Included original features after review", 
    "RT-resolved analytes", "One-per-metabolite features", "Metabolites with >1 retained RT cluster", 
    "Manual review-status overrides active", "RT tolerance (min)", "Treatment statistic used for selection"), 
    value = c(nrow(metadata_aligned), sum(metadata_aligned$Organ == "cecum"), sum(metadata_aligned$Organ == 
        "feces"), dplyr::n_distinct(metadata_aligned$Donor), dplyr::n_distinct(metadata_aligned$Donor[metadata_aligned$Separation == 
        "ctrl"]), dplyr::n_distinct(metadata_aligned$Donor[metadata_aligned$Separation == "ppi"]), ncol(X_all), 
        nrow(included_features), ncol(rt_resolved_all), ncol(one_per_metabolite_all), dplyr::n_distinct(multiple_rt_cluster_choices$corrected_metabolite), 
        sum(review$status_was_overridden), rt_tolerance_minutes, "No"))
audit_workbook_path <- file.path(out_dir, "Cecum_Feces_LCMS_collaborator_review_audit.xlsx")
writexl::write_xlsx(list(Summary = summary_table, Metadata = metadata_export, Collaborator_review = review, 
    Feature_decisions = feature_decisions, Name_harmonization = name_harmonization_audit, Included_features = included_features, 
    Analyte_cluster_audit = representative_table, Selected_RT_analytes = selected_analytes, Metabolite_selection = metabolite_representative_audit, 
    Selected_one_per_met = selected_metabolites, Multiple_RT_choices = multiple_rt_cluster_choices, Manual_cluster_audit = manual_cluster_audit, 
    Post_rename_duplicates = post_rename_duplicates, Rank2_conflicts = rank2_annotation_conflicts, Excluded_features = excluded_features), 
    path = ppi_output(audit_workbook_path))
readr::write_csv(matrix_export(rt_cecum$matrix), ppi_output(file.path(out_dir, "Cecum_LCMS_cleaned_RT_resolved_matrix.csv"), 
    name = "caecum_lcms_cleaned_rt_resolved_matrix.csv"))
readr::write_csv(matrix_export(rt_feces$matrix), ppi_output(file.path(out_dir, "Feces_LCMS_cleaned_RT_resolved_matrix.csv"), 
    name = "faeces_lcms_cleaned_rt_resolved_matrix.csv"))
readr::write_csv(matrix_export(met_cecum$matrix), ppi_output(file.path(out_dir, "Cecum_LCMS_cleaned_one_per_metabolite_matrix.csv"), 
    name = "caecum_lcms_cleaned_one_per_metabolite_matrix.csv"))
readr::write_csv(matrix_export(met_feces$matrix), ppi_output(file.path(out_dir, "Feces_LCMS_cleaned_one_per_metabolite_matrix.csv"), 
    name = "faeces_lcms_cleaned_one_per_metabolite_matrix.csv"))
readr::write_csv(metadata_export, ppi_output(file.path(out_dir, "Cecum_Feces_LCMS_limma_metadata.csv"), 
    name = "caecum_faeces_lcms_limma_metadata.csv"))
readr::write_csv(selected_analytes, ppi_output(file.path(out_dir, "Cecum_Feces_LCMS_selected_RT_resolved_analytes.csv"), 
    name = "caecum_faeces_lcms_selected_rt_resolved_analytes.csv"))
readr::write_csv(selected_metabolites, ppi_output(file.path(out_dir, "Cecum_Feces_LCMS_selected_one_per_metabolite.csv"), 
    name = "caecum_faeces_lcms_selected_one_per_metabolite.csv"))
main_rds_path <- file.path(out_dir, "Cecum_Feces_LCMS_cleaned_analysis_objects.rds")
saveRDS(list(settings = list(rt_tolerance_minutes = rt_tolerance_minutes, exclude_unresolved_rank2_conflicts = exclude_unresolved_rank2_conflicts, 
    manual_exclude_metabolites = manual_exclude_metabolites, manual_feature_cluster_overrides = manual_feature_cluster_overrides, 
    note = "No treatment-association statistic was used for feature selection. Metabolites were grouped using canonical biochemical names."), 
    review = review, feature_decisions = feature_decisions, included_features = included_features, representative_table = representative_table, 
    selected_analytes = selected_analytes, metabolite_representative_audit = metabolite_representative_audit, 
    selected_metabolites = selected_metabolites, metadata_all = metadata_aligned, X_all_original_features = X_all, 
    rt_resolved_all = rt_resolved_all, one_per_metabolite_all = one_per_metabolite_all, rt_resolved = list(cecum = rt_cecum$matrix, 
        feces = rt_feces$matrix), one_per_metabolite = list(cecum = met_cecum$matrix, feces = met_feces$matrix), 
    metadata = list(cecum = rt_cecum$metadata, feces = rt_feces$metadata), summary = summary_table), 
    ppi_output(main_rds_path))
limma_rds_path <- file.path(out_dir, "Cecum_Feces_LCMS_limma_input_objects.rds")
saveRDS(list(rt_resolved = list(cecum = rt_cecum$matrix, feces = rt_feces$matrix), one_per_metabolite = list(cecum = met_cecum$matrix, 
    feces = met_feces$matrix), metadata = list(cecum = rt_cecum$metadata, feces = rt_feces$metadata), 
    selected_analytes = selected_analytes, selected_metabolites = selected_metabolites, input_notes = c("Rows are recipient-mouse samples; columns are curated LC-MS analytes/metabolites.", 
        "rt_resolved retains distinct same-name RT clusters and is recommended as the primary analytical representation.", 
        "one_per_metabolite selects one RT/analyte cluster per canonical biochemical name and is a secondary representation.", 
        "No log transformation is applied in this RDS.", "The limma script applies feature-specific log2(intensity + half-minimum-positive).", 
        "Donor, Mouse_ID, Donor_Mouse_ID, Organ, and Separation are explicit metadata columns.", "Separation is coded ctrl then ppi.", 
        "No sample-weight normalization is performed.")), ppi_output(limma_rds_path))
message("\nCaecum/faeces LC-MS cleaning completed.\n", "Samples: ", nrow(metadata_aligned), " (", sum(metadata_aligned$Organ == 
    "cecum"), " caecum; ", sum(metadata_aligned$Organ == "feces"), " faeces)", "\nIndependent donors: ", 
    dplyr::n_distinct(metadata_aligned$Donor), "\nOriginal Rank1/Rank2 features reconstructed: ", ncol(X_all), 
    "\nOriginal features retained after review: ", nrow(included_features), "\nRT-resolved analytes: ", 
    ncol(rt_resolved_all), "\nOne-per-metabolite features: ", ncol(one_per_metabolite_all), "\nMetabolites with multiple retained RT clusters: ", 
    dplyr::n_distinct(multiple_rt_cluster_choices$corrected_metabolite), "\nManual review-status overrides active: ", 
    sum(review$status_was_overridden), "\n\nAudit workbook:\n", audit_workbook_path, "\n\nMain RDS:\n", 
    main_rds_path, "\n\nLimma input RDS:\n", limma_rds_path)
