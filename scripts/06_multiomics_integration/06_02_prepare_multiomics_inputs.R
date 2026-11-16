# SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0
# Prepare inputs
# Run through run.R so input and output paths resolve within the package.
# Inputs: see data/input_manifest.csv and docs/RUN_ORDER.md.
# Writes generated files under results/; bundled inputs are read-only.

suppressPackageStartupMessages({
    library(dplyr)
    library(tidyr)
    library(purrr)
    library(readr)
    library(readxl)
    library(stringr)
    library(tibble)
    library(phyloseq)
    library(writexl)
})
ppi_setwd(ppi_path("integration"))
script_path <- tryCatch(ppi_normalize_path(sys.frame(1)$ofile, winslash = "/", mustWork = FALSE), error = function(e) NA_character_)
script_dir <- if (!is.na(script_path) && nzchar(script_path)) {
    dirname(script_path)
} else {
    getwd()
}
config_path <- ppi_code("scripts/06_multiomics_integration/06_01_configure_multiomics.R")
if (!ppi_file_exists(config_path)) {
    config_path <- "06_01_configure_multiomics.R"
}
ppi_source(config_path)
set.seed(20260715)
dir.create(ppi_output(association_root, directory = TRUE), recursive = TRUE, showWarnings = FALSE)
prep_dir <- file.path(association_root, "01_prepared_inputs")
dir.create(ppi_output(prep_dir, directory = TRUE), recursive = TRUE, showWarnings = FALSE)
prepared_output_file <- file.path(prep_dir, "PPI_multiomics_prepared_inputs.rds")
if (ppi_file_exists(prepared_output_file)) {
    unlink(prepared_output_file)
}
`%||%` <- function(x, y) {
    if (is.null(x) || length(x) == 0L) 
        y
    else x
}
normalize_search_basename <- function(x) {
    ext <- tolower(tools::file_ext(x))
    stem <- tools::file_path_sans_ext(basename(x))
    stem <- stringr::str_remove(stem, " ?\\([0-9]+\\)$")
    paste0(tolower(stem), ifelse(nzchar(ext), paste0(".", ext), ""))
}
resolve_input_file <- function(path, search_root = project_dir, required = TRUE) {
    if (!is.null(path) && length(path) == 1L && !is.na(path) && ppi_file_exists(path)) {
        return(ppi_normalize_path(path, winslash = "/", mustWork = TRUE))
    }
    if (is.null(path) || length(path) != 1L || is.na(path) || !nzchar(path)) {
        if (required) 
            stop("A required input path was not supplied.")
        return(NA_character_)
    }
    if (isTRUE(strict_input_paths)) {
        if (required) {
            stop("Required input file not found at the exact configured path:\n", path, "\n\nEdit 06_01_configure_multiomics.R so that it ", 
                "points to the final reviewed file. Recursive filename fallback is ", "disabled to prevent an older copy being selected silently.")
        }
        return(NA_character_)
    }
    target <- basename(path)
    all_files <- ppi_list_files(search_root, recursive = TRUE, full.names = TRUE, include.dirs = FALSE)
    hits <- all_files[normalize_search_basename(all_files) == normalize_search_basename(target)]
    if (length(hits) == 0L) {
        if (required) {
            stop("Required input file not found: ", path, "\nNo matching basename was found recursively under: ", 
                search_root)
        }
        return(NA_character_)
    }
    if (length(hits) > 1L) {
        stop("Multiple files matched ", target, ". Set the exact final path in the ", "configuration file. Matches:\n", 
            paste(hits, collapse = "\n"))
    }
    ppi_normalize_path(hits[[1]], winslash = "/", mustWork = TRUE)
}
first_existing_column <- function(data, candidates, required = TRUE) {
    hit <- candidates[candidates %in% names(data)]
    if (length(hit) > 0L) 
        return(hit[[1]])
    if (required) {
        stop("None of the required columns were found: ", paste(candidates, collapse = ", "), "\nAvailable columns: ", 
            paste(names(data), collapse = ", "))
    }
    NULL
}
normalize_group <- function(x) {
    x <- str_to_lower(str_squish(as.character(x)))
    case_when(x %in% c("ctrl", "control", "ctlr", "c", "0", "no") ~ "ctrl", x %in% c("ppi", "p", "1", 
        "yes") ~ "ppi", str_detect(x, "^ctrl|^control|^ctlr|^c[0-9]+$") ~ "ctrl", str_detect(x, "^ppi|^p[0-9]+$") ~ 
        "ppi", TRUE ~ NA_character_)
}
normalize_donor <- function(x) {
    raw <- str_to_lower(str_squish(as.character(x)))
    raw <- str_replace_all(raw, "control|ctlr", "ctrl")
    m <- str_match(raw, "(ctrl|ppi|c|p)[_-]?([0-9]+)")
    prefix <- case_when(m[, 2] %in% c("ctrl", "c") ~ "ctrl", m[, 2] %in% c("ppi", "p") ~ "ppi", TRUE ~ 
        NA_character_)
    if_else(is.na(prefix), NA_character_, paste0(prefix, m[, 3]))
}
parse_mouse_identity <- function(x) {
    raw <- str_to_lower(str_squish(as.character(x)))
    raw <- str_replace_all(raw, "control|ctlr", "ctrl")
    m <- str_match(raw, "(?:^|.*?)(ctrl|ppi)[_-]?([0-9]+)[_-]m[_-]?([0-9]+)(?:$|[^0-9].*$)")
    donor <- if_else(is.na(m[, 2]), NA_character_, paste0(m[, 2], m[, 3]))
    mouse_number <- suppressWarnings(as.integer(m[, 4]))
    tibble(Group = factor(if_else(is.na(donor), NA_character_, if_else(str_starts(donor, "ppi"), "ppi", 
        "ctrl")), levels = c("ctrl", "ppi")), Donor = donor, Mouse_Number = mouse_number, Mouse_Key = if_else(is.na(donor) | 
        is.na(mouse_number), NA_character_, paste0(donor, "_mouse", mouse_number)))
}
parse_culture_identity <- function(x) {
    raw <- str_to_lower(str_squish(as.character(x)))
    raw <- str_replace_all(raw, "control|ctlr", "ctrl")
    m <- str_match(raw, "(?:^|.*?)(ctrl|ppi)[_-]?([0-9]+)[_-]cp[_-]?d[_-]?([0-9]+)(?:$|[^0-9].*$)")
    donor <- if_else(is.na(m[, 2]), NA_character_, paste0(m[, 2], m[, 3]))
    day <- suppressWarnings(as.integer(m[, 4]))
    tibble(Group = factor(if_else(is.na(donor), NA_character_, if_else(str_starts(donor, "ppi"), "ppi", 
        "ctrl")), levels = c("ctrl", "ppi")), Donor = donor, Day = day, DayF = factor(day, levels = settings$exvivo_days), 
        Culture_Key = if_else(is.na(donor) | is.na(day), NA_character_, paste0(donor, "_d", day)))
}
parse_original_donor_identity <- function(x) {
    donor <- normalize_donor(x)
    tibble(Group = factor(if_else(is.na(donor), NA_character_, if_else(str_starts(donor, "ppi"), "ppi", 
        "ctrl")), levels = c("ctrl", "ppi")), Donor = donor, Donor_Key = if_else(is.na(donor), NA_character_, 
        paste0(donor, "_original")))
}
calculate_half_minimum <- function(x) {
    positive <- x[is.finite(x) & x > 0]
    if (length(positive) == 0L) 
        return(NA_real_)
    min(positive, na.rm = TRUE)/2
}
log2_half_min_matrix <- function(X, dataset_id) {
    X <- as.matrix(X)
    storage.mode(X) <- "numeric"
    X[!is.finite(X)] <- NA_real_
    X[X < 0] <- 0
    feature_names <- colnames(X)
    if (is.null(feature_names)) 
        feature_names <- paste0("feature_", seq_len(ncol(X)))
    n_positive <- colSums(X > 0, na.rm = TRUE)
    missing_fraction <- colMeans(is.na(X))
    pseudocount <- vapply(seq_len(ncol(X)), function(j) calculate_half_minimum(X[, j]), numeric(1))
    minimum_positive_required <- max(as.integer(settings$metabolite_minimum_positive_samples), as.integer(ceiling(nrow(X) * 
        settings$metabolite_minimum_positive_fraction)))
    preliminary_keep <- n_positive >= minimum_positive_required & missing_fraction <= settings$metabolite_maximum_missing_fraction & 
        is.finite(pseudocount) & pseudocount > 0
    X_log <- matrix(NA_real_, nrow = nrow(X), ncol = sum(preliminary_keep), dimnames = list(rownames(X), 
        feature_names[preliminary_keep]))
    source_columns <- which(preliminary_keep)
    for (k in seq_along(source_columns)) {
        j <- source_columns[[k]]
        values <- log2(X[, j] + pseudocount[[j]])
        if (anyNA(values)) {
            replacement <- median(values, na.rm = TRUE)
            if (!is.finite(replacement)) 
                replacement <- 0
            values[is.na(values)] <- replacement
        }
        X_log[, k] <- values
    }
    variance_log <- if (ncol(X_log) > 0L) {
        apply(X_log, 2, var, na.rm = TRUE)
    }
    else {
        numeric()
    }
    variance_keep <- is.finite(variance_log) & variance_log > 0
    variance_by_original_feature <- setNames(rep(NA_real_, length(feature_names)), feature_names)
    if (length(variance_log) > 0L) {
        variance_by_original_feature[names(variance_log)] <- as.numeric(variance_log)
    }
    X_log <- X_log[, variance_keep, drop = FALSE]
    audit <- tibble(dataset_id = dataset_id, metabolite = feature_names, n_samples = nrow(X), minimum_positive_required = minimum_positive_required, 
        n_positive = n_positive, detected_fraction = n_positive/nrow(X), missing_fraction = missing_fraction, 
        pseudocount = pseudocount, variance = unname(variance_by_original_feature[feature_names]), passes_usability_filter = preliminary_keep, 
        retained_final = feature_names %in% colnames(X_log), exclusion_reason = case_when(n_positive < 
            minimum_positive_required ~ "too_few_positive_samples", missing_fraction > settings$metabolite_maximum_missing_fraction ~ 
            "too_many_missing_values", !is.finite(pseudocount) | pseudocount <= 0 ~ "no_positive_measurement", 
            preliminary_keep & !(feature_names %in% colnames(X_log)) ~ "zero_variance_after_transform", 
            TRUE ~ NA_character_))
    list(matrix = X_log, audit = audit)
}
prepare_score_matrix <- function(X, dataset_id) {
    X <- as.matrix(X)
    storage.mode(X) <- "numeric"
    feature_names <- colnames(X)
    if (is.null(feature_names)) 
        feature_names <- paste0("feature_", seq_len(ncol(X)))
    missing_fraction <- colMeans(!is.finite(X))
    for (j in seq_len(ncol(X))) {
        bad <- !is.finite(X[, j])
        if (any(bad)) {
            replacement <- median(X[, j][!bad], na.rm = TRUE)
            if (!is.finite(replacement)) 
                replacement <- 0
            X[bad, j] <- replacement
        }
    }
    v <- apply(X, 2, var, na.rm = TRUE)
    keep <- missing_fraction <= settings$metabolite_maximum_missing_fraction & is.finite(v) & v > 0
    audit <- tibble(dataset_id = dataset_id, metabolite = feature_names, missing_fraction = missing_fraction, 
        variance = v, retained_final = keep, exclusion_reason = case_when(missing_fraction > settings$metabolite_maximum_missing_fraction ~ 
            "too_many_missing_values", !is.finite(v) | v <= 0 ~ "zero_variance", TRUE ~ NA_character_))
    list(matrix = X[, keep, drop = FALSE], audit = audit)
}
align_matrix_metadata <- function(X, metadata, sample_column) {
    X <- as.matrix(X)
    metadata <- as.data.frame(metadata, check.names = FALSE)
    sample_ids <- as.character(metadata[[sample_column]])
    if (!is.null(rownames(X)) && all(rownames(X) %in% sample_ids)) {
        metadata <- metadata[match(rownames(X), sample_ids), , drop = FALSE]
    }
    else if (nrow(X) == nrow(metadata)) {
        rownames(X) <- sample_ids
    }
    else {
        stop("Could not align matrix and metadata using column: ", sample_column)
    }
    metadata$Sample_ID <- rownames(X)
    rownames(metadata) <- metadata$Sample_ID
    list(matrix = X, metadata = metadata)
}
resolved_paths <- list(phyloseq = resolve_input_file(input_paths$phyloseq), pathway_scores = resolve_input_file(input_paths$pathway_scores), 
    rna_loader_script = resolve_input_file(input_paths$rna_loader_script), rna_metadata_table = resolve_input_file(input_paths$rna_metadata_table, 
        required = FALSE), serum_lcms = resolve_input_file(input_paths$serum_lcms), cecum_feces_lcms = resolve_input_file(input_paths$cecum_feces_lcms), 
    cecum_nmr = resolve_input_file(input_paths$cecum_nmr), exvivo_nmr_raw = resolve_input_file(input_paths$exvivo_nmr_raw), 
    exvivo_nmr_mapping = resolve_input_file(input_paths$exvivo_nmr_mapping), exvivo_hilic = resolve_input_file(input_paths$exvivo_hilic))
resolved_gsea_paths <- purrr::map_chr(gsea_paths, resolve_input_file)
path_audit <- bind_rows(enframe(unlist(resolved_paths), name = "input", value = "resolved_path"), enframe(resolved_gsea_paths, 
    name = "input", value = "resolved_path") %>% mutate(input = paste0("GSEA_", input))) %>% mutate(file_exists = vapply(resolved_path, 
    function(p) !is.na(p) && nzchar(p) && ppi_file_exists(p), logical(1)), size_bytes = vapply(resolved_path, 
    function(p) {
        if (is.na(p) || !nzchar(p) || !ppi_file_exists(p)) 
            return(NA_real_)
        as.numeric(file.info(ppi_input(p))$size)
    }, numeric(1)), modified_time = vapply(resolved_path, function(p) {
    if (is.na(p) || !nzchar(p) || !ppi_file_exists(p)) 
        return(NA_character_)
    format(file.info(ppi_input(p))$mtime, "%Y-%m-%d %H:%M:%S")
}, character(1)), md5 = vapply(resolved_path, function(p) {
    if (is.na(p) || !nzchar(p) || !ppi_file_exists(p)) 
        return(NA_character_)
    unname(as.character(tools::md5sum(p)))
}, character(1)))
write_csv(path_audit, ppi_output(file.path(prep_dir, "resolved_input_paths.csv"), name = "resolved_input_paths.csv"))
collect_rna_metadata_from_loader <- function(loader_path) {
    env <- new.env(parent = globalenv())
    ppi_sys_source(loader_path, envir = env)
    candidate_objects <- c("meta", "meta_liver", "meta_colon", "meta_brain")
    frames <- list()
    for (nm in candidate_objects) {
        if (exists(nm, envir = env, inherits = FALSE)) {
            value <- get(nm, envir = env, inherits = FALSE)
            if (is.data.frame(value)) 
                frames[[nm]] <- as.data.frame(value)
        }
    }
    if (length(frames) == 0L) {
        stop("The RNA loader ran, but none of these metadata objects were found: ", paste(candidate_objects, 
            collapse = ", "))
    }
    bind_rows(frames, .id = "metadata_source")
}
rna_meta_raw <- if (!is.na(resolved_paths$rna_metadata_table)) {
    if (tolower(tools::file_ext(resolved_paths$rna_metadata_table)) %in% c("xlsx", "xls")) {
        read_excel(ppi_input(resolved_paths$rna_metadata_table))
    }
    else {
        read_csv(ppi_input(resolved_paths$rna_metadata_table), show_col_types = FALSE)
    }
} else {
    collect_rna_metadata_from_loader(resolved_paths$rna_loader_script)
}
rna_sample_col <- first_existing_column(rna_meta_raw, c("name", "sample", "sample_id", "Sample"))
rna_external_col <- first_existing_column(rna_meta_raw, c("external_name", "external.name", "external", 
    "mouse_external_name"))
rna_organ_col <- first_existing_column(rna_meta_raw, c("organ", "Organ", "tissue"), required = FALSE)
rna_group_col <- first_existing_column(rna_meta_raw, c("group", "Group", "Separation"), required = FALSE)
rna_donor_col <- first_existing_column(rna_meta_raw, c("donor", "Donor"), required = FALSE)
rna_identity <- parse_mouse_identity(rna_meta_raw[[rna_external_col]])
rna_bridge <- tibble(RNA_Sample = as.character(rna_meta_raw[[rna_sample_col]]), RNA_External_Name = as.character(rna_meta_raw[[rna_external_col]]), 
    Organ = if (!is.null(rna_organ_col)) as.character(rna_meta_raw[[rna_organ_col]]) else NA_character_, 
    Group_loader = if (!is.null(rna_group_col)) normalize_group(rna_meta_raw[[rna_group_col]]) else NA_character_, 
    Donor_loader = if (!is.null(rna_donor_col)) normalize_donor(rna_meta_raw[[rna_donor_col]]) else NA_character_) %>% 
    bind_cols(rna_identity) %>% mutate(Organ = str_to_title(Organ), Group = factor(coalesce(as.character(Group), 
    Group_loader), levels = c("ctrl", "ppi")), Donor = coalesce(Donor, Donor_loader)) %>% select(RNA_Sample, 
    RNA_External_Name, Organ, Group, Donor, Mouse_Number, Mouse_Key) %>% filter(!is.na(RNA_Sample), RNA_Sample != 
    "") %>% distinct(RNA_Sample, .keep_all = TRUE)
if (anyNA(rna_bridge$Mouse_Key)) {
    warning(sum(is.na(rna_bridge$Mouse_Key)), " RNA samples could not be mapped to a mouse key.")
}
pathway_scores_raw <- read_csv(ppi_input(resolved_paths$pathway_scores), show_col_types = FALSE)
score_sample_col <- first_existing_column(pathway_scores_raw, c("sample", "name", "RNA_Sample"))
score_organ_col <- first_existing_column(pathway_scores_raw, c("organ", "Organ", "tissue"))
score_group_col <- first_existing_column(pathway_scores_raw, c("group", "Group"), required = FALSE)
score_donor_col <- first_existing_column(pathway_scores_raw, c("donor", "Donor"), required = FALSE)
pathway_scores_raw <- pathway_scores_raw %>% mutate(RNA_Sample = as.character(.data[[score_sample_col]]), 
    Organ_Score = str_to_title(as.character(.data[[score_organ_col]])), Group_Score = if (!is.null(score_group_col)) normalize_group(.data[[score_group_col]]) else NA_character_, 
    Donor_Score = if (!is.null(score_donor_col)) normalize_donor(.data[[score_donor_col]]) else NA_character_) %>% 
    left_join(rna_bridge, by = "RNA_Sample") %>% mutate(Organ = coalesce(Organ, Organ_Score), Group = factor(coalesce(as.character(Group), 
    Group_Score), levels = c("ctrl", "ppi")), Donor = coalesce(Donor, Donor_Score))
read_gsea_selection <- function(path, label) {
    parts <- str_split_fixed(label, "_", 2)
    organ <- str_to_title(parts[, 1])
    library <- toupper(parts[, 2])
    dat <- read_csv(ppi_input(path), show_col_types = FALSE)
    pathway_col <- first_existing_column(dat, c("pathway", "ID", "Description"))
    padj_col <- first_existing_column(dat, c("padj", "p.adjust", "FDR"))
    nes_col <- first_existing_column(dat, c("NES", "nes"), required = FALSE)
    pval_col <- first_existing_column(dat, c("pval", "pvalue", "p.value"), required = FALSE)
    tibble(Organ = organ, Library = library, Pathway = as.character(dat[[pathway_col]]), GSEA_NES = if (!is.null(nes_col)) 
        as.numeric(dat[[nes_col]])
    else NA_real_, GSEA_P = if (!is.null(pval_col)) 
        as.numeric(dat[[pval_col]])
    else NA_real_, GSEA_FDR = as.numeric(dat[[padj_col]])) %>% filter(is.finite(GSEA_FDR), GSEA_FDR < 
        settings$gsea_fdr_cutoff) %>% distinct(Organ, Library, Pathway, .keep_all = TRUE)
}
pathway_selection <- imap_dfr(as.list(resolved_gsea_paths), ~read_gsea_selection(.x, .y)) %>% mutate(Present_in_score_table = Pathway %in% 
    names(pathway_scores_raw))
if (any(!pathway_selection$Present_in_score_table)) {
    warning(sum(!pathway_selection$Present_in_score_table), " significant GSEA pathways were absent from the ssGSEA score table.")
}
pathway_selection_present <- pathway_selection %>% filter(Present_in_score_table)
build_pathway_dataset <- function(organ, library) {
    selected <- pathway_selection_present %>% filter(.data$Organ == .env$organ, .data$Library == .env$library)
    pathway_names <- intersect(selected$Pathway, names(pathway_scores_raw))
    meta <- pathway_scores_raw %>% filter(.data$Organ == .env$organ, !is.na(.data$Mouse_Key)) %>% transmute(Sample_ID = RNA_Sample, 
        RNA_External_Name, Organ, Library = library, Group, Donor, Mouse_Number, Mouse_Key)
    X <- pathway_scores_raw %>% filter(.data$Organ == .env$organ, !is.na(.data$Mouse_Key)) %>% select(all_of(pathway_names)) %>% 
        as.matrix()
    storage.mode(X) <- "numeric"
    rownames(X) <- meta$Mouse_Key
    keep <- vapply(seq_len(ncol(X)), function(j) {
        x <- X[, j]
        sum(is.finite(x)) >= settings$minimum_samples_mouse && is.finite(var(x, na.rm = TRUE)) && var(x, 
            na.rm = TRUE) > 0
    }, logical(1))
    X <- X[, keep, drop = FALSE]
    for (j in seq_len(ncol(X))) {
        bad <- !is.finite(X[, j])
        if (any(bad)) 
            X[bad, j] <- median(X[, j][!bad], na.rm = TRUE)
    }
    list(matrix = X, metadata = meta, selection = selected %>% filter(Pathway %in% colnames(X)), key_col = "Mouse_Key", 
        organ = organ, library = library)
}
pathway_datasets <- list()
for (organ in c("Liver", "Colon", "Brain")) {
    for (library in c("HALLMARK", "REACTOME")) {
        id <- paste(organ, library, sep = "__")
        pathway_datasets[[id]] <- build_pathway_dataset(organ, library)
    }
}
write_csv(rna_bridge, ppi_output(file.path(prep_dir, "RNA_sample_to_mouse_bridge.csv"), name = "rna_sample_to_mouse_bridge.csv"))
write_csv(pathway_selection, ppi_output(file.path(prep_dir, "significant_pathway_selection_and_matching.csv"), 
    name = "significant_pathway_selection_and_matching.csv"))
load_phyloseq_any <- function(path) {
    ext <- tolower(tools::file_ext(path))
    if (ext == "rds") {
        obj <- readRDS(ppi_input(path))
        if (inherits(obj, "phyloseq")) 
            return(obj)
        if (is.list(obj)) {
            hit <- which(vapply(obj, inherits, logical(1), what = "phyloseq"))
            if (length(hit) >= 1L) 
                return(obj[[hit[[1]]]])
        }
        stop("The RDS file did not contain a phyloseq object.")
    }
    env <- new.env(parent = emptyenv())
    loaded <- load(ppi_input(path), envir = env)
    hit <- loaded[vapply(loaded, function(nm) inherits(get(nm, envir = env), "phyloseq"), logical(1))]
    if (length(hit) == 0L) 
        stop("No phyloseq object was found in: ", path)
    get(hit[[1]], envir = env)
}
ps <- load_phyloseq_any(resolved_paths$phyloseq)
otu <- as(otu_table(ps), "matrix")
if (taxa_are_rows(ps)) otu <- t(otu)
storage.mode(otu) <- "numeric"
sample_meta_16s <- as(sample_data(ps), "data.frame") %>% rownames_to_column("Phyloseq_Sample_ID")
otu <- otu[match(sample_meta_16s$Phyloseq_Sample_ID, rownames(otu)), , drop = FALSE]
rownames(otu) <- sample_meta_16s$Phyloseq_Sample_ID
rank_name <- settings$taxonomic_rank
tax <- as(tax_table(ps), "matrix")
if (!rank_name %in% colnames(tax)) {
    stop("Taxonomic rank not found in phyloseq taxonomy table: ", rank_name)
}
tax <- tax[match(colnames(otu), rownames(tax)), , drop = FALSE]
genus <- str_squish(as.character(tax[, rank_name]))
valid_genus <- !is.na(genus) & genus != "" & !str_detect(str_to_lower(genus), "uncultured|unclassified|unknown|unassigned|^na$")
genus_counts <- rowsum(t(otu[, valid_genus, drop = FALSE]), group = genus[valid_genus], reorder = FALSE) %>% 
    t()
sample_name_col_16s <- first_existing_column(sample_meta_16s, c("sample_name", "Sample_Name", "sample", 
    "Sample", "biological_sample"))
experiment_col_16s <- first_existing_column(sample_meta_16s, c("experiment", "Experiment", "sample_type", 
    "Sample_Type"), required = FALSE)
sample_name_16s <- as.character(sample_meta_16s[[sample_name_col_16s]])
experiment_16s <- if (!is.null(experiment_col_16s)) {
    str_to_lower(str_squish(as.character(sample_meta_16s[[experiment_col_16s]])))
} else {
    rep(NA_character_, nrow(sample_meta_16s))
}
mouse_identity_16s <- parse_mouse_identity(sample_name_16s)
culture_identity_16s <- parse_culture_identity(sample_name_16s)
original_donor_identity_16s <- parse_original_donor_identity(sample_name_16s)
both_identity <- !is.na(mouse_identity_16s$Mouse_Key) & !is.na(culture_identity_16s$Culture_Key)
if (any(both_identity)) {
    stop("Some 16S sample names parsed simultaneously as mouse and culture: ", paste(sample_name_16s[both_identity], 
        collapse = ", "))
}
context_from_metadata <- case_when(experiment_16s == "culture" ~ "culture", experiment_16s %in% c("mouse", 
    "mice", "recipient", "recipient_mouse", "in_vivo") ~ "mouse", experiment_16s %in% c("og_sample", 
    "og sample", "original_sample", "original sample", "original_donor", "original donor", "donor") ~ 
    "donor", TRUE ~ NA_character_)
context_from_name <- case_when(!is.na(culture_identity_16s$Culture_Key) ~ "culture", !is.na(mouse_identity_16s$Mouse_Key) ~ 
    "mouse", TRUE ~ "unmatched")
metadata_name_disagreement <- context_from_metadata %in% c("mouse", "culture") & context_from_name != 
    "unmatched" & context_from_metadata != context_from_name
if (any(metadata_name_disagreement)) {
    stop("The 16S experiment metadata disagrees with the biological sample name for: ", paste(sample_name_16s[metadata_name_disagreement], 
        collapse = ", "))
}
sample_identity_16s <- sample_meta_16s %>% transmute(Phyloseq_Sample_ID, Original_16S_Sample_Name = sample_name_16s, 
    Experiment = experiment_16s, Context_from_metadata = context_from_metadata, Context_from_name = context_from_name) %>% 
    bind_cols(mouse_identity_16s %>% rename_with(~paste0("Mouse_", .x)), culture_identity_16s %>% rename_with(~paste0("Culture_", 
        .x)), original_donor_identity_16s %>% rename_with(~paste0("OriginalDonor_", .x))) %>% mutate(Biological_Context = case_when(context_from_metadata == 
    "donor" ~ "donor", context_from_metadata == "culture" & !is.na(Culture_Culture_Key) ~ "culture", 
    context_from_metadata == "mouse" & !is.na(Mouse_Mouse_Key) ~ "mouse", context_from_name == "culture" ~ 
        "culture", context_from_name == "mouse" ~ "mouse", TRUE ~ "unmatched"), Biological_Key = case_when(Biological_Context == 
    "mouse" ~ Mouse_Mouse_Key, Biological_Context == "culture" ~ Culture_Culture_Key, Biological_Context == 
    "donor" ~ dplyr::coalesce(OriginalDonor_Donor_Key, paste0("original_", Phyloseq_Sample_ID)), TRUE ~ 
    NA_character_), Group = factor(case_when(Biological_Context == "mouse" ~ as.character(Mouse_Group), 
    Biological_Context == "culture" ~ as.character(Culture_Group), Biological_Context == "donor" ~ as.character(OriginalDonor_Group), 
    TRUE ~ NA_character_), levels = c("ctrl", "ppi")), Donor = case_when(Biological_Context == "mouse" ~ 
    Mouse_Donor, Biological_Context == "culture" ~ Culture_Donor, Biological_Context == "donor" ~ OriginalDonor_Donor, 
    TRUE ~ NA_character_), Mouse_Number = if_else(Biological_Context == "mouse", Mouse_Mouse_Number, 
    NA_integer_), Day = if_else(Biological_Context == "culture", Culture_Day, NA_integer_), DayF = factor(Day, 
    levels = settings$exvivo_days))
prepare_genus_context <- function(context_name, keep_rows) {
    counts <- genus_counts[keep_rows, , drop = FALSE]
    meta <- sample_identity_16s[keep_rows, , drop = FALSE]
    valid_key <- !is.na(meta$Biological_Key) & meta$Biological_Key != ""
    counts <- counts[valid_key, , drop = FALSE]
    meta <- meta[valid_key, , drop = FALSE]
    counts_agg <- rowsum(counts, group = meta$Biological_Key, reorder = FALSE)
    key_order <- rownames(counts_agg)
    meta_agg <- meta %>% filter(Biological_Key %in% key_order) %>% arrange(match(Biological_Key, key_order)) %>% 
        distinct(Biological_Key, .keep_all = TRUE)
    total <- rowSums(counts_agg)
    keep_depth <- is.finite(total) & total > 0
    counts_agg <- counts_agg[keep_depth, , drop = FALSE]
    meta_agg <- meta_agg[keep_depth, , drop = FALSE]
    total <- total[keep_depth]
    relative <- sweep(counts_agg, 1, total, "/")
    prevalence <- colMeans(counts_agg > 0)
    mean_relative <- colMeans(relative, na.rm = TRUE)
    donor_presence <- vapply(seq_len(ncol(counts_agg)), function(j) {
        n_distinct(meta_agg$Donor[counts_agg[, j] > 0])
    }, integer(1))
    keep_genus <- prevalence >= settings$genus_prevalence_fraction & mean_relative >= settings$genus_mean_relative_abundance & 
        donor_presence >= settings$genus_minimum_donors
    counts_filtered <- counts_agg[, keep_genus, drop = FALSE]
    log_counts <- log(counts_filtered + settings$genus_pseudocount)
    clr <- sweep(log_counts, 1, rowMeans(log_counts), "-")
    meta_out <- meta_agg %>% transmute(Sample_ID = Biological_Key, Biological_Context = context_name, 
        Biological_Key, Group, Donor, Mouse_Number, Mouse_Key = Mouse_Mouse_Key, Day, DayF, Culture_Key = Culture_Culture_Key, 
        n_technical_repeats = as.integer(table(meta$Biological_Key)[Biological_Key]))
    rownames(clr) <- meta_out$Sample_ID
    audit <- tibble(context = context_name, Genus = colnames(counts_agg), prevalence = prevalence, mean_relative_abundance = mean_relative, 
        n_donors_present = donor_presence, retained = keep_genus, exclusion_reason = case_when(prevalence < 
            settings$genus_prevalence_fraction ~ "low_prevalence", mean_relative < settings$genus_mean_relative_abundance ~ 
            "low_mean_relative_abundance", donor_presence < settings$genus_minimum_donors ~ "too_few_donors", 
            TRUE ~ NA_character_))
    list(matrix = clr, counts = counts_filtered, metadata = meta_out, audit = audit, key_col = if (context_name == 
        "mouse") "Mouse_Key" else "Culture_Key")
}
mouse_keep <- sample_identity_16s$Biological_Context == "mouse"
culture_keep <- sample_identity_16s$Biological_Context == "culture" & sample_identity_16s$Day %in% settings$exvivo_days
microbiome_datasets <- list(mouse_genus = prepare_genus_context("mouse", mouse_keep), culture_genus = prepare_genus_context("culture", 
    culture_keep))
write_csv(sample_identity_16s, ppi_output(file.path(prep_dir, "phyloseq_sample_identity_audit.csv"), 
    name = "phyloseq_sample_identity_audit.csv"))
context_count_audit <- sample_identity_16s %>% count(Biological_Context, name = "n_phyloseq_rows")
write_csv(context_count_audit, ppi_output(file.path(prep_dir, "phyloseq_context_counts.csv"), name = "phyloseq_context_counts.csv"))
write_csv(bind_rows(microbiome_datasets$mouse_genus$audit, microbiome_datasets$culture_genus$audit), 
    ppi_output(file.path(prep_dir, "genus_filtering_audit.csv"), name = "genus_filtering_audit.csv"))
load_serum_lcms <- function(path) {
    obj <- readRDS(ppi_input(path))
    if (is.null(obj$serum_clean_bio) || is.null(obj$metadata_clean_bio)) {
        stop("Serum input object must contain serum_clean_bio and ", "metadata_clean_bio.")
    }
    X <- as.matrix(obj$serum_clean_bio)
    meta <- as.data.frame(obj$metadata_clean_bio, check.names = FALSE)
    required_columns <- c("Sample_Name", "Separation", "Donor")
    missing_columns <- setdiff(required_columns, names(meta))
    if (length(missing_columns) > 0L) {
        stop("Serum metadata is missing: ", paste(missing_columns, collapse = ", "))
    }
    aligned <- align_matrix_metadata(X, meta, "Sample_Name")
    serum_meta <- aligned$metadata
    sample_name_lower <- stringr::str_to_lower(stringr::str_squish(as.character(serum_meta$Sample_Name)))
    donor <- normalize_donor(serum_meta$Donor)
    donor_from_name <- normalize_donor(sample_name_lower)
    donor <- dplyr::coalesce(donor, donor_from_name)
    mouse_match <- stringr::str_match(sample_name_lower, "(?:ctrl|ppi)[_-]?[0-9]+[_-]([0-9]+)")
    mouse_number <- suppressWarnings(as.integer(mouse_match[, 2]))
    group <- factor(normalize_group(serum_meta$Separation), levels = c("ctrl", "ppi"))
    mouse_key <- dplyr::if_else(is.na(donor) | is.na(mouse_number), NA_character_, paste0(donor, "_mouse", 
        mouse_number))
    parsing_audit <- tibble::tibble(Sample_Name = serum_meta$Sample_Name, Group = group, Donor = donor, 
        Mouse_Number = mouse_number, Mouse_Key = mouse_key)
    if (anyNA(parsing_audit$Group)) {
        stop("Failed to identify the PPI/control group for these serum samples: ", paste(parsing_audit$Sample_Name[is.na(parsing_audit$Group)], 
            collapse = "; "))
    }
    if (anyNA(parsing_audit$Donor)) {
        stop("Failed to identify the inoculum donor for these serum samples: ", paste(parsing_audit$Sample_Name[is.na(parsing_audit$Donor)], 
            collapse = "; "))
    }
    if (anyNA(parsing_audit$Mouse_Number)) {
        stop("Failed to extract the recipient-mouse number from these serum samples: ", paste(parsing_audit$Sample_Name[is.na(parsing_audit$Mouse_Number)], 
            collapse = "; "))
    }
    if (anyDuplicated(parsing_audit$Mouse_Key)) {
        duplicated_keys <- unique(parsing_audit$Mouse_Key[duplicated(parsing_audit$Mouse_Key) | duplicated(parsing_audit$Mouse_Key, 
            fromLast = TRUE)])
        stop("Duplicated serum mouse keys were generated: ", paste(duplicated_keys, collapse = "; "))
    }
    meta_out <- serum_meta %>% dplyr::transmute(Sample_ID = .data$Sample_ID, Original_Sample_Name = .data$Sample_Name) %>% 
        dplyr::bind_cols(parsing_audit %>% dplyr::select(Group, Donor, Mouse_Number, Mouse_Key)) %>% 
        dplyr::mutate(Dataset = "serum_lcms")
    prep <- log2_half_min_matrix(aligned$matrix, "serum_lcms")
    list(matrix = prep$matrix, metadata = meta_out, audit = prep$audit, key_col = "Mouse_Key", context = "mouse")
}
load_cecum_feces_lcms <- function(path) {
    obj <- readRDS(ppi_input(path))
    if (is.null(obj$one_per_metabolite) || is.null(obj$metadata)) {
        stop("Caecum/faeces LC-MS object must contain one_per_metabolite and metadata.")
    }
    out <- list()
    for (compartment in c("cecum", "feces")) {
        X <- as.matrix(obj$one_per_metabolite[[compartment]])
        meta <- as.data.frame(obj$metadata[[compartment]])
        aligned <- align_matrix_metadata(X, meta, "Sample_Name")
        meta_out <- aligned$metadata %>% transmute(Sample_ID, Original_Sample_Name = Sample_Name, Group = factor(normalize_group(Separation), 
            levels = c("ctrl", "ppi")), Donor = normalize_donor(Donor), Mouse_Number = as.integer(Mouse_Number), 
            Mouse_Key = paste0(normalize_donor(Donor), "_mouse", as.integer(Mouse_Number)), Dataset = paste0(compartment, 
                "_lcms"))
        prep <- log2_half_min_matrix(aligned$matrix, paste0(compartment, "_lcms"))
        out[[paste0(compartment, "_lcms")]] <- list(matrix = prep$matrix, metadata = meta_out, audit = prep$audit, 
            key_col = "Mouse_Key", context = "mouse")
    }
    out
}
load_cecum_nmr <- function(path) {
    long <- read_excel(ppi_input(path), sheet = "Mouse_level_values") %>% mutate(sample_id = as.character(sample_id), 
        analysis_group = as.character(analysis_group), score = as.numeric(score))
    meta <- long %>% distinct(sample_id, inoculum_donor, group, mouse) %>% transmute(Sample_ID = sample_id, 
        Original_Sample_Name = sample_id, Group = factor(normalize_group(group), levels = c("ctrl", "ppi")), 
        Donor = normalize_donor(inoculum_donor), Mouse_Number = as.integer(mouse), Mouse_Key = paste0(normalize_donor(inoculum_donor), 
            "_mouse", as.integer(mouse)), Dataset = "cecum_nmr")
    wide <- long %>% select(sample_id, analysis_group, score) %>% pivot_wider(names_from = analysis_group, 
        values_from = score)
    X <- wide %>% select(-sample_id) %>% as.matrix()
    rownames(X) <- wide$sample_id
    meta <- meta[match(rownames(X), meta$Sample_ID), , drop = FALSE]
    prep <- prepare_score_matrix(X, "cecum_nmr")
    list(matrix = prep$matrix, metadata = meta, audit = prep$audit, key_col = "Mouse_Key", context = "mouse")
}
load_exvivo_nmr <- function(raw_path, mapping_path) {
    raw <- read_csv(ppi_input(raw_path), show_col_types = FALSE, name_repair = "unique")
    mapping <- read_csv(ppi_input(mapping_path), show_col_types = FALSE) %>% mutate(feature_id = as.character(feature_id), 
        analysis_group = as.character(analysis_group), include_group_level = as.logical(include_group_level)) %>% 
        filter(include_group_level, !is.na(analysis_group), analysis_group != "")
    classy_col <- first_existing_column(raw, c("CLASSYID", "feature_id"))
    sample_columns <- names(raw)[str_detect(names(raw), "^[1-4]_[CP][1-3]$")]
    info <- tibble(Sample_ID = sample_columns) %>% extract(Sample_ID, into = c("Day", "Short_Donor"), 
        regex = "^([1-4])_([CP][1-3])$", remove = FALSE, convert = TRUE) %>% mutate(Donor = if_else(str_starts(Short_Donor, 
        "C"), paste0("ctrl", str_extract(Short_Donor, "[0-9]+")), paste0("ppi", str_extract(Short_Donor, 
        "[0-9]+"))), Group = factor(if_else(str_starts(Short_Donor, "C"), "ctrl", "ppi"), levels = c("ctrl", 
        "ppi")), DayF = factor(Day, levels = settings$exvivo_days), Culture_Key = paste0(Donor, "_d", 
        Day), Dataset = "exvivo_nmr")
    feature_table <- raw %>% mutate(feature_id = as.character(.data[[classy_col]])) %>% inner_join(mapping %>% 
        select(feature_id, analysis_group), by = "feature_id")
    X_raw <- feature_table %>% select(all_of(sample_columns)) %>% as.matrix()
    storage.mode(X_raw) <- "numeric"
    rownames(X_raw) <- feature_table$feature_id
    X_raw[!is.finite(X_raw)] <- NA_real_
    X_raw <- pmax(X_raw, 0)
    pc <- apply(X_raw, 1, calculate_half_minimum)
    keep_feature <- is.finite(pc) & pc > 0
    X_raw <- X_raw[keep_feature, , drop = FALSE]
    pc <- pc[keep_feature]
    feature_table <- feature_table[keep_feature, , drop = FALSE]
    X_log <- X_raw
    for (i in seq_len(nrow(X_log))) X_log[i, ] <- log2(X_raw[i, ] + pc[[i]])
    X_z <- t(scale(t(X_log)))
    X_z[!is.finite(X_z)] <- 0
    groups <- split(feature_table$feature_id, feature_table$analysis_group)
    group_matrix <- vapply(groups, function(ids) {
        ids <- intersect(ids, rownames(X_z))
        if (length(ids) == 1L) 
            as.numeric(X_z[ids, ])
        else apply(X_z[ids, , drop = FALSE], 2, median, na.rm = TRUE)
    }, numeric(ncol(X_z)))
    X <- as.matrix(group_matrix)
    rownames(X) <- sample_columns
    prep <- prepare_score_matrix(X, "exvivo_nmr")
    list(matrix = prep$matrix, metadata = info %>% rename(Original_Sample_Name = Sample_ID), audit = prep$audit, 
        key_col = "Culture_Key", context = "culture")
}
load_exvivo_hilic <- function(path) {
    obj <- readRDS(ppi_input(path))
    if (is.null(obj$primary_matrix_all_samples) || is.null(obj$metadata_all)) {
        stop("Ex vivo HILIC object must contain primary_matrix_all_samples and metadata_all.")
    }
    X <- as.matrix(obj$primary_matrix_all_samples)
    meta <- as.data.frame(obj$metadata_all)
    required_meta <- c("Sample_Name", "is_analysis_sample", "is_medium", "is_day0", "Day", "Donor", "Separation")
    missing_meta <- setdiff(required_meta, names(meta))
    if (length(missing_meta) > 0L) {
        stop("HILIC metadata is missing: ", paste(missing_meta, collapse = ", "))
    }
    aligned <- align_matrix_metadata(X, meta, "Sample_Name")
    keep <- dplyr::coalesce(as.logical(aligned$metadata$is_analysis_sample), FALSE) & !dplyr::coalesce(as.logical(aligned$metadata$is_medium), 
        FALSE) & !dplyr::coalesce(as.logical(aligned$metadata$is_day0), FALSE)
    X <- aligned$matrix[keep, , drop = FALSE]
    meta <- aligned$metadata[keep, , drop = FALSE]
    day <- as.integer(str_extract(as.character(meta$Day), "[0-9]+"))
    donor <- normalize_donor(meta$Donor)
    meta_out <- tibble(Sample_ID = rownames(X), Original_Sample_Name = as.character(meta$Sample_Name), 
        Group = factor(normalize_group(meta$Separation), levels = c("ctrl", "ppi")), Donor = donor, Day = day, 
        DayF = factor(day, levels = settings$exvivo_days), Culture_Key = paste0(donor, "_d", day), Dataset = "exvivo_hilic")
    prep <- log2_half_min_matrix(X, "exvivo_hilic")
    list(matrix = prep$matrix, metadata = meta_out, audit = prep$audit, key_col = "Culture_Key", context = "culture")
}
metabolomics_datasets <- list()
metabolomics_datasets$serum_lcms <- load_serum_lcms(resolved_paths$serum_lcms)
metabolomics_datasets <- c(metabolomics_datasets, load_cecum_feces_lcms(resolved_paths$cecum_feces_lcms))
metabolomics_datasets$cecum_nmr <- load_cecum_nmr(resolved_paths$cecum_nmr)
metabolomics_datasets$exvivo_nmr <- load_exvivo_nmr(resolved_paths$exvivo_nmr_raw, resolved_paths$exvivo_nmr_mapping)
metabolomics_datasets$exvivo_hilic <- load_exvivo_hilic(resolved_paths$exvivo_hilic)
metabolite_audit <- bind_rows(map(metabolomics_datasets, "audit"))
write_csv(metabolite_audit, ppi_output(file.path(prep_dir, "metabolite_usability_audit.csv"), name = "metabolite_usability_audit.csv"))
make_matching_row <- function(layer, dataset, target, source_keys, target_keys) {
    source_keys <- sort(unique(na.omit(as.character(source_keys))))
    target_keys <- sort(unique(na.omit(as.character(target_keys))))
    matched <- intersect(source_keys, target_keys)
    tibble(layer = layer, dataset = dataset, target = target, n_source_keys = length(source_keys), n_target_keys = length(target_keys), 
        n_matched_keys = length(matched), matched_keys = paste(matched, collapse = ";"), missing_from_target = paste(setdiff(source_keys, 
            target_keys), collapse = ";"), extra_in_target = paste(setdiff(target_keys, source_keys), 
            collapse = ";"))
}
matching_audit <- list()
for (dataset_id in mouse_metabolomics_datasets) {
    ds <- metabolomics_datasets[[dataset_id]]
    met_keys <- unique(na.omit(ds$metadata$Mouse_Key))
    genus_keys <- unique(na.omit(microbiome_datasets$mouse_genus$metadata$Mouse_Key))
    matching_audit[[paste(dataset_id, "mouse_genus", sep = "__")]] <- make_matching_row("metabolite_genus", 
        dataset_id, "mouse_genus", met_keys, genus_keys)
    for (path_id in names(pathway_datasets)) {
        path_keys <- unique(na.omit(pathway_datasets[[path_id]]$metadata$Mouse_Key))
        matching_audit[[paste(dataset_id, path_id, sep = "__")]] <- make_matching_row("metabolite_pathway", 
            dataset_id, path_id, met_keys, path_keys)
    }
}
for (dataset_id in exvivo_metabolomics_datasets) {
    ds <- metabolomics_datasets[[dataset_id]]
    met_keys <- unique(na.omit(ds$metadata$Culture_Key))
    genus_keys <- unique(na.omit(microbiome_datasets$culture_genus$metadata$Culture_Key))
    matching_audit[[paste(dataset_id, "culture_genus", sep = "__")]] <- make_matching_row("metabolite_genus", 
        dataset_id, "culture_genus", met_keys, genus_keys)
}
mouse_genus_keys <- unique(na.omit(microbiome_datasets$mouse_genus$metadata$Mouse_Key))
for (path_id in names(pathway_datasets)) {
    path_keys <- unique(na.omit(pathway_datasets[[path_id]]$metadata$Mouse_Key))
    matching_audit[[paste("mouse_genus", path_id, sep = "__")]] <- make_matching_row("genus_pathway", 
        "mouse_genus", path_id, mouse_genus_keys, path_keys)
}
matching_audit <- bind_rows(matching_audit)
write_csv(matching_audit, ppi_output(file.path(prep_dir, "sample_matching_audit.csv"), name = "sample_matching_audit.csv"))
dataset_summary <- bind_rows(imap_dfr(metabolomics_datasets, function(ds, id) {
    tibble(data_type = "metabolomics", dataset_id = id, context = ds$context, n_samples = nrow(ds$matrix), 
        n_variables = ncol(ds$matrix), n_donors = n_distinct(ds$metadata$Donor))
}), imap_dfr(microbiome_datasets, function(ds, id) {
    tibble(data_type = "microbiome", dataset_id = id, context = if_else(id == "mouse_genus", "mouse", 
        "culture"), n_samples = nrow(ds$matrix), n_variables = ncol(ds$matrix), n_donors = n_distinct(ds$metadata$Donor))
}), imap_dfr(pathway_datasets, function(ds, id) {
    tibble(data_type = "pathway", dataset_id = id, context = "mouse", n_samples = nrow(ds$matrix), n_variables = ncol(ds$matrix), 
        n_donors = n_distinct(ds$metadata$Donor))
}))
write_csv(dataset_summary, ppi_output(file.path(prep_dir, "prepared_dataset_summary.csv"), name = "prepared_dataset_summary.csv"))
validation_issues <- character()
n_mouse_16s <- nrow(microbiome_datasets$mouse_genus$matrix)
n_culture_16s <- nrow(microbiome_datasets$culture_genus$matrix)
n_mouse_raw <- sum(sample_identity_16s$Biological_Context == "mouse", na.rm = TRUE)
n_culture_raw <- sum(sample_identity_16s$Biological_Context == "culture", na.rm = TRUE)
n_donor_raw <- sum(sample_identity_16s$Biological_Context == "donor", na.rm = TRUE)
n_unmatched_raw <- sum(sample_identity_16s$Biological_Context == "unmatched", na.rm = TRUE)
if (n_mouse_raw != settings$expected_mouse_16s_samples) {
    validation_issues <- c(validation_issues, paste0("phyloseq contains ", n_mouse_raw, " mouse rows; expected ", 
        settings$expected_mouse_16s_samples))
}
if (n_culture_raw != settings$expected_culture_16s_samples_all_days) {
    validation_issues <- c(validation_issues, paste0("phyloseq contains ", n_culture_raw, " culture rows across days 0-4; expected ", 
        settings$expected_culture_16s_samples_all_days))
}
if (n_donor_raw != settings$expected_original_donor_16s_samples) {
    validation_issues <- c(validation_issues, paste0("phyloseq contains ", n_donor_raw, " original donor rows; expected ", 
        settings$expected_original_donor_16s_samples))
}
if (n_unmatched_raw != settings$expected_unmatched_16s_samples) {
    validation_issues <- c(validation_issues, paste0("phyloseq contains ", n_unmatched_raw, " unmatched rows; expected ", 
        settings$expected_unmatched_16s_samples))
}
if (n_mouse_16s != settings$expected_mouse_16s_samples) {
    validation_issues <- c(validation_issues, paste0("mouse_genus contains ", n_mouse_16s, " biological samples; expected ", 
        settings$expected_mouse_16s_samples))
}
if (n_culture_16s != settings$expected_culture_16s_samples_days1_4) {
    validation_issues <- c(validation_issues, paste0("culture_genus contains ", n_culture_16s, " donor-day samples after excluding day 0; expected ", 
        settings$expected_culture_16s_samples_days1_4))
}
mouse_16s_keys <- sort(unique(na.omit(microbiome_datasets$mouse_genus$metadata$Mouse_Key)))
for (dataset_id in mouse_metabolomics_datasets) {
    dataset_keys <- sort(unique(na.omit(metabolomics_datasets[[dataset_id]]$metadata$Mouse_Key)))
    if (!setequal(mouse_16s_keys, dataset_keys)) {
        validation_issues <- c(validation_issues, paste0(dataset_id, " keys do not match mouse 16S keys. Missing from 16S: ", 
            paste(setdiff(dataset_keys, mouse_16s_keys), collapse = ";"), ". Extra in 16S: ", paste(setdiff(mouse_16s_keys, 
                dataset_keys), collapse = ";")))
    }
    if (anyDuplicated(metabolomics_datasets[[dataset_id]]$metadata$Mouse_Key)) {
        validation_issues <- c(validation_issues, paste0(dataset_id, " contains duplicated Mouse_Key values"))
    }
}
culture_16s_keys <- sort(unique(na.omit(microbiome_datasets$culture_genus$metadata$Culture_Key)))
for (dataset_id in exvivo_metabolomics_datasets) {
    dataset_keys <- sort(unique(na.omit(metabolomics_datasets[[dataset_id]]$metadata$Culture_Key)))
    if (!setequal(culture_16s_keys, dataset_keys)) {
        validation_issues <- c(validation_issues, paste0(dataset_id, " keys do not match culture 16S keys. Missing from 16S: ", 
            paste(setdiff(dataset_keys, culture_16s_keys), collapse = ";"), ". Extra in 16S: ", paste(setdiff(culture_16s_keys, 
                dataset_keys), collapse = ";")))
    }
    if (anyDuplicated(metabolomics_datasets[[dataset_id]]$metadata$Culture_Key)) {
        validation_issues <- c(validation_issues, paste0(dataset_id, " contains duplicated Culture_Key values"))
    }
}
expected_rna_by_organ <- rna_bridge %>% filter(!is.na(Mouse_Key), !is.na(Organ)) %>% distinct(Organ, 
    Mouse_Key) %>% count(Organ, name = "expected_n")
for (path_id in names(pathway_datasets)) {
    ds <- pathway_datasets[[path_id]]
    expected_n <- expected_rna_by_organ$expected_n[match(ds$organ, expected_rna_by_organ$Organ)]
    actual_n <- nrow(ds$matrix)
    if (length(expected_n) != 1L || is.na(expected_n) || actual_n != expected_n) {
        validation_issues <- c(validation_issues, paste0(path_id, " contains ", actual_n, " samples; expected ", 
            ifelse(length(expected_n) == 1L, expected_n, NA_integer_), " for organ ", ds$organ))
    }
    if (anyDuplicated(ds$metadata$Mouse_Key)) {
        validation_issues <- c(validation_issues, paste0(path_id, " contains duplicated Mouse_Key values"))
    }
    if (any(as.character(ds$metadata$Organ) != ds$organ)) {
        validation_issues <- c(validation_issues, paste0(path_id, " contains samples from another organ"))
    }
}
validation_table <- tibble(status = if (length(validation_issues) == 0L) "PASS" else "FAIL", issue = if (length(validation_issues) == 
    0L) "All structural checks passed" else validation_issues)
write_csv(validation_table, ppi_output(file.path(prep_dir, "structural_validation.csv"), name = "structural_validation.csv"))
writexl::write_xlsx(list(Structural_validation = validation_table, Dataset_summary = dataset_summary, 
    Matching_audit = matching_audit, Pathway_selection = pathway_selection, RNA_bridge = rna_bridge, 
    Phyloseq_context_counts = context_count_audit, Genus_filtering = bind_rows(microbiome_datasets$mouse_genus$audit, 
        microbiome_datasets$culture_genus$audit), Metabolite_filtering = metabolite_audit, Resolved_paths = path_audit), 
    path = ppi_output(file.path(prep_dir, "PPI_multiomics_input_audit.xlsx"), name = "ppi_multiomics_input_audit.xlsx"))
if (length(validation_issues) > 0L && isTRUE(settings$fail_on_audit_mismatch)) {
    stop("Input preparation failed structural validation:\n- ", paste(validation_issues, collapse = "\n- "), 
        "\n\nThe audit CSV/XLSX files were written, but no prepared RDS bundle was saved. ", "Correct the inputs or parsing and rerun script 01.")
}
prepared_bundle <- list(metabolomics = metabolomics_datasets, microbiome = microbiome_datasets, pathways = pathway_datasets, 
    pathway_selection = pathway_selection_present, rna_bridge = rna_bridge, phyloseq_sample_identity = sample_identity_16s, 
    phyloseq_context_counts = context_count_audit, matching_audit = matching_audit, dataset_summary = dataset_summary, 
    resolved_paths = path_audit, settings = settings, run_layers = run_layers, mouse_metabolomics_datasets = mouse_metabolomics_datasets, 
    exvivo_metabolomics_datasets = exvivo_metabolomics_datasets)
saveRDS(prepared_bundle, ppi_output(prepared_output_file), compress = "xz")
message("Preparation complete. Outputs written to: ", prep_dir)
