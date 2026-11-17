# SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0
# Ex vivo NMR

source("scripts/helpers/setup.R")
invisible(run_analysis(
  section = "04_mouse_metabolomics",
  directory = "metabolomics",
  script = "04_07_analyse_ex_vivo_nmr.R",
  restricted = FALSE,
  code = {
required_packages <- c("dplyr", "tidyr", "stringr", "tibble", "readr", "limma", "statmod", "writexl")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]
if (length(missing_packages) > 0) {
    stop("Install these packages first: ", paste(missing_packages, collapse = ", "))
}
raw_csv <- project_path("metabolomics/NMR_data/NMR_NEW_raw_data.csv")
out_dir <- project_path("metabolomics/NMR_data/NMR_annotated_limma_dupCor")
dir.create(project_output(out_dir, directory = TRUE), recursive = TRUE, showWarnings = FALSE)
keep_days <- c("1", "2", "3", "4")
group_levels <- c("ctrl", "ppi")
clean_metabolite_name <- function(x) {
    x <- stringr::str_squish(as.character(x))
    x <- dplyr::na_if(x, "")
    xl <- stringr::str_to_lower(x)
    dplyr::case_when(is.na(x) ~ NA_character_, xl %in% c("lactaate", "lactate2") ~ "lactate", xl == "his" ~
        "histidine", TRUE ~ xl)
}
assign_analysis_group <- function(feature_id, canonical_metabolite) {
    dplyr::case_when(canonical_metabolite == "tsp" ~ "reference__TSP", feature_id %in% c("CLASSY611",
        "CLASSY612") ~ "pyroglutamate_4.18", feature_id %in% c("CLASSY253", "CLASSY256") ~ "pyroglutamate_1.78",
        feature_id %in% c("CLASSY267", "CLASSY268", "CLASSY270") ~ "pyroglutamate_1.86", feature_id ==
            "CLASSY355" ~ "succinate_2.41", feature_id == "CLASSY254" ~ "ambiguous_1.78_overlap", TRUE ~
            canonical_metabolite)
}
assign_group_label <- function(analysis_group) {
    dplyr::case_when(analysis_group == "pyroglutamate_4.18" ~ "pyroglutamate, 4.18-ppm resonance", analysis_group ==
        "pyroglutamate_1.78" ~ "pyroglutamate, 1.78-ppm resonance", analysis_group == "pyroglutamate_1.86" ~
        "pyroglutamate, 1.86-ppm resonance", analysis_group == "succinate_2.41" ~ "succinate, 2.41-ppm resonance",
        analysis_group == "ambiguous_1.78_overlap" ~ "ambiguous 1.78-ppm overlap", analysis_group ==
            "reference__TSP" ~ "TSP", TRUE ~ analysis_group)
}
assign_group_type <- function(analysis_group) {
    dplyr::case_when(stringr::str_starts(analysis_group, "pyroglutamate_") ~ "pyroglutamate_resonance_group",
        analysis_group == "succinate_2.41" ~ "metabolite_resonance", analysis_group == "ambiguous_1.78_overlap" ~
            "ambiguous_overlap", analysis_group == "reference__TSP" ~ "reference_compound", TRUE ~ "metabolite")
}
assign_confidence <- function(analysis_group) {
    dplyr::case_when(analysis_group == "pyroglutamate_4.18" ~ "high", analysis_group == "pyroglutamate_1.78" ~
        "moderate", analysis_group == "pyroglutamate_1.86" ~ "low", analysis_group == "succinate_2.41" ~
        "moderate", analysis_group == "ambiguous_1.78_overlap" ~ "low", analysis_group == "reference__TSP" ~
        "not applicable", TRUE ~ "standard")
}
assign_interpretation_note <- function(analysis_group) {
    dplyr::case_when(analysis_group == "pyroglutamate_4.18" ~ paste("Highly concordant pair and the strongest-supported",
        "pyroglutamate resonance group."), analysis_group == "pyroglutamate_1.78" ~ paste("Highly concordant local pair, but located in an",
        "overlap-prone spectral region."), analysis_group == "pyroglutamate_1.86" ~ paste("Candidate local multiplet; member features show only",
        "moderate concordance. Interpret cautiously."), analysis_group == "succinate_2.41" ~ paste("Retained as the primary succinate resonance and not",
        "combined with CLASSY254."), analysis_group == "ambiguous_1.78_overlap" ~ paste("Excluded from primary group-level analysis. Labelled",
        "succinate in the export but highly concordant with the", "local pyroglutamate 1.78-ppm envelope."),
        analysis_group == "reference__TSP" ~ "Excluded: NMR reference/calibration compound.", TRUE ~
            "")
}
calculate_feature_pseudocount <- function(x) {
    positive <- x[is.finite(x) & x > 0]
    if (length(positive) == 0) {
        return(NA_real_)
    }
    min(positive, na.rm = TRUE)/2
}
calculate_group_concordance <- function(group_name, member_ids, X_log2) {
    member_ids <- intersect(member_ids, rownames(X_log2))
    if (length(member_ids) == 1) {
        return(tibble::tibble(analysis_group = group_name, n_CLASSY_features = 1L, member_CLASSY_IDs = member_ids,
            median_pairwise_spearman = NA_real_, minimum_pairwise_spearman = NA_real_, maximum_pairwise_spearman = NA_real_,
            low_concordance_flag = FALSE))
    }
    R <- stats::cor(t(X_log2[member_ids, , drop = FALSE]), method = "spearman", use = "pairwise.complete.obs")
    pairwise_rho <- R[upper.tri(R)]
    tibble::tibble(analysis_group = group_name, n_CLASSY_features = length(member_ids), member_CLASSY_IDs = paste(member_ids,
        collapse = "; "), median_pairwise_spearman = stats::median(pairwise_rho, na.rm = TRUE), minimum_pairwise_spearman = min(pairwise_rho,
        na.rm = TRUE), maximum_pairwise_spearman = max(pairwise_rho, na.rm = TRUE), low_concordance_flag = stats::median(pairwise_rho,
        na.rm = TRUE) < 0.5 | min(pairwise_rho, na.rm = TRUE) < 0)
}
make_group_score <- function(member_ids, X_z) {
    member_ids <- intersect(member_ids, rownames(X_z))
    if (length(member_ids) == 1) {
        return(as.numeric(X_z[member_ids, , drop = TRUE]))
    }
    apply(X_z[member_ids, , drop = FALSE], 2, stats::median, na.rm = TRUE)
}
nmr_raw <- readr::read_csv(project_input(raw_csv), show_col_types = FALSE)
required_columns <- c("metabolite", "CLASSYID", "ppm")
missing_columns <- setdiff(required_columns, names(nmr_raw))
if (length(missing_columns) > 0) {
    stop("Missing required columns: ", paste(missing_columns, collapse = ", "))
}
all_sample_cols <- names(nmr_raw)[stringr::str_detect(names(nmr_raw), "^[0-9]+_[CP][0-9]+$")]
sample_info <- tibble::tibble(sample_col = all_sample_cols) %>% tidyr::extract(sample_col, into = c("day",
    "donor"), regex = "^([0-9]+)_([CP][0-9]+)$", remove = FALSE) %>% dplyr::mutate(group = dplyr::case_when(stringr::str_starts(donor,
    "C") ~ "ctrl", stringr::str_starts(donor, "P") ~ "ppi", TRUE ~ NA_character_), grp = factor(group,
    levels = group_levels), dayF = factor(day, levels = keep_days), donor = factor(donor)) %>% dplyr::filter(day %in%
    keep_days)
if (nrow(sample_info) != 24) {
    warning("Expected 24 biological samples from days 1-4, but found ", nrow(sample_info), ". Check the sample-column parsing.")
}
feature_info <- nmr_raw %>% dplyr::transmute(feature_id = as.character(CLASSYID), original_metabolite = stringr::str_squish(as.character(metabolite)),
    canonical_metabolite = clean_metabolite_name(metabolite), ppm = as.numeric(ppm)) %>% dplyr::filter(!is.na(canonical_metabolite)) %>%
    dplyr::mutate(analysis_group = assign_analysis_group(feature_id, canonical_metabolite), group_label = assign_group_label(analysis_group),
        group_type = assign_group_type(analysis_group), confidence = assign_confidence(analysis_group),
        interpretation_note = assign_interpretation_note(analysis_group), include_feature_level = analysis_group !=
            "reference__TSP", include_group_level = !analysis_group %in% c("reference__TSP", "ambiguous_1.78_overlap"))
if (anyDuplicated(feature_info$feature_id)) {
    stop("CLASSY IDs are not unique.")
}
expected_manual_features <- c("CLASSY611", "CLASSY612", "CLASSY253", "CLASSY256", "CLASSY267", "CLASSY268",
    "CLASSY270", "CLASSY355", "CLASSY254")
missing_manual_features <- setdiff(expected_manual_features, feature_info$feature_id)
if (length(missing_manual_features) > 0) {
    stop("Expected manually curated CLASSY features are missing: ", paste(missing_manual_features, collapse = ", "))
}
feature_level_info <- feature_info %>% dplyr::filter(include_feature_level)
X_raw <- nmr_raw %>% dplyr::filter(CLASSYID %in% feature_level_info$feature_id) %>% dplyr::arrange(match(CLASSYID,
    feature_level_info$feature_id)) %>% dplyr::select(dplyr::all_of(sample_info$sample_col)) %>% as.matrix()
rownames(X_raw) <- feature_level_info$feature_id
storage.mode(X_raw) <- "numeric"
negative_summary <- tibble::tibble(feature_id = rownames(X_raw), n_negative = rowSums(X_raw < 0, na.rm = TRUE),
    proportion_negative = rowMeans(X_raw < 0, na.rm = TRUE)) %>% dplyr::left_join(feature_level_info,
    by = "feature_id")
X_nonnegative <- pmax(X_raw, 0)
feature_pseudocount <- apply(X_nonnegative, 1, calculate_feature_pseudocount)
valid_feature <- is.finite(feature_pseudocount)
if (any(!valid_feature)) {
    warning(sum(!valid_feature), " annotated CLASSY features had no positive values and were removed.")
}
X_nonnegative <- X_nonnegative[valid_feature, , drop = FALSE]
feature_pseudocount <- feature_pseudocount[valid_feature]
feature_level_info <- feature_level_info %>% dplyr::filter(feature_id %in% rownames(X_nonnegative))
X_log2 <- X_nonnegative
for (i in seq_len(nrow(X_log2))) {
    X_log2[i, ] <- log2(X_nonnegative[i, ] + feature_pseudocount[i])
}
group_level_info <- feature_level_info %>% dplyr::filter(include_group_level)
group_members <- split(group_level_info$feature_id, group_level_info$analysis_group)
X_z <- t(scale(t(X_log2)))
Y_group <- vapply(group_members, make_group_score, X_z = X_z, FUN.VALUE = numeric(ncol(X_z)))
Y_group <- t(Y_group)
colnames(Y_group) <- colnames(X_z)
group_variance <- apply(Y_group, 1, stats::var, na.rm = TRUE)
Y_group <- Y_group[is.finite(group_variance) & group_variance > 0, , drop = FALSE]
group_concordance <- dplyr::bind_rows(lapply(names(group_members), function(group_name) {
    calculate_group_concordance(group_name = group_name, member_ids = group_members[[group_name]], X_log2 = X_log2)
}))
group_summary <- group_level_info %>% dplyr::filter(analysis_group %in% rownames(Y_group)) %>% dplyr::group_by(analysis_group,
    group_label, group_type, confidence, interpretation_note) %>% dplyr::summarise(canonical_metabolite = paste(sort(unique(canonical_metabolite)),
    collapse = "; "), n_CLASSY_features = dplyr::n_distinct(feature_id), member_CLASSY_IDs = paste(feature_id,
    collapse = "; "), original_names = paste(sort(unique(original_metabolite)), collapse = "; "), ppm_values = paste(sort(unique(round(ppm,
    5))), collapse = "; "), ppm_min = min(ppm, na.rm = TRUE), ppm_max = max(ppm, na.rm = TRUE), .groups = "drop") %>%
    dplyr::left_join(group_concordance, by = c("analysis_group", "n_CLASSY_features", "member_CLASSY_IDs"))
design <- stats::model.matrix(~grp + dayF, data = sample_info)
colnames(design) <- make.names(colnames(design))
group_coefficient <- "grpppi"
if (!group_coefficient %in% colnames(design)) {
    stop("Could not find treatment coefficient ", group_coefficient, ". Design columns are: ", paste(colnames(design),
        collapse = ", "))
}
dupcor_group <- limma::duplicateCorrelation(Y_group, design = design, block = sample_info$donor)
group_consensus_correlation <- dupcor_group$consensus.correlation
if (!is.finite(group_consensus_correlation)) {
    warning("Group-level duplicateCorrelation was not finite; ", "using correlation = 0.")
    group_consensus_correlation <- 0
}
fit_group <- limma::lmFit(Y_group, design = design, block = sample_info$donor, correlation = group_consensus_correlation)
fit_group <- limma::eBayes(fit_group, robust = TRUE)
group_results <- limma::topTable(fit_group, coef = group_coefficient, number = Inf, sort.by = "P") %>%
    tibble::rownames_to_column("analysis_group") %>% dplyr::rename(score_difference = logFC, average_score = AveExpr,
    t_statistic = t, p_value = P.Value, FDR = adj.P.Val, log_odds = B)
coef_index <- match(group_coefficient, colnames(fit_group$coefficients))
moderated_SE <- fit_group$stdev.unscaled[, coef_index] * sqrt(fit_group$s2.post)
df_total <- fit_group$df.total
critical_value <- stats::qt(0.975, df = df_total)
confidence_table <- tibble::tibble(analysis_group = rownames(fit_group$coefficients), moderated_SE = moderated_SE,
    df_total = df_total, conf_low = fit_group$coefficients[, coef_index] - critical_value * moderated_SE,
    conf_high = fit_group$coefficients[, coef_index] + critical_value * moderated_SE)
group_results <- group_results %>% dplyr::left_join(confidence_table, by = "analysis_group") %>% dplyr::left_join(group_summary,
    by = "analysis_group") %>% dplyr::mutate(direction = dplyr::case_when(score_difference > 0 ~ "higher in PPI",
    score_difference < 0 ~ "lower in PPI", TRUE ~ "no difference"), effect_scale = "difference in standardized group score") %>%
    dplyr::arrange(p_value)
dupcor_feature <- limma::duplicateCorrelation(X_log2, design = design, block = sample_info$donor)
feature_consensus_correlation <- dupcor_feature$consensus.correlation
if (!is.finite(feature_consensus_correlation)) {
    warning("Feature-level duplicateCorrelation was not finite; ", "using correlation = 0.")
    feature_consensus_correlation <- 0
}
fit_feature <- limma::lmFit(X_log2, design = design, block = sample_info$donor, correlation = feature_consensus_correlation)
fit_feature <- limma::eBayes(fit_feature, robust = TRUE)
feature_results <- limma::topTable(fit_feature, coef = group_coefficient, number = Inf, sort.by = "P") %>%
    tibble::rownames_to_column("feature_id") %>% dplyr::rename(log2_difference = logFC, p_value = P.Value,
    FDR = adj.P.Val) %>% dplyr::left_join(feature_level_info, by = "feature_id") %>% dplyr::mutate(direction = dplyr::case_when(log2_difference >
    0 ~ "higher in PPI", log2_difference < 0 ~ "lower in PPI", TRUE ~ "no difference"), primary_group_level_feature = include_group_level) %>%
    dplyr::arrange(p_value)
manual_region_features <- c("CLASSY611", "CLASSY612", "CLASSY253", "CLASSY254", "CLASSY256", "CLASSY267",
    "CLASSY268", "CLASSY270", "CLASSY355")
manual_region_matrix <- X_log2[intersect(manual_region_features, rownames(X_log2)), , drop = FALSE]
manual_region_correlations <- stats::cor(t(manual_region_matrix), method = "spearman", use = "pairwise.complete.obs") %>%
    as.data.frame() %>% tibble::rownames_to_column("feature_id")
Y_group <- Y_group[, sample_info$sample_col, drop = FALSE]
donor_levels <- levels(sample_info$donor)
Y_donor <- vapply(donor_levels, function(d) {
    rowMeans(Y_group[, sample_info$donor == d, drop = FALSE], na.rm = TRUE)
}, FUN.VALUE = numeric(nrow(Y_group)))
rownames(Y_donor) <- rownames(Y_group)
colnames(Y_donor) <- donor_levels
donor_info <- sample_info %>% dplyr::distinct(donor, grp) %>% dplyr::arrange(match(donor, donor_levels))
design_donor <- model.matrix(~grp, data = donor_info)
colnames(design_donor) <- make.names(colnames(design_donor))
fit_donor <- limma::lmFit(Y_donor, design_donor)
fit_donor <- limma::eBayes(fit_donor, robust = TRUE)
donor_results <- limma::topTable(fit_donor, coef = "grpppi", number = Inf, sort.by = "P") %>% tibble::rownames_to_column("analysis_group")
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(stringr)
group_long <- as.data.frame(t(Y_group)) %>% tibble::rownames_to_column("sample_col") %>% dplyr::left_join(sample_info %>%
    dplyr::select(sample_col, day, dayF, donor, group, grp), by = "sample_col") %>% tidyr::pivot_longer(cols = -c(sample_col,
    day, dayF, donor, group, grp), names_to = "analysis_group", values_to = "score")
plot_metabolite_trajectory <- function(metabolite, data = group_long) {
    plot_df <- data %>% dplyr::filter(analysis_group == metabolite)
    if (nrow(plot_df) == 0) {
        stop("Metabolite not found: ", metabolite)
    }
    ggplot(plot_df, aes(x = dayF, y = score, group = donor, colour = grp)) + geom_hline(yintercept = 0,
        linetype = "dashed", linewidth = 0.4) + geom_line(alpha = 0.65, linewidth = 0.7) + geom_point(size = 2.8) +
        stat_summary(aes(group = grp), fun = mean, geom = "line", linewidth = 1.4) + stat_summary(aes(group = grp),
        fun = mean, geom = "point", size = 3.8, shape = 18) + scale_colour_discrete(labels = c(ctrl = "Control",
        ppi = "PPI")) + labs(title = stringr::str_replace_all(metabolite, "_", " "), subtitle = "Individual culture trajectories; diamonds show group means",
        x = "Culture day", y = "Standardized metabolite score", colour = NULL) + theme_bw() + theme(plot.title = element_text(face = "bold"),
        legend.position = "top")
}
plot_metabolite_trajectory("trimethylamine")
plot_metabolite_trajectory("tryptophan")
plot_metabolite_trajectory("isocaproate")
plot_metabolite_trajectory("pyroglutamate_4.18")
donor_average_long <- group_long %>% dplyr::group_by(analysis_group, donor, grp) %>% dplyr::summarise(donor_mean_score = mean(score,
    na.rm = TRUE), .groups = "drop")
plot_metabolite_donor_average <- function(metabolite, data = donor_average_long, results = donor_results) {
    plot_df <- data %>% dplyr::filter(analysis_group == metabolite)
    if (nrow(plot_df) == 0) {
        stop("Metabolite not found: ", metabolite)
    }
    result_row <- results %>% dplyr::filter(analysis_group == metabolite)
    subtitle_text <- NULL
    if (nrow(result_row) == 1) {
        subtitle_text <- paste0("Donor-averaged limma: p = ", format.pval(result_row$P.Value, digits = 3,
            eps = 0.001), "; FDR = ", format.pval(result_row$adj.P.Val, digits = 3, eps = 0.001))
    }
    ggplot(plot_df, aes(x = grp, y = donor_mean_score, colour = grp)) + geom_hline(yintercept = 0, linetype = "dashed",
        linewidth = 0.4) + geom_point(size = 3.5, position = position_jitter(width = 0.05, height = 0)) +
        stat_summary(fun = mean, geom = "crossbar", width = 0.45, linewidth = 0.7, colour = "black") +
        scale_x_discrete(labels = c(ctrl = "Control", ppi = "PPI")) + guides(colour = "none") + labs(title = stringr::str_replace_all(metabolite,
        "_", " "), subtitle = subtitle_text, x = NULL, y = "Mean standardized score across days 1–4") +
        theme_bw() + theme(plot.title = element_text(face = "bold"))
}
plot_metabolite_donor_average("trimethylamine")
plot_metabolite_donor_average("tryptophan")
plot_metabolite_donor_average("isocaproate")
metabolites_to_plot <- c("trimethylamine", "tryptophan", "isocaproate", "alpha-aminobutyrate", "pyroglutamate_4.18",
    "succinate_2.41")
group_long %>% dplyr::filter(analysis_group %in% metabolites_to_plot) %>% ggplot(aes(x = dayF, y = score,
    group = donor, colour = grp)) + geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3) +
    geom_line(alpha = 0.55) + geom_point(size = 1.8) + stat_summary(aes(group = grp), fun = mean, geom = "line",
    linewidth = 1.2) + facet_wrap(~analysis_group, scales = "free_y", labeller = labeller(analysis_group = function(x) {
    stringr::str_replace_all(x, "_", " ")
})) + scale_colour_discrete(labels = c(ctrl = "Control", ppi = "PPI")) + labs(x = "Culture day", y = "Standardized metabolite score",
    colour = NULL) + theme_bw() + theme(legend.position = "top", strip.text = element_text(face = "bold"))
excluded_or_ambiguous <- feature_info %>% dplyr::filter(!include_group_level)
analysis_summary <- tibble::tibble(item = c("Annotated CLASSY features retained for feature analysis",
    "Primary metabolite/resonance groups tested", "Samples", "Donors", "Days", "Group-level consensus correlation",
    "Feature-level consensus correlation", "Primary model", "Group-level effect scale", "Pyroglutamate treatment",
    "Succinate treatment"), value = c(nrow(feature_level_info), nrow(Y_group), ncol(Y_group), dplyr::n_distinct(sample_info$donor),
    paste(keep_days, collapse = ", "), group_consensus_correlation, feature_consensus_correlation, "group score ~ group + categorical day; block = donor",
    "standardized score difference", paste("Three separate resonance groups:", "4.18 ppm, 1.78 ppm and 1.86 ppm"),
    paste("CLASSY355 retained as succinate_2.41;", "CLASSY254 retained only in feature-level sensitivity",
        "as ambiguous_1.78_overlap")))
settings <- tibble::tibble(setting = c("raw_csv", "annotated_only", "negative_values", "transformation",
    "general_duplicate_combination", "lactate_cleanup", "pyroglutamate_cleanup", "succinate_cleanup",
    "multiple_testing"), value = c(raw_csv, "TRUE", "truncated to zero", paste("feature-specific log2 pseudocount,",
    "then z-score before group combination"), "median standardized intensity for same analysis_group",
    "lactaate + lactate2 -> lactate", paste("split into pyroglutamate_4.18,", "pyroglutamate_1.78 and pyroglutamate_1.86"),
    paste("CLASSY355 -> succinate_2.41;", "CLASSY254 -> ambiguous_1.78_overlap", "and excluded from primary group testing"),
    "BH FDR across primary metabolite/resonance groups"))
output_xlsx <- file.path(out_dir, "NMR_exvivo_annotated_limma_duplicateCorrelation.xlsx")
writexl::write_xlsx(list(Group_level_results = group_results, Feature_level_results = feature_results,
    Cleaned_feature_mapping = feature_info, Group_summary = group_summary, Group_concordance = group_concordance,
    Manual_region_correlations = manual_region_correlations, Excluded_or_ambiguous = excluded_or_ambiguous,
    Negative_value_summary = negative_summary, Sample_metadata = sample_info, Analysis_summary = analysis_summary,
    Settings = settings), project_output(output_xlsx))
readr::write_csv(group_results, project_output(file.path(out_dir, "NMR_exvivo_annotated_group_results_v2.csv"),
    name = "nmr_exvivo_annotated_group_results.csv"))
readr::write_csv(feature_results, project_output(file.path(out_dir, "NMR_exvivo_annotated_feature_results_v2.csv"),
    name = "nmr_exvivo_annotated_feature_results.csv"))
readr::write_csv(feature_info, project_output(file.path(out_dir, "NMR_exvivo_cleaned_feature_mapping_v2.csv"),
    name = "nmr_exvivo_cleaned_feature_mapping.csv"))
message("\nAnalysis complete.\n", "Main workbook:\n", output_xlsx, "\n\nGroup-level consensus correlation: ",
    round(group_consensus_correlation, 4), "\n")
print(group_results %>% dplyr::select(analysis_group, group_label, group_type, confidence, n_CLASSY_features,
    score_difference, conf_low, conf_high, p_value, FDR, direction, low_concordance_flag) %>% head(20))
}))
