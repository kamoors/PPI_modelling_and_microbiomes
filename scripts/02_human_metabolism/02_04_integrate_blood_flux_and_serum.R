# SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0
# WBM/serum representative tables

source("scripts/helpers/setup.R")
invisible(run_analysis(
  section = "02_human_metabolism",
  directory = "human_cohorts/blood_flux/metacardis",
  script = "02_04_integrate_blood_flux_and_serum.R",
  restricted = FALSE,
  code = {
required_packages <- c("dplyr", "tidyr", "readr", "stringr", "forcats", "ggplot2", "ggrepel", "patchwork",
    "scales")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0L) {
    stop("Install the following packages before running the script: ", paste(missing_packages, collapse = ", "))
}
results_dir <- project_path("human_cohorts/results")
blood_dir <- project_path("human_cohorts/blood_flux/metacardis/Blood_FBA_Metacardis/resultFlux_Metacardis_v2/analysis_outputs/Kiel_Metacardis_comparison")
wbm_file_default <- file.path(blood_dir, "Blood_FBA_Kiel_MetaCardis_v2_fixed_effect_meta_analysis.csv")
metabolomics_file_default <- file.path(results_dir, "Metacardis_metabolomics_results_fullFormula_preBoot_NoPPIexluded.csv")
output_dir_default <- file.path(results_dir, "MetaCardis_WBM_metabolomics_integrated_figure")
wbm_file <- project_path("human_cohorts/blood_flux/metacardis/Blood_FBA_Metacardis/resultFlux_Metacardis_v2/analysis_outputs/Kiel_Metacardis_comparison/Blood_FBA_Kiel_MetaCardis_v2_fixed_effect_meta_analysis.csv")
metabolomics_file <- project_path("human_cohorts/results/Metacardis_metabolomics_fullFormula_totalMedicationAdjusted_noPPIExcludedByMetabolomics_preBoot.csv")
output_dir <- output_dir_default
if (!project_dir_exists(results_dir)) {
    stop("The configured Kiel–MetaCardis results directory does not exist:\n", results_dir)
}
if (!project_file_exists(wbm_file)) {
    stop("The WBM meta-analysis file was not found:\n", wbm_file, "\nExpected default location:\n", wbm_file_default)
}
if (!project_file_exists(metabolomics_file)) {
    stop("The MetaCardis metabolomics results file was not found:\n", metabolomics_file, "\nExpected default location:\n",
        metabolomics_file_default, "\nIf the file is stored elsewhere, either change ", "metabolomics_file_default or pass its path as the second command-line argument.")
}
dir.create(project_output(output_dir, directory = TRUE), recursive = TRUE, showWarnings = FALSE)
project_setwd(results_dir)
manual_curation_file <- file.path(output_dir, "07_manual_curation_template_v2.csv")
fdr_cutoff <- 0.1
wbm_max_i2 <- 25
wbm_min_q_pvalue <- 0.05
max_wbm_representatives_per_theme <- 6L
max_measured_metabolites_per_theme <- 7L
collapse_wbm_near_identical_estimates <- TRUE
wbm_cluster_digits <- 7L
exclude_medications_from_main <- TRUE
exclude_unidentified_from_main <- TRUE
minimum_abs_effect_for_main_wbm <- 1e-06
figure_width_mm <- 270
figure_height_mm <- 270
figure_dpi <- 600
write_csv_safe <- function(x, filename) {
    readr::write_csv(x, project_output(file.path(output_dir, filename)), na = "")
}
assert_columns <- function(x, required, object_name) {
    missing <- setdiff(required, names(x))
    if (length(missing) > 0L) {
        stop(object_name, " is missing required columns: ", paste(missing, collapse = ", "))
    }
}
as_logical_safe <- function(x) {
    if (is.logical(x)) {
        return(x)
    }
    z <- tolower(trimws(as.character(x)))
    out <- rep(NA, length(z))
    out[z %in% c("true", "t", "1", "yes", "y")] <- TRUE
    out[z %in% c("false", "f", "0", "no", "n")] <- FALSE
    out
}
canonical_name <- function(x) {
    stringr::str_replace_all(stringr::str_replace_all(stringr::str_replace_all(stringr::str_replace_all(stringr::str_to_lower(as.character(x)),
        "\\([^)]*\\)", ""), "alpha", "a"), "gamma", "g"), "[^a-z0-9]+", "")
}
clean_display_name <- function(x) {
    stringr::str_replace(stringr::str_replace(stringr::str_replace(stringr::str_replace(stringr::str_replace(stringr::str_replace(stringr::str_replace(stringr::str_replace(stringr::str_replace(stringr::str_squish(stringr::str_replace_all(stringr::str_replace_all(stringr::str_replace_all(stringr::str_replace_all(as.character(x),
        "\\.", " "), "_", " "), "\\*", ""), "\\s*\\([^)]*\\)\\s*$", "")), "^Alpha-N-Phenylacetyl-L-glutamine$",
        "Phenylacetylglutamine"), "^3-Methyl-2-Oxovaleric Acid$", "3-Methyl-2-oxovalerate"), "^alpha-Ketoisovaleric acid$",
        "α-Ketoisovalerate"), "^4-Hydroxyphenylpyruvic acid$", "4-Hydroxyphenylpyruvate"), "^Imidazoleacetic acid$",
        "Imidazoleacetate"), "^Methylimidazoleacetic acid$", "Methylimidazoleacetate"), "^3,4-Dihydroxybenzeneacetic acid$",
        "3,4-Dihydroxyphenylacetate"), "^Phenyllactate PLA$", "Phenyllactate"), "^Dopamine 3-O-sulfate$",
        "Dopamine 3-O-sulfate")
}
assign_metabolic_theme <- function(name) {
    z <- stringr::str_to_lower(as.character(name))
    dplyr::case_when(stringr::str_detect(z, paste0("phenyl|phenol|hydroquinone|tyros|tyramine|dopamine|",
        "norepinephrine|epinephrine|homovanil|hydroxycinnam|", "hydroxybenz|hydroxyphenyl|methoxytyros|l-dopa|ldopa|",
        "homogentis|indoxyl|indole|norsalsolinol|normetanephrine")) ~ "Aromatic amino-acid metabolism",
        stringr::str_detect(z, "histid|histamin|urocan|imidazole|carnosine|formiminoglut") ~ "Histidine and imidazole metabolism",
        stringr::str_detect(z, paste0("isoleuc|leucine|leucyl|valine|ketoisovaler|", "methylcrotonoyl|methyl-2-oxovaler|isobutyryl|",
            "ketoleuc|acetylisoleuc")) ~ "Branched-chain amino-acid metabolism", stringr::str_detect(z,
            "lysine|pipecol|aminoadip|oxoadip|hydroxyadip|homocitrull") ~ "Lysine and aminoadipate metabolism",
        stringr::str_detect(z, paste0("threonine|acetyl-threonine|acetylthreonine|", "acetylserine|methylthio|methion")) ~
            "Other amino-acid metabolism", stringr::str_detect(z, paste0("cholate|cholic|deoxychol|ursodeoxychol|chenodeoxychol|",
            "taurocholen|glycoursodeoxychol|bile|xol24")) ~ "Bile-acid metabolism", stringr::str_detect(z,
            "androster|epiandro|androstan|androstenediol|steroid hormone") ~ "Steroid hormone conjugates",
        stringr::str_detect(z, paste0("ascorbic|niacin|nicotin|nicotinamide|methylnicotinamide|", "pyridone.*carboxamide|pyridox|thiamine|vitamin|",
            "ubiquinone|qh2|ergothioneine|superoxide")) ~ "Vitamins, cofactors and redox", stringr::str_detect(z,
            paste0("^sulfate$|^sulfite$|cysteic|cysteine-s-sulfate|", "chondroitin|heparan|proteoglycan")) ~
            "Sulfur and glycosaminoglycan metabolism", stringr::str_detect(z, "citrate|aconitate|acetoacetate") ~
            "Central carbon and ketone metabolism", stringr::str_detect(z, "choline|glycerophosphoinositol|gpe|phosphatid|plasmalogen") ~
            "Phospholipid and choline metabolism", stringr::str_detect(z, "cytosine|acetylglucosamine|acetylgalactosamine|amino sugar") ~
            "Nucleotide and amino-sugar metabolism", TRUE ~ "Other or unclassified")
}
is_medication <- function(x) {
    stringr::str_detect(stringr::str_to_lower(as.character(x)), "omeprazole|pantoprazole|dexlansoprazole|olmesartan")
}
is_unidentified <- function(x) {
    stringr::str_detect(trimws(as.character(x)), "^[Xx]-[0-9]+$")
}
is_uninformative_wbm_name <- function(x) {
    stringr::str_detect(stringr::str_to_lower(as.character(x)), "^xol[0-9]")
}
representative_name_priority <- function(name, measured_canonical_names) {
    canon <- canonical_name(name)
    direct_match <- canon %in% measured_canonical_names
    simple_name_bonus <- dplyr::case_when(stringr::str_detect(stringr::str_to_lower(name), "phenyllactate|dopamine 3-o-sulfate|phenylacetylglutamine|l-phenylalanine") ~
        6, stringr::str_detect(stringr::str_to_lower(name), "l-isoleucine|l-leucine|l-valine|l-histidine|histamine|l-lysine|sulfate|niacinamide|thiamine") ~
        4, stringr::str_detect(name, "Proteoglycan|^xol") ~ -5, TRUE ~ 0)
    10 * as.numeric(direct_match) + simple_name_bonus - nchar(as.character(name))/1000
}
build_ordered_labels <- function(df) {
    if (nrow(df) == 0L) {
        return(df)
    }
    dplyr::mutate(dplyr::arrange(df, factor(.data$theme, levels = theme_order), .data$effect, .data$FDR,
        .data$display_label), theme = factor(.data$theme, levels = theme_order), plot_label = make.unique(as.character(.data$display_label)),
        plot_label = factor(.data$plot_label, levels = rev(unique(.data$plot_label))))
}
select_top_per_theme <- function(df, max_per_theme) {
    if (nrow(df) == 0L) {
        return(df)
    }
    dplyr::ungroup(dplyr::filter(dplyr::mutate(dplyr::arrange(dplyr::group_by(df, .data$theme), dplyr::desc(.data$force_plot),
        .data$FDR, dplyr::desc(abs(.data$effect)), .by_group = TRUE), theme_rank = dplyr::row_number()),
        .data$force_plot | .data$theme_rank <= max_per_theme))
}
if (!project_file_exists(wbm_file)) {
    stop("WBM input file not found: ", wbm_file)
}
if (!project_file_exists(metabolomics_file)) {
    stop("Metabolomics input file not found: ", metabolomics_file)
}
wbm_raw <- readr::read_csv(project_input(wbm_file), show_col_types = FALSE)
met_raw <- readr::read_csv(project_input(metabolomics_file), show_col_types = FALSE)
assert_columns(wbm_raw, c("reaction", "met_abbr", "metabolite_name", "same_direction", "beta_kiel", "beta_metacardis",
    "beta_meta", "lower_meta", "upper_meta", "Q_pvalue", "I2", "FDR_meta"), "WBM results")
assert_columns(met_raw, c("metabolite", "hc3_coef", "hc3_ci_low", "hc3_ci_high", "hc3_fdr"), "Measured metabolomics results")
measured_significant_names <- dplyr::pull(dplyr::filter(met_raw, is.finite(.data$hc3_fdr), .data$hc3_fdr <
    fdr_cutoff), .data$metabolite)
measured_canonical_names <- unique(canonical_name(measured_significant_names))
wbm_annotated <- dplyr::mutate(wbm_raw, same_direction = as_logical_safe(.data$same_direction), source = "Predicted blood flux",
    original_name = as.character(.data$metabolite_name), display_label = clean_display_name(.data$metabolite_name),
    effect = as.numeric(.data$beta_meta), CI_low = as.numeric(.data$lower_meta), CI_high = as.numeric(.data$upper_meta),
    FDR = as.numeric(.data$FDR_meta), theme = assign_metabolic_theme(.data$metabolite_name), significant = is.finite(.data$FDR) &
        .data$FDR < fdr_cutoff, robust_wbm = (.data$significant & .data$same_direction & is.finite(.data$Q_pvalue) &
        .data$Q_pvalue >= wbm_min_q_pvalue & is.finite(.data$I2) & .data$I2 <= wbm_max_i2), medication = is_medication(.data$metabolite_name),
    unidentified = is_unidentified(.data$metabolite_name), cryptic_wbm_name = is_uninformative_wbm_name(.data$metabolite_name),
    auto_include_main = (.data$robust_wbm & .data$theme != "Other or unclassified" & abs(.data$effect) >=
        minimum_abs_effect_for_main_wbm & !.data$medication & !.data$unidentified & !.data$cryptic_wbm_name),
    force_plot = canonical_name(.data$metabolite_name) %in% canonical_name(c("Phenyllactate", "Dopamine 3-O-sulfate")),
    cross_platform_exact_name = canonical_name(.data$metabolite_name) %in% measured_canonical_names)
met_annotated <- dplyr::mutate(met_raw, source = "Measured serum abundance", original_name = as.character(.data$metabolite),
    display_label = clean_display_name(.data$metabolite), effect = as.numeric(.data$hc3_coef), CI_low = as.numeric(.data$hc3_ci_low),
    CI_high = as.numeric(.data$hc3_ci_high), FDR = as.numeric(.data$hc3_fdr), theme = assign_metabolic_theme(.data$metabolite),
    significant = is.finite(.data$FDR) & .data$FDR < fdr_cutoff, medication = is_medication(.data$metabolite),
    unidentified = is_unidentified(.data$metabolite), auto_include_main = (.data$significant & .data$theme !=
        "Other or unclassified" & (!exclude_medications_from_main | !.data$medication) & (!exclude_unidentified_from_main |
        !.data$unidentified)), force_plot = canonical_name(.data$metabolite) %in% canonical_name(c("phenyllactate.(PLA)",
        "dopamine.3-O-sulfate", "4-hydroxyphenylacetate", "4-hydroxyphenylacetylglutamine", "ursodeoxycholate",
        "cholate")), cross_platform_exact_name = canonical_name(.data$metabolite) %in% canonical_name(wbm_annotated$metabolite_name[wbm_annotated$robust_wbm]))
curation_template <- dplyr::arrange(dplyr::distinct(dplyr::bind_rows(dplyr::transmute(dplyr::filter(wbm_annotated,
    .data$robust_wbm), source, original_name, auto_theme = .data$theme, theme = .data$theme, auto_display_label = .data$display_label,
    display_label = .data$display_label, auto_include_main, include_main = .data$auto_include_main, force_plot,
    notes = ""), dplyr::transmute(dplyr::filter(met_annotated, .data$significant), source, original_name,
    auto_theme = .data$theme, theme = .data$theme, auto_display_label = .data$display_label, display_label = .data$display_label,
    auto_include_main, include_main = .data$auto_include_main, force_plot, notes = "")), .data$source,
    .data$original_name, .keep_all = TRUE), .data$source, .data$theme, .data$original_name)
if (!project_file_exists(manual_curation_file)) {
    readr::write_csv(curation_template, project_output(manual_curation_file), na = "")
    message("Created optional manual curation template: ", manual_curation_file, ". The automatic analysis will continue. Edit and rerun to override themes, labels, or inclusion.")
}
manual_curation <- readr::read_csv(project_input(manual_curation_file), show_col_types = FALSE)
assert_columns(manual_curation, c("source", "original_name", "theme", "display_label", "include_main",
    "force_plot"), "Manual curation file")
manual_curation <- dplyr::transmute(manual_curation, source = as.character(.data$source), original_name = as.character(.data$original_name),
    manual_theme = as.character(.data$theme), manual_display_label = as.character(.data$display_label),
    manual_include_main = as_logical_safe(.data$include_main), manual_force_plot = as_logical_safe(.data$force_plot))
apply_manual_curation <- function(x) {
    dplyr::select(dplyr::mutate(dplyr::left_join(x, manual_curation, by = c("source", "original_name")),
        theme = dplyr::if_else(!is.na(.data$manual_theme) & nzchar(.data$manual_theme), .data$manual_theme,
            .data$theme), display_label = dplyr::if_else(!is.na(.data$manual_display_label) & nzchar(.data$manual_display_label),
            .data$manual_display_label, .data$display_label), include_main = dplyr::coalesce(.data$manual_include_main,
            .data$auto_include_main), force_plot = dplyr::coalesce(.data$manual_force_plot, .data$force_plot)),
        -dplyr::starts_with("manual_"))
}
wbm_annotated <- apply_manual_curation(wbm_annotated)
met_annotated <- apply_manual_curation(met_annotated)
theme_order <- c("Aromatic amino-acid metabolism", "Histidine and imidazole metabolism", "Branched-chain amino-acid metabolism",
    "Lysine and aminoadipate metabolism", "Other amino-acid metabolism", "Bile-acid metabolism", "Steroid hormone conjugates",
    "Sulfur and glycosaminoglycan metabolism", "Vitamins, cofactors and redox", "Central carbon and ketone metabolism",
    "Phospholipid and choline metabolism", "Nucleotide and amino-sugar metabolism", "Other or unclassified")
extra_themes <- setdiff(unique(c(wbm_annotated$theme, met_annotated$theme)), theme_order)
theme_order <- c(theme_order, sort(extra_themes))
wbm_robust <- dplyr::filter(wbm_annotated, .data$robust_wbm)
if (collapse_wbm_near_identical_estimates) {
    wbm_robust <- dplyr::mutate(wbm_robust, estimate_cluster_key = paste(.data$theme, round(.data$effect,
        wbm_cluster_digits), round(.data$CI_low, wbm_cluster_digits), round(.data$CI_high, wbm_cluster_digits),
        sep = "|"))
} else {
    wbm_robust <- dplyr::mutate(wbm_robust, estimate_cluster_key = .data$reaction)
}
wbm_cluster_members <- dplyr::ungroup(dplyr::mutate(dplyr::group_by(wbm_robust, .data$estimate_cluster_key),
    cluster_size = dplyr::n(), representative_priority = representative_name_priority(.data$original_name,
        measured_canonical_names), representative_rank = rank(-(.data$representative_priority - log10(pmax(.data$FDR,
        1e-300))/1000), ties.method = "first")))
wbm_representatives_all <- dplyr::mutate(dplyr::ungroup(dplyr::slice(dplyr::arrange(dplyr::group_by(wbm_cluster_members,
    .data$estimate_cluster_key), dplyr::desc(.data$force_plot), dplyr::desc(.data$cross_platform_exact_name),
    dplyr::desc(.data$representative_priority), .data$FDR, .by_group = TRUE), 1L)), cluster_members = vapply(.data$estimate_cluster_key,
    function(k) {
        paste(unique(wbm_cluster_members$display_label[wbm_cluster_members$estimate_cluster_key == k]),
            collapse = "; ")
    }, character(1)))
wbm_plot_data <- build_ordered_labels(select_top_per_theme(dplyr::filter(wbm_representatives_all, .data$include_main),
    max_wbm_representatives_per_theme))
met_plot_data <- build_ordered_labels(select_top_per_theme(dplyr::filter(met_annotated, .data$significant,
    .data$include_main), max_measured_metabolites_per_theme))
if (nrow(wbm_plot_data) == 0L) {
    stop("No WBM results remained for the main figure after filtering/curation.")
}
if (nrow(met_plot_data) == 0L) {
    stop("No measured metabolomics results remained for the main figure after filtering/curation.")
}
technical_plot_data <- dplyr::mutate(dplyr::filter(wbm_annotated, is.finite(.data$beta_kiel), is.finite(.data$beta_metacardis)),
    technical_status = dplyr::if_else(.data$robust_wbm, "Robust meta-analysis signal", "Other tested reaction"),
    same_cohort_direction = dplyr::case_when(.data$beta_kiel == 0 | .data$beta_metacardis == 0 ~ NA,
        TRUE ~ sign(.data$beta_kiel) == sign(.data$beta_metacardis)))
technical_complete <- dplyr::filter(technical_plot_data, is.finite(.data$beta_kiel), is.finite(.data$beta_metacardis))
technical_concordance <- tibble::tibble(metric = c("Number of reactions", "Pearson correlation", "Spearman correlation",
    "Same-direction fraction", "Robust reactions", "Robust same-direction fraction"), value = c(nrow(technical_complete),
    stats::cor(technical_complete$beta_kiel, technical_complete$beta_metacardis, method = "pearson"),
    stats::cor(technical_complete$beta_kiel, technical_complete$beta_metacardis, method = "spearman"),
    mean(technical_complete$same_cohort_direction, na.rm = TRUE), sum(technical_complete$robust_wbm,
        na.rm = TRUE), mean(technical_complete$same_cohort_direction[technical_complete$robust_wbm],
        na.rm = TRUE)))
technical_pearson <- technical_concordance$value[technical_concordance$metric == "Pearson correlation"]
technical_spearman <- technical_concordance$value[technical_concordance$metric == "Spearman correlation"]
technical_same_direction <- technical_concordance$value[technical_concordance$metric == "Same-direction fraction"]
source_colors <- c(`Higher in PPI` = "#5A4E6D", `Lower in PPI` = "#4C7A9F")
pal_sig <- c(`Robust meta-analysis signal` = "#5A4E6D", `Other tested reaction` = "#D3D3D3")
main_theme <- ggplot2::theme_bw(base_size = 14, base_family = "Arial") + ggplot2::theme(panel.grid = ggplot2::element_blank(),
    plot.title = ggplot2::element_blank(), plot.subtitle = ggplot2::element_blank(), axis.text = ggplot2::element_text(colour = "#2F2F2F"),
    axis.title = ggplot2::element_text(colour = "#2F2F2F"), legend.position = "bottom", legend.title = ggplot2::element_blank())
forest_theme <- main_theme + ggplot2::theme(panel.grid = ggplot2::element_blank(), strip.background = ggplot2::element_rect(fill = "#D9D9D9",
    colour = "#3A3A3A", linewidth = 0.4), strip.text.y.left = ggplot2::element_text(angle = 0, face = "bold",
    size = 7.5, lineheight = 0.9, colour = "#2F2F2F"), axis.text.x = ggplot2::element_text(size = 8),
    axis.text.y = ggplot2::element_text(size = 7.2), axis.title.x = ggplot2::element_text(size = 9),
    axis.title.y = ggplot2::element_blank(), legend.text = ggplot2::element_text(size = 8), panel.spacing.y = grid::unit(1.5,
        "mm"), plot.margin = ggplot2::margin(5.5, 5.5, 5.5, 5.5))
theme_strip_labeller <- ggplot2::as_labeller(c(`Aromatic amino-acid metabolism` = "Aromatic\namino acids",
    `Histidine and imidazole metabolism` = "Histidine /\nimidazole", `Branched-chain amino-acid metabolism` = "BCAA",
    `Lysine and aminoadipate metabolism` = "Lysine /\naminoadipate", `Other amino-acid metabolism` = "Other AA",
    `Bile-acid metabolism` = "Bile acids", `Steroid hormone conjugates` = "Steroid hormone\nconjugates",
    `Sulfur and glycosaminoglycan metabolism` = "Sulfur /\nglycosaminoglycan", `Vitamins, cofactors and redox` = "Vitamins, cofactors\nand redox",
    `Central carbon and ketone metabolism` = "Central carbon /\nketones", `Phospholipid and choline metabolism` = "Phospholipid/choline",
    `Nucleotide and amino-sugar metabolism` = "Nucleotide /\namino sugar", `Other or unclassified` = "Other"),
    default = ggplot2::label_value)
p_wbm <- ggplot2::ggplot(wbm_plot_data, ggplot2::aes(x = .data$effect, y = .data$plot_label, colour = dplyr::if_else(.data$effect >
    0, "Higher in PPI", "Lower in PPI"))) + ggplot2::geom_vline(xintercept = 0, linetype = "dashed",
    linewidth = 0.45, colour = "#8C8C8C") + ggplot2::geom_errorbar(ggplot2::aes(xmin = .data$CI_low,
    xmax = .data$CI_high), orientation = "y", width = 0, linewidth = 0.45) + ggplot2::geom_point(ggplot2::aes(shape = .data$cross_platform_exact_name,
    size = pmin(-log10(pmax(.data$FDR, 1e-300)), 8)), alpha = 0.9, stroke = 0.7) + ggplot2::facet_grid(rows = ggplot2::vars(theme),
    scales = "free_y", space = "free_y", switch = "y", labeller = ggplot2::labeller(theme = theme_strip_labeller)) +
    ggplot2::scale_colour_manual(values = source_colors, breaks = c("Higher in PPI", "Lower in PPI"),
        labels = c("Higher in PPI users", "Lower in PPI users"), drop = FALSE) + ggplot2::scale_shape_manual(values = c(`FALSE` = 16,
    `TRUE` = 18), labels = c(`FALSE` = "Layer-specific", `TRUE` = "Exact metabolite-name match")) + ggplot2::scale_size_continuous(range = c(1.8,
    3.2), guide = "none") + ggplot2::labs(x = "Blood WBM meta-analysis effect size (95% CI)", y = NULL,
    colour = NULL, shape = NULL) + forest_theme + ggplot2::theme(plot.margin = ggplot2::margin(t = 5.5,
    r = 3, b = 5.5, l = 0), axis.text.y = ggplot2::element_text(size = 7.2, margin = ggplot2::margin(r = 1)))
p_met <- ggplot2::ggplot(met_plot_data, ggplot2::aes(x = .data$effect, y = .data$plot_label, colour = dplyr::if_else(.data$effect >
    0, "Higher in PPI", "Lower in PPI"))) + ggplot2::geom_vline(xintercept = 0, linetype = "dashed",
    linewidth = 0.45, colour = "#8C8C8C") + ggplot2::geom_errorbar(ggplot2::aes(xmin = .data$CI_low,
    xmax = .data$CI_high), orientation = "y", width = 0, linewidth = 0.45) + ggplot2::geom_point(ggplot2::aes(shape = .data$cross_platform_exact_name,
    size = pmin(-log10(pmax(.data$FDR, 1e-300)), 8)), alpha = 0.9, stroke = 0.7) + ggplot2::facet_grid(rows = ggplot2::vars(theme),
    scales = "free_y", space = "free_y", switch = "y", labeller = ggplot2::labeller(theme = theme_strip_labeller)) +
    ggplot2::scale_colour_manual(values = source_colors, breaks = c("Higher in PPI", "Lower in PPI"),
        labels = c("Higher in PPI users", "Lower in PPI users"), drop = FALSE) + ggplot2::scale_shape_manual(values = c(`FALSE` = 16,
    `TRUE` = 18), labels = c(`FALSE` = "Layer-specific", `TRUE` = "Exact metabolite-name match")) + ggplot2::scale_size_continuous(range = c(1.8,
    3.2), guide = "none") + ggplot2::labs(x = "Adjusted serum-metabolite effect size (95% CI)", y = NULL,
    colour = NULL, shape = NULL) + forest_theme
technical_plot_data <- dplyr::mutate(technical_plot_data, technical_status = factor(.data$technical_status,
    levels = c("Other tested reaction", "Robust meta-analysis signal")))
technical_labels <- dplyr::slice_head(dplyr::arrange(dplyr::filter(technical_plot_data, .data$robust_wbm),
    .data$FDR, dplyr::desc(abs(.data$beta_meta))), n = 3L)
p_technical <- ggplot2::ggplot(technical_plot_data, ggplot2::aes(x = .data$beta_kiel, y = .data$beta_metacardis)) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.45, colour = "#9E9E9E") +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.45, colour = "#9E9E9E") +
    ggplot2::geom_abline(slope = 1, intercept = 0, linewidth = 0.5, colour = "#BDBDBD") + ggplot2::geom_point(ggplot2::aes(colour = .data$technical_status),
    size = 2.2, alpha = 0.9) + ggplot2::geom_point(data = dplyr::filter(technical_plot_data, .data$robust_wbm,
    .data$cross_platform_exact_name), shape = 23, fill = "#5A4E6D", colour = "#2F2F2F", stroke = 0.55,
    size = 3, show.legend = FALSE) + ggrepel::geom_text_repel(data = technical_labels, ggplot2::aes(label = stringr::str_trunc(.data$display_label,
    35)), size = 3, nudge_y = 0.012, force = 8, box.padding = grid::unit(1.5, "mm"), point.padding = grid::unit(2.5,
    "mm"), min.segment.length = 0, max.overlaps = Inf, max.time = 2, seed = 123, segment.colour = "#B0B0B0",
    segment.size = 0.35, colour = "#2F2F2F", family = "Arial", show.legend = FALSE) + ggplot2::annotate("text",
    x = min(technical_plot_data$beta_kiel, na.rm = TRUE), y = max(technical_plot_data$beta_metacardis,
        na.rm = TRUE), hjust = 0, vjust = 1, size = 4, family = "Arial", colour = "#2F2F2F", label = paste0("Spearman ρ = ",
        formatC(technical_spearman, digits = 2, format = "f"), "\nConcordance = ", formatC(100 * technical_same_direction,
            digits = 1, format = "f"), "%")) + ggplot2::scale_colour_manual(values = pal_sig, breaks = c("Robust meta-analysis signal",
    "Other tested reaction"), labels = c("Meta-analysis significant", "Not significant")) + ggplot2::labs(x = "Kiel blood WBM effect size",
    y = "MetaCardis WBM effect size", colour = NULL) + main_theme + ggplot2::theme(legend.position = "bottom",
    legend.text = ggplot2::element_text(size = 9), axis.title.x = ggplot2::element_text(size = 10), axis.title.y = ggplot2::element_text(size = 10),
    plot.margin = ggplot2::margin(5.5, 5.5, 5.5, 5.5))
forest_block <- ((p_wbm | p_met)/patchwork::guide_area()) + patchwork::plot_layout(heights = c(1, 0.12),
    guides = "collect") & ggplot2::theme(legend.position = "bottom", legend.justification = "center",
    legend.box.just = "center", legend.box = "horizontal", legend.direction = "horizontal", legend.text = ggplot2::element_text(size = 8),
    legend.spacing.x = grid::unit(5, "mm"), legend.margin = ggplot2::margin(1, 0, 1, 0))
technical_row <- patchwork::free(p_technical)
combined_figure <- (technical_row/forest_block) + patchwork::plot_layout(heights = c(2,
    4.2), guides = "keep") + patchwork::plot_annotation(tag_levels = "A", theme = ggplot2::theme(plot.tag = ggplot2::element_text(family = "Arial",
    face = "bold", size = 14, colour = "#2F2F2F")))
write_csv_safe(dplyr::arrange(dplyr::filter(wbm_annotated, .data$robust_wbm), .data$FDR), "01_WBM_robust_annotated.csv")
write_csv_safe(dplyr::arrange(dplyr::filter(met_annotated, .data$significant), .data$FDR), "02_metabolomics_significant_annotated.csv")
write_csv_safe(dplyr::arrange(wbm_cluster_members, .data$theme, .data$estimate_cluster_key, .data$FDR),
    "03_WBM_collapsed_cluster_members.csv")
write_csv_safe(dplyr::select(wbm_plot_data, theme, display_label, reaction, met_abbr, effect, CI_low,
    CI_high, FDR, cluster_size, cluster_members, cross_platform_exact_name, force_plot), "04_WBM_plot_representatives.csv")
write_csv_safe(dplyr::select(met_plot_data, theme, display_label, effect, CI_low, CI_high, FDR, cross_platform_exact_name,
    force_plot), "05_metabolomics_plot_representatives.csv")
write_csv_safe(dplyr::select(technical_plot_data, reaction, met_abbr, metabolite_name, beta_kiel, beta_metacardis,
    beta_meta, lower_meta, upper_meta, FDR_kiel, FDR_metacardis, FDR_meta, Q_pvalue, I2, same_direction,
    robust_wbm, theme, cross_platform_exact_name, technical_status), "06_WBM_cohort_coefficient_plot_data.csv")
write_csv_safe(technical_concordance, "06_WBM_cohort_coefficient_concordance.csv")
excluded_audit <- dplyr::bind_rows(dplyr::transmute(dplyr::filter(wbm_annotated, .data$robust_wbm, !.data$include_main),
    source, original_name, display_label, theme, FDR, effect, medication, unidentified, reason = dplyr::case_when(.data$medication ~
        "Medication", .data$unidentified ~ "Unidentified feature", .data$cryptic_wbm_name ~ "Cryptic WBM metabolite identifier",
        .data$theme == "Other or unclassified" ~ "Unclassified theme", abs(.data$effect) < minimum_abs_effect_for_main_wbm ~
            "Numerically negligible WBM effect", TRUE ~ "Excluded by manual curation")), dplyr::transmute(dplyr::filter(met_annotated,
    .data$significant, !.data$include_main), source, original_name, display_label, theme, FDR, effect,
    medication, unidentified, reason = dplyr::case_when(.data$medication ~ "Medication", .data$unidentified ~
        "Unidentified feature", .data$theme == "Other or unclassified" ~ "Unclassified theme", TRUE ~
        "Excluded by manual curation")))
write_csv_safe(dplyr::arrange(excluded_audit, .data$source, .data$reason, .data$FDR), "08_unclassified_or_excluded_results.csv")
ggplot2::ggsave(filename = project_output(file.path(output_dir, "MetaCardis_WBM_metabolomics_integrated_theme_figure.pdf"),
    name = "meta_cardis_wbm_metabolomics_integrated_theme_figure.pdf"), plot = combined_figure, width = 350,
    height = 350, units = "mm", device = grDevices::cairo_pdf, limitsize = FALSE)
ggplot2::ggsave(filename = project_output(file.path(output_dir, "MetaCardis_WBM_metabolomics_integrated_theme_figure.png"),
    name = "meta_cardis_wbm_metabolomics_integrated_theme_figure.png"), plot = combined_figure, width = 350,
    height = 350, units = "mm", dpi = figure_dpi, bg = "white", limitsize = FALSE)
analysis_summary <- tibble::tibble(item = c("WBM input rows", "Robust WBM rows before collapsing", "WBM coupled estimate clusters",
    "WBM representatives shown", "Measured metabolomics input rows", "Measured metabolites at FDR cutoff",
    "Measured metabolites shown", "Exact cross-platform names among shown WBM representatives", "Exact cross-platform names among shown measured metabolites",
    "Technical coefficient pairs", "Technical Pearson correlation", "Technical Spearman correlation",
    "Technical same-direction fraction"), value = c(nrow(wbm_raw), nrow(wbm_robust), dplyr::n_distinct(wbm_cluster_members$estimate_cluster_key),
    nrow(wbm_plot_data), nrow(met_raw), sum(met_annotated$significant, na.rm = TRUE), nrow(met_plot_data),
    sum(wbm_plot_data$cross_platform_exact_name, na.rm = TRUE), sum(met_plot_data$cross_platform_exact_name,
        na.rm = TRUE), nrow(technical_complete), technical_pearson, technical_spearman, technical_same_direction))
write_csv_safe(analysis_summary, "09_analysis_summary.csv")
message("Analysis complete. Outputs written to: ", project_normalize_path(output_dir))
print(analysis_summary)
}))
