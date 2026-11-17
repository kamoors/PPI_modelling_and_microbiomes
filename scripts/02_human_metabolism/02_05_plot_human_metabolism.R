# SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0
# Figure 2: blood-flux and serum-metabolite associations.

source("scripts/helpers/setup.R")
invisible(run_analysis(
  section = "02_human_metabolism",
  directory = "figures/scripts",
  script = "02_05_plot_human_metabolism.R",
  restricted = FALSE,
  code = {
required_packages <- c("dplyr", "readr", "stringr", "forcats", "ggplot2", "ggrepel", "cowplot", "grid")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) {
    stop("Install required packages first: ", paste(missing_packages, collapse = ", "))
}
project_check_plot_versions(required_packages)
FIG_ROOT <- project_path("figures")
DATA_DIR <- file.path(FIG_ROOT, "source_data", "WBM_Serum_Metabolomics")
OUT_DIR <- file.path(FIG_ROOT, "final")
dir.create(project_output(DATA_DIR, directory = TRUE), recursive = TRUE, showWarnings = FALSE)
dir.create(project_output(OUT_DIR, directory = TRUE), recursive = TRUE, showWarnings = FALSE)
REFRESH_SOURCE_DATA <- FALSE
PPI_ROOT <- project_path("")
LEGACY_DIR <- file.path(PPI_ROOT, "Results", "Kiel_Metacardis", "MetaCardis_WBM_metabolomics_integrated_figure")
stage_input <- function(filename, candidates) {
    target <- file.path(DATA_DIR, filename)
    if (REFRESH_SOURCE_DATA || !project_file_exists(target)) {
        hits <- candidates[project_file_exists(candidates)]
        if (!length(hits)) {
            stop("Could not stage ", filename, ".\nChecked:\n", paste(candidates, collapse = "\n"))
        }
        if (!file.copy(project_input(hits[[1]]), project_output(target), overwrite = TRUE)) {
            stop("Failed to copy ", hits[[1]])
        }
        message("Staged: ", filename)
    }
    target
}
wbm_plot_file <- stage_input("04_WBM_plot_representatives.csv", c(file.path(LEGACY_DIR, "04_WBM_plot_representatives.csv")))
met_plot_file <- stage_input("05_metabolomics_plot_representatives.csv", c(file.path(LEGACY_DIR, "05_metabolomics_plot_representatives.csv")))
technical_file <- stage_input("06_WBM_cohort_coefficient_plot_data.csv", c(file.path(LEGACY_DIR, "06_WBM_cohort_coefficient_plot_data.csv")))
concordance_file <- stage_input("06_WBM_cohort_coefficient_concordance.csv", c(file.path(LEGACY_DIR,
    "06_WBM_cohort_coefficient_concordance.csv")))
readr::write_csv(data.frame(file = basename(c(wbm_plot_file, met_plot_file, technical_file, concordance_file)),
    centralized_path = c(wbm_plot_file, met_plot_file, technical_file, concordance_file)), project_output(file.path(DATA_DIR,
    "WBM_Serum_Metabolomics_source_manifest.csv"), name = "wbm_serum_metabolomics_source_manifest.csv"))
BASE_SIZE <- 18
FIG_W_IN <- 18
FIG_H_IN <- 18
DPI <- 600
source_colors <- c(`Higher in PPI` = "#5A4E6D", `Lower in PPI` = "#4C7A9F")
pal_sig <- c(`Robust meta-analysis signal` = "#5A4E6D", `Other tested reaction` = "#D3D3D3")
main_theme <- ggplot2::theme_bw(base_size = BASE_SIZE, base_family = "Arial") + ggplot2::theme(panel.grid = ggplot2::element_blank(),
    plot.title = ggplot2::element_blank(), axis.text = ggplot2::element_text(size = 17, colour = "#2F2F2F"),
    axis.title = ggplot2::element_text(size = 18, colour = "#2F2F2F"), legend.position = "bottom", legend.title = ggplot2::element_blank(),
    legend.text = ggplot2::element_text(size = 15), legend.key = ggplot2::element_blank(), legend.background = ggplot2::element_blank(),
    plot.margin = ggplot2::margin(7, 7, 7, 7))
forest_theme <- main_theme + ggplot2::theme(strip.background = ggplot2::element_rect(fill = "#D9D9D9",
    colour = "#3A3A3A", linewidth = 0.35), strip.text.y.left = ggplot2::element_text(angle = 0, face = "bold",
    size = 13, lineheight = 0.85), axis.text.x = ggplot2::element_text(size = 10), axis.text.y = ggplot2::element_text(size = 15),
    axis.title.x = ggplot2::element_text(size = 17), axis.title.y = ggplot2::element_blank(), panel.spacing.y = grid::unit(1.2,
        "mm"))
theme_order <- c("Aromatic amino-acid metabolism", "Histidine and imidazole metabolism", "Branched-chain amino-acid metabolism",
    "Lysine and aminoadipate metabolism", "Other amino-acid metabolism", "Bile-acid metabolism", "Steroid hormone conjugates",
    "Sulfur and glycosaminoglycan metabolism", "Vitamins, cofactors and redox", "Central carbon and ketone metabolism",
    "Phospholipid and choline metabolism", "Nucleotide and amino-sugar metabolism", "Other or unclassified")
theme_strip_labeller <- ggplot2::as_labeller(c(`Aromatic amino-acid metabolism` = "Aromatic\namino acids",
    `Histidine and imidazole metabolism` = "Histidine /\nimidazole", `Branched-chain amino-acid metabolism` = "BCAA",
    `Lysine and aminoadipate metabolism` = "Lysine /\naminoadipate", `Other amino-acid metabolism` = "Other AA",
    `Bile-acid metabolism` = "Bile acids", `Steroid hormone conjugates` = "Steroid hormone\nconjugates",
    `Sulfur and glycosaminoglycan metabolism` = "Sulfur /\nglycosaminoglycan", `Vitamins, cofactors and redox` = "Vitamins, cofactors\nand redox",
    `Central carbon and ketone metabolism` = "Central carbon /\nketones", `Phospholipid and choline metabolism` = "Phospholipid/choline",
    `Nucleotide and amino-sugar metabolism` = "Nucleotide /\namino sugar", `Other or unclassified` = "Other"),
    default = ggplot2::label_value)
as_logical_safe <- function(x) {
    if (is.logical(x))
        return(x)
    tolower(as.character(x)) %in% c("true", "t", "1", "yes")
}
capitalize_first_letter <- function(x) {
    sub("([A-Za-z])", "\\U\\1", x, perl = TRUE)
}
prepare_forest <- function(df) {
    dplyr::mutate(dplyr::arrange(dplyr::mutate(df, theme = factor(theme, levels = theme_order), cross_platform_exact_name = as_logical_safe(cross_platform_exact_name),
        direction = ifelse(effect > 0, "Higher in PPI", "Lower in PPI"), plot_label_chr = as.character(display_label)),
        theme, effect, FDR, plot_label_chr), plot_label_chr = make.unique(plot_label_chr), plot_label = factor(plot_label_chr,
        levels = rev(unique(plot_label_chr))))
}
wbm_plot_data <- prepare_forest(dplyr::mutate(readr::read_csv(project_input(wbm_plot_file), show_col_types = FALSE),
    display_label = dplyr::case_when(stringr::str_detect(display_label, "^Chondroitin Sulfate E") ~ "Chondroitin sulfate E proteoglycan",
        stringr::str_detect(display_label, "^Chondroitin Sulfate D") ~ "Chondroitin sulfate D proteoglycan",
        stringr::str_detect(display_label, "^Chondroitin Sulfate C") ~ "Chondroitin sulfate C proteoglycan",
        TRUE ~ display_label)))
met_plot_data <- prepare_forest(dplyr::mutate(readr::read_csv(project_input(met_plot_file), show_col_types = FALSE),
    display_label = capitalize_first_letter(display_label), theme = dplyr::if_else(theme == "Lysine and aminoadipate metabolism",
        "Other amino-acid metabolism", theme)))
make_forest <- function(df, xlab) {
    ggplot2::ggplot(df, ggplot2::aes(x = effect, y = plot_label, colour = direction)) + ggplot2::geom_vline(xintercept = 0,
        linetype = "dashed", linewidth = 0.4, colour = "#8C8C8C") + ggplot2::geom_errorbar(ggplot2::aes(xmin = CI_low,
        xmax = CI_high), orientation = "y", width = 0, linewidth = 0.42) + ggplot2::geom_point(ggplot2::aes(shape = cross_platform_exact_name,
        size = pmin(-log10(pmax(FDR, 1e-300)), 8)), alpha = 0.9, stroke = 0.6) + ggplot2::facet_grid(rows = ggplot2::vars(theme),
        scales = "free_y", space = "free_y", switch = "y", labeller = ggplot2::labeller(theme = theme_strip_labeller)) +
        ggplot2::scale_colour_manual(values = source_colors, breaks = c("Higher in PPI", "Lower in PPI"),
            labels = c("Higher in PPI users", "Lower in PPI users")) + ggplot2::scale_shape_manual(values = c(`FALSE` = 16,
        `TRUE` = 18), breaks = c(FALSE, TRUE), labels = c(`FALSE` = "Layer-specific", `TRUE` = "Exact metabolite-name match")) +
        ggplot2::scale_size_continuous(range = c(1.7, 3), guide = "none") + ggplot2::guides(shape = ggplot2::guide_legend(order = 1,
        nrow = 1), colour = ggplot2::guide_legend(order = 2, nrow = 1, override.aes = list(size = 2.4,
        alpha = 1))) + ggplot2::labs(x = xlab, y = NULL, colour = NULL, shape = NULL) + forest_theme
}
p_wbm <- make_forest(wbm_plot_data, "Meta-analysis effect size (95% CI)") + ggplot2::theme(legend.position = "none")
p_met <- make_forest(met_plot_data, "Effect size (95% CI)") + ggplot2::theme(legend.position = "none")
forest_legend <- ggplot2::ggplot() + ggplot2::coord_cartesian(xlim = c(0, 12.5), ylim = c(0, 1), clip = "off") +
    ggplot2::geom_point(ggplot2::aes(x = 0.35, y = 0.5), shape = 16, size = 3.2, colour = "#2F2F2F") +
    ggplot2::annotate("text", x = 0.65, y = 0.5, label = "Layer-specific", hjust = 0, vjust = 0.5, size = 5,
        family = "Arial", colour = "#2F2F2F") + ggplot2::geom_point(ggplot2::aes(x = 3.05, y = 0.5),
    shape = 18, size = 3.5, colour = "#2F2F2F") + ggplot2::annotate("text", x = 3.35, y = 0.5, label = "Exact metabolite-name match",
    hjust = 0, vjust = 0.5, size = 5, family = "Arial", colour = "#2F2F2F") + ggplot2::geom_segment(ggplot2::aes(x = 7.15,
    xend = 7.65, y = 0.5, yend = 0.5), linewidth = 0.55, colour = "#5A4E6D") + ggplot2::geom_point(ggplot2::aes(x = 7.4,
    y = 0.5), shape = 16, size = 3.2, colour = "#5A4E6D") + ggplot2::annotate("text", x = 7.8, y = 0.5,
    label = "Higher in PPI users", hjust = 0, vjust = 0.5, size = 5, family = "Arial", colour = "#2F2F2F") +
    ggplot2::geom_segment(ggplot2::aes(x = 10.35, xend = 10.85, y = 0.5, yend = 0.5), linewidth = 0.55,
        colour = "#4C7A9F") + ggplot2::geom_point(ggplot2::aes(x = 10.6, y = 0.5), shape = 16, size = 3.2,
    colour = "#4C7A9F") + ggplot2::annotate("text", x = 11, y = 0.5, label = "Lower in PPI users", hjust = 0,
    vjust = 0.5, size = 5, family = "Arial", colour = "#2F2F2F") + ggplot2::theme_void() + ggplot2::theme(plot.margin = ggplot2::margin(0,
    10, 0, 10))
technical_plot_data <- dplyr::mutate(readr::read_csv(project_input(technical_file), show_col_types = FALSE),
    robust_wbm = as_logical_safe(robust_wbm), cross_platform_exact_name = as_logical_safe(cross_platform_exact_name),
    technical_status = factor(technical_status, levels = c("Other tested reaction", "Robust meta-analysis signal")))
conc <- readr::read_csv(project_input(concordance_file), show_col_types = FALSE)
get_metric <- function(name) {
    hit <- conc$value[conc$metric == name]
    if (!length(hit)) {
        NA_real_
    }
    else {
        as.numeric(hit[[1]])
    }
}
technical_spearman <- get_metric("Spearman correlation")
technical_same_direction <- get_metric("Same-direction fraction")
technical_labels <- dplyr::slice_head(dplyr::arrange(dplyr::filter(technical_plot_data, robust_wbm),
    FDR_meta, dplyr::desc(abs(beta_meta))), n = 3)
p_technical <- ggplot2::ggplot(technical_plot_data, ggplot2::aes(beta_kiel, beta_metacardis)) + ggplot2::geom_hline(yintercept = 0,
    linetype = "dashed", colour = "#9E9E9E", linewidth = 0.4) + ggplot2::geom_vline(xintercept = 0, linetype = "dashed",
    colour = "#9E9E9E", linewidth = 0.4) + ggplot2::geom_abline(slope = 1, intercept = 0, linewidth = 0.4,
    colour = "#BDBDBD") + ggplot2::geom_point(ggplot2::aes(colour = technical_status), size = 2.6, alpha = 0.9) +
    ggplot2::geom_point(data = dplyr::filter(technical_plot_data, robust_wbm, cross_platform_exact_name),
        shape = 23, fill = "#5A4E6D", colour = "#2F2F2F", stroke = 0.7, size = 3.6, show.legend = FALSE) +
    ggrepel::geom_text_repel(data = technical_labels, ggplot2::aes(label = stringr::str_trunc(metabolite_name,
        30)), size = 5.5, force = 5, max.overlaps = Inf, box.padding = 0.35, point.padding = 0.25, segment.colour = "#B0B0B0",
        colour = "#2F2F2F", family = "Arial", show.legend = FALSE) + ggplot2::annotate("text", x = min(technical_plot_data$beta_kiel,
    na.rm = TRUE), y = max(technical_plot_data$beta_metacardis, na.rm = TRUE), hjust = 0, vjust = 1,
    size = 5.5, family = "Arial", label = paste0("Spearman ρ = ", formatC(technical_spearman, digits = 2,
        format = "f"), "\nConcordance = ", formatC(100 * technical_same_direction, digits = 1, format = "f"),
        "%")) + ggplot2::scale_colour_manual(values = pal_sig, breaks = c("Robust meta-analysis signal",
    "Other tested reaction"), labels = c("Meta-analysis significant", "Not significant")) + ggplot2::labs(x = "Kiel blood WBM effect size",
    y = "MetaCardis blood WBM effect size") + main_theme
tag_theme <- ggplot2::theme(plot.tag = ggplot2::element_text(family = "Arial", face = "bold", size = 24),
    plot.tag.position = c(0, 1))
p_technical_tag <- p_technical + ggplot2::labs(tag = "A") + tag_theme
p_wbm_tag <- p_wbm + ggplot2::labs(tag = "B") + tag_theme
p_met_tag <- p_met + ggplot2::labs(tag = "C") + tag_theme
top_full <- patchwork::wrap_elements(full = p_technical_tag, clip = FALSE)
bottom_panels <- p_wbm_tag | p_met_tag
bottom_block <- bottom_panels/forest_legend + patchwork::plot_layout(heights = c(1, 0.08))
final_figure <- top_full/bottom_block + patchwork::plot_layout(heights = c(1, 2.2))
print(final_figure)
pdf_file <- file.path(OUT_DIR, "WBM_Serum_Metabolomics.pdf")
png_file <- file.path(OUT_DIR, "WBM_Serum_Metabolomics.png")
ggplot2::ggsave(project_output(pdf_file), final_figure, width = FIG_W_IN, height = FIG_H_IN, units = "in",
    device = grDevices::cairo_pdf, limitsize = FALSE)
if (requireNamespace("ragg", quietly = TRUE)) {
    ragg::agg_png(project_output(png_file), width = FIG_W_IN, height = FIG_H_IN, units = "in", res = DPI,
        background = "white")
    print(final_figure)
    grDevices::dev.off()
} else {
    ggplot2::ggsave(project_output(png_file), final_figure, width = FIG_W_IN, height = FIG_H_IN, units = "in",
        dpi = DPI, bg = "white", limitsize = FALSE)
}
message("WBM + serum metabolomics figure written to: ", OUT_DIR)
}))
