# SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0
# Figure 1: cohort taxonomy, oxygen-associated reactions and secretion.

source("scripts/helpers/setup.R")
invisible(run_analysis(
  section = "01_human_microbiome",
  directory = "figures/scripts",
  script = "01_09_plot_human_microbiome.R",
  restricted = FALSE,
  code = {
required_packages <- c("dplyr", "data.table", "readr", "tidyr", "ggplot2", "ggrepel", "forcats", "stringr",
    "clusterProfiler", "patchwork", "ggsignif", "grid", "scales")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) {
    stop("Install required packages first: ", paste(missing_packages, collapse = ", "))
}
project_check_plot_versions(required_packages)
FIG_ROOT <- project_path("figures")
DATA_DIR <- file.path(FIG_ROOT, "source_data", "Cohort_Microbiome_Function")
OUT_DIR <- file.path(FIG_ROOT, "final")
dir.create(project_output(DATA_DIR, directory = TRUE), recursive = TRUE, showWarnings = FALSE)
dir.create(project_output(OUT_DIR, directory = TRUE), recursive = TRUE, showWarnings = FALSE)
REFRESH_SOURCE_DATA <- FALSE
stage_input <- function(filename, candidates) {
    target <- file.path(DATA_DIR, filename)
    if (REFRESH_SOURCE_DATA || !project_file_exists(target)) {
        hits <- candidates[project_file_exists(candidates)]
        if (!length(hits)) {
            stop("Could not stage ", filename, ".\nChecked:\n", paste(candidates, collapse = "\n"))
        }
        ok <- file.copy(project_input(hits[[1]]), project_output(target), overwrite = TRUE)
        if (!ok)
            stop("Failed to copy: ", hits[[1]])
        message("Staged: ", filename)
    }
    target
}
PPI_ROOT <- project_path("")
aldex_file <- stage_input("Kiel_species_ALDEx2_FULL_PPI_clean_results.csv", c(file.path(PPI_ROOT, "Results",
    "Aldex2", "Kiel_species_ALDEx2_FULL_PPI_clean_results.csv")))
maaslin_file <- stage_input("MetaCardis_MaAsLin2_all_results.tsv", c(file.path(PPI_ROOT, "Results", "Kiel_Metacardis",
    "Maaslin2_FULL", "all_results.tsv")))
mg_file <- stage_input("Kiel_MetaCardis_predicted_secretion_meta.csv", c(file.path(PPI_ROOT, "Results",
    "Kiel_Metacardis", "Kiel_Metacardis_mg_sec_comparison_with_meta_focus_newPPI.csv")))
hmdb_map_file <- stage_input("VMH_HMDB_map.csv", c(file.path(PPI_ROOT, "Results", "Kiel_Metacardis",
    "figs_PPI_Paper", "figs", "data", "VMH_HMDB_map.csv")))
hmdb_class_file <- stage_input("hmdb_class_annotations_collapsed.tsv", c(file.path(PPI_ROOT, "Results",
    "Kiel_Metacardis", "figs_PPI_Paper", "figs", "data", "hmdb_class_annotations_collapsed.tsv")))
anaerobe_file <- stage_input("anaerobe_plot.RDS", c(file.path(PPI_ROOT, "Results", "Kiel_Metacardis",
    "figs", "anaerobe_plot.RDS")))
reaction_file <- stage_input("reaction_abundance_meta_analysis.csv", c(file.path(PPI_ROOT, "Results",
    "Kiel_Metacardis", "Focus_metacardis_react_ab_meta_analysis_focusNewPPI.csv")))
reaction_map_file <- stage_input("unique_panModel_reactions.csv", c(file.path(PPI_ROOT, "Bram", "transfer_130444_files_9d5afaf0",
    "mgpipRes", "totalresults", "unique_panModel_reactions.csv")))
manifest <- data.frame(file = basename(c(aldex_file, maaslin_file, mg_file, hmdb_map_file, hmdb_class_file,
    anaerobe_file, reaction_file, reaction_map_file)), centralized_path = c(aldex_file, maaslin_file,
    mg_file, hmdb_map_file, hmdb_class_file, anaerobe_file, reaction_file, reaction_map_file))
readr::write_csv(manifest, project_output(file.path(DATA_DIR, "Cohort_Microbiome_Function_source_manifest.csv"),
    name = "cohort_microbiome_function_source_manifest.csv"))
sig_threshold <- 0.1
BASE_SIZE <- 18
FIG_W_IN <- 18
FIG_H_IN <- 18
DPI <- 600
pal_direction <- c(`Higher in PPI users` = "#5A4E6D", `Lower in PPI users` = "#4C7A9F")
pal_taxa <- c(`Higher in PPI users` = "#5A4E6D", `Lower in PPI users` = "#4C7A9F", `Not significant` = "#D3D3D3")
pal_agree <- c(Concordant = "#5A4E6D", Discordant = "#B85C5C")
main_theme <- ggplot2::theme_bw(base_size = BASE_SIZE, base_family = "Arial") + ggplot2::theme(panel.grid = ggplot2::element_blank(),
    plot.title = ggplot2::element_blank(), axis.text = ggplot2::element_text(size = 17, colour = "#2F2F2F"),
    axis.title = ggplot2::element_text(size = 18, colour = "#2F2F2F"), legend.position = "bottom", legend.title = ggplot2::element_blank(),
    legend.text = ggplot2::element_text(size = 16), strip.text = ggplot2::element_text(size = 16), plot.margin = ggplot2::margin(4,
        4, 2, 4))
clean_species <- function(x) {
    x <- as.character(x)
    x <- stringr::str_replace_all(x, "\\.", "_")
    x <- stringr::str_replace_all(x, "[^A-Za-z0-9_]", "_")
    x <- stringr::str_replace_all(x, "_+", "_")
    x <- stringr::str_replace(x, "^_", "")
    x <- stringr::str_replace(x, "_$", "")
    vapply(strsplit(x, "_"), function(parts) {
        if (length(parts) >= 2)
            paste(parts[1], parts[2], sep = "_")
        else NA_character_
    }, character(1))
}
pretty_species <- function(x) stringr::str_replace_all(x, "_", " ")
pad_hmdb <- function(x) {
    x <- trimws(x)
    num <- sub("^HMDB", "", x, ignore.case = TRUE)
    num <- sub("^0+", "", num)
    num <- sprintf("%07d", as.integer(num))
    paste0("HMDB", num)
}
build_class_membership <- function(hmdb_annotation) {
    hmdb_annotation[, `:=`(groups_list, strsplit(class, ";", fixed = TRUE))]
    hmdb_long <- hmdb_annotation[, .(class_label = trimws(unlist(groups_list))[-1]), by = hmdb]
    hmdb_long <- hmdb_long[!is.na(class_label) & class_label != ""]
    unique(hmdb_long, by = c("hmdb", "class_label"))
}
read_csv_clean_names <- function(path) {
    x <- utils::read.csv(project_input(path), check.names = FALSE, stringsAsFactors = FALSE)
    nms <- names(x)
    bad <- is.na(nms) | trimws(nms) == ""
    if (any(bad)) {
        nms[bad] <- paste0("unnamed_", seq_len(sum(bad)))
    }
    names(x) <- make.unique(nms, sep = "_dup")
    x <- x[, !grepl("^unnamed_[0-9]+$", names(x)), drop = FALSE]
    x
}
ald <- read_csv_clean_names(aldex_file)
ald_clean <- dplyr::filter(dplyr::mutate(ald, qval = p.adjust(pval, method = "BH"), species = clean_species(feature),
    label = pretty_species(species), direction = dplyr::case_when(qval < sig_threshold & coef > 0 ~ "Higher in PPI users",
        qval < sig_threshold & coef < 0 ~ "Lower in PPI users", TRUE ~ "Not significant"), neglog10_fdr = -log10(pmax(qval,
        .Machine$double.xmin))), !is.na(species), is.finite(coef), is.finite(qval))
maas <- readr::read_tsv(project_input(maaslin_file), show_col_types = FALSE)
maas_clean <- dplyr::filter(dplyr::transmute(dplyr::filter(maas, metadata == "Proton.pump.inhibitor",
    value == "Yes"), species = clean_species(feature), coef = coef, pval = pval, qval = qval, label = pretty_species(species),
    direction = dplyr::case_when(qval < sig_threshold & coef > 0 ~ "Higher in PPI users", qval < sig_threshold &
        coef < 0 ~ "Lower in PPI users", TRUE ~ "Not significant"), neglog10_fdr = -log10(pmax(qval,
        .Machine$double.xmin))), !is.na(species), is.finite(coef), is.finite(qval))
make_taxa_volcano <- function(df, label_n = 8) {
    label_df <- dplyr::slice_head(dplyr::arrange(dplyr::filter(df, qval < sig_threshold), qval, dplyr::desc(abs(coef))),
        n = label_n)
    ggplot2::ggplot(df, ggplot2::aes(coef, neglog10_fdr)) + ggplot2::geom_vline(xintercept = 0, linetype = "dashed",
        colour = "#9E9E9E", linewidth = 0.4) + ggplot2::geom_hline(yintercept = -log10(sig_threshold),
        linetype = "dashed", colour = "#9E9E9E", linewidth = 0.4) + ggplot2::geom_point(ggplot2::aes(colour = direction),
        size = 2.1, alpha = 0.85) + ggrepel::geom_text_repel(data = label_df, ggplot2::aes(label = stringr::str_trunc(label,
        34)), size = 5.5, max.overlaps = Inf, box.padding = 0.35, point.padding = 0.2, segment.colour = "#B0B0B0",
        show.legend = FALSE) + ggplot2::scale_colour_manual(values = pal_taxa) + ggplot2::labs(x = "Effect size",
        y = expression(-log[10]("FDR"))) + main_theme
}
p_taxa_kiel <- make_taxa_volcano(ald_clean)
p_taxa_meta <- make_taxa_volcano(maas_clean)
anaerobe_plot <- readRDS(project_input(anaerobe_file))
ana_df <- as.data.frame(anaerobe_plot$data)
ana_df$class_label <- factor(ana_df$class_label, levels = c("aerobe\n(n=247)", "facultative\n(n=419)",
    "anaerobe\n(n=1412)", "mixed\n(n=4279)"), labels = c("Aerobe\n(n=247)", "Facultative\n(n=419)", "Anaerobe\n(n=1412)",
    "Mixed\n(n=4279)"))
pair_tests <- list(c("Aerobe\n(n=247)", "Facultative\n(n=419)"), c("Aerobe\n(n=247)", "Anaerobe\n(n=1412)"),
    c("Facultative\n(n=419)", "Anaerobe\n(n=1412)"), c("Aerobe\n(n=247)", "Mixed\n(n=4279)"), c("Facultative\n(n=419)",
        "Mixed\n(n=4279)"), c("Anaerobe\n(n=1412)", "Mixed\n(n=4279)"))
test_results <- dplyr::mutate(dplyr::bind_rows(lapply(pair_tests, function(pair) {
    g1 <- ana_df$beta_meta[ana_df$class_label == pair[1]]
    g2 <- ana_df$beta_meta[ana_df$class_label == pair[2]]
    data.frame(p = stats::wilcox.test(g1, g2)$p.value)
})), fdr = p.adjust(p, method = "fdr"), label = paste0("FDR=", format(fdr, scientific = TRUE, digits = 2)))
p_anaerobe <- ggplot2::ggplot(ana_df, ggplot2::aes(class_label, beta_meta)) + ggplot2::geom_hline(yintercept = 0,
    linetype = "dashed", colour = "#8C8C8C", linewidth = 0.4) + ggplot2::geom_boxplot(outlier.shape = NA,
    width = 0.6, fill = "#D9D9D9", colour = "#3A3A3A", linewidth = 0.4) + ggplot2::geom_jitter(width = 0.15,
    size = 0.8, alpha = 0.25, colour = "#5A4E6D") + ggsignif::geom_signif(comparisons = pair_tests, annotations = test_results$label,
    y_position = c(3.65, 4.05, 4.45, 4.85, 5.25, 5.65), tip_length = 0.01, textsize = 4, vjust = 0.25) +
    ggplot2::coord_cartesian(ylim = c(-1.2, 5.8), clip = "off") + ggplot2::labs(x = NULL, y = "Meta-analysis effect size") +
    main_theme + ggplot2::theme(legend.position = "none", axis.text.x = ggplot2::element_text(angle = 30,
    hjust = 1, vjust = 1, size = 14, lineheight = 0.9), plot.margin = ggplot2::margin(7, 7, 12, 7))
react_df <- read_csv_clean_names(reaction_file)
panModel_rxns <- read_csv_clean_names(reaction_map_file)
rxn_map <- dplyr::distinct(dplyr::filter(tidyr::drop_na(data.frame(subsystem = panModel_rxns$Subsystem,
    rxn_name = panModel_rxns$ReactionID)), subsystem != ""))
react_ab_data <- dplyr::left_join(dplyr::select(react_df, -dplyr::any_of("subsystem")), rxn_map, by = "rxn_name")
sig_df <- dplyr::filter(react_ab_data, fdr_meta < sig_threshold)
term2gene <- dplyr::distinct(dplyr::filter(react_ab_data, !is.na(subsystem), subsystem != ""), subsystem,
    rxn_name)
enrich_pos <- clusterProfiler::enricher(gene = sig_df$rxn_name[sig_df$beta_meta > 0], TERM2GENE = term2gene,
    maxGSSize = 500, minGSSize = 1, qvalueCutoff = 0.1, universe = react_ab_data$rxn_name)@result
enrich_neg <- clusterProfiler::enricher(gene = sig_df$rxn_name[sig_df$beta_meta < 0], TERM2GENE = term2gene,
    maxGSSize = 500, minGSSize = 1, qvalueCutoff = 0.1, universe = react_ab_data$rxn_name)@result
enrich_pos$up_down <- "Increased"
enrich_neg$up_down <- "Decreased"
enriched_terms <- dplyr::mutate(dplyr::distinct(dplyr::bind_rows(dplyr::filter(enrich_pos, qvalue < sig_threshold),
    dplyr::filter(enrich_neg, qvalue < sig_threshold)), Description, up_down, qvalue), enrich_direction = up_down)
box_df <- dplyr::mutate(dplyr::filter(dplyr::left_join(dplyr::filter(react_ab_data, !is.na(subsystem),
    subsystem != "", subsystem %in% enriched_terms$Description, is.finite(beta_meta), fdr_meta < sig_threshold),
    enriched_terms, by = c(subsystem = "Description")), (up_down == "Increased" & beta_meta > 0) | (up_down ==
    "Decreased" & beta_meta < 0)), subsystem_plot = forcats::fct_reorder(subsystem, beta_meta, .fun = median),
    direction = ifelse(beta_meta > 0, "Higher in PPI users", "Lower in PPI users"), enrich_direction = factor(enrich_direction,
        levels = c("Decreased", "Increased")))
p_react_enrich <- ggplot2::ggplot(box_df, ggplot2::aes(beta_meta, subsystem_plot)) + ggplot2::geom_vline(xintercept = 0,
    linetype = "dashed", linewidth = 0.4, colour = "#8C8C8C") + ggplot2::geom_boxplot(outlier.shape = NA,
    width = 0.55, fill = "white", colour = "#3A3A3A", linewidth = 0.35) + ggplot2::geom_jitter(ggplot2::aes(colour = direction),
    height = 0.13, width = 0, size = 1.2, alpha = 0.9) + ggplot2::facet_grid(enrich_direction ~ ., scales = "free_y",
    space = "free_y") + ggplot2::scale_colour_manual(values = pal_direction) + ggplot2::labs(x = "Meta-analysis effect size",
    y = NULL) + main_theme + ggplot2::theme(axis.text.y = ggplot2::element_text(size = 13.5), strip.background = ggplot2::element_rect(fill = "#D9D9D9",
    colour = "#3A3A3A"), strip.text.y = ggplot2::element_text(face = "bold", angle = 270, size = 13.5))
react_sig_scatter <- dplyr::mutate(dplyr::filter(react_df, fdr_meta < sig_threshold, is.finite(hc3_coef.f),
    is.finite(hc3_coef.m)), agreement = ifelse(sign(hc3_coef.f) == sign(hc3_coef.m), "Concordant", "Discordant"))
p_react_scatter <- ggplot2::ggplot(react_sig_scatter, ggplot2::aes(hc3_coef.f, hc3_coef.m)) + ggplot2::geom_hline(yintercept = 0,
    linetype = "dashed", colour = "#9E9E9E", linewidth = 0.4) + ggplot2::geom_vline(xintercept = 0, linetype = "dashed",
    colour = "#9E9E9E", linewidth = 0.4) + ggplot2::geom_abline(slope = 1, intercept = 0, linewidth = 0.4,
    colour = "#BDBDBD") + ggplot2::geom_point(ggplot2::aes(colour = agreement), size = 1.6, alpha = 0.9) +
    ggplot2::scale_colour_manual(values = pal_agree) + ggplot2::labs(x = "Kiel effect size", y = "MetaCardis effect size") +
    main_theme
mg_df <- read_csv_clean_names(mg_file)
mg_plot_df <- dplyr::mutate(mg_df, label = ifelse(!is.na(fullName) & fullName != "", fullName, metabolite),
    direction = dplyr::case_when(fdr_meta < sig_threshold & beta_meta > 0 ~ "Higher in PPI users", fdr_meta <
        sig_threshold & beta_meta < 0 ~ "Lower in PPI users", TRUE ~ "Not significant"))
top_mg <- dplyr::mutate(dplyr::slice_head(dplyr::arrange(dplyr::filter(mg_plot_df, fdr_meta < sig_threshold),
    fdr_meta, dplyr::desc(abs(beta_meta))), n = 30), label_short = stringr::str_trunc(label, 52), label_plot = forcats::fct_reorder(label_short,
    beta_meta))
p_mg_top30 <- ggplot2::ggplot(top_mg, ggplot2::aes(beta_meta, label_plot, fill = direction)) + ggplot2::geom_col(width = 0.75) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", colour = "#3A3A3A", linewidth = 0.4) + ggplot2::scale_fill_manual(values = pal_direction) +
    ggplot2::labs(x = "Meta-analysis effect size", y = NULL) + main_theme + ggplot2::theme(axis.text.y = ggplot2::element_text(size = 13))
performSecretionEnrichmentPlot <- function(excr_file, map_file, hmdb_file, min_class_size = 5, fdr_cutoff = 0.1) {
    excr <- data.table::fread(project_input(excr_file))
    map <- data.table::fread(project_input(map_file), drop = 1, header = TRUE)
    map <- map[HMDB != "HMDB0000216"]
    hmdb_annotation <- data.table::fread(project_input(hmdb_file))
    excr[, `:=`(metabr, trimws(metabr))]
    map[, `:=`(met_id, trimws(met_id))]
    map[, `:=`(HMDB, trimws(HMDB))]
    map <- unique(map, by = c("met_id", "HMDB"))
    map[HMDB == "", `:=`(HMDB, NA_character_)]
    data.table::setnames(map, "met_id", "metabr")
    excr <- merge(excr, map, by = "metabr", all.x = TRUE, allow.cartesian = TRUE)
    excr <- excr[!is.na(HMDB)]
    excr[, `:=`(HMDB, pad_hmdb(HMDB))]
    global_median <- median(excr$beta_meta, na.rm = TRUE)
    hmdb_long <- build_class_membership(hmdb_annotation)
    hmdb_long[, `:=`(hmdb, pad_hmdb(hmdb))]
    hmdb_long <- hmdb_long[hmdb %in% excr[, HMDB]]
    class_sizes <- table(hmdb_long$class_label)
    class_sizes <- class_sizes[class_sizes >= min_class_size]
    res <- data.table::rbindlist(lapply(names(class_sizes), function(testClass) {
        hmdbs <- unique(hmdb_long[class_label == testClass, hmdb])
        thisClass <- excr[HMDB %in% hmdbs, beta_meta]
        otherClass <- excr[!(HMDB %in% hmdbs), beta_meta]
        data.table::data.table(class = testClass, meanThis = mean(thisClass, na.rm = TRUE), meanOther = mean(otherClass,
            na.rm = TRUE), pval = stats::wilcox.test(thisClass, otherClass)$p.value)
    }))
    res[, `:=`(padj, p.adjust(pval, method = "fdr"))]
    res <- res[padj <= fdr_cutoff]
    res[, `:=`(effect, meanThis - meanOther)]
    data.table::setorder(res, effect)
    plot_dt <- data.table::rbindlist(lapply(res$class, function(testClass) {
        hmdbs <- unique(hmdb_long[class_label == testClass, hmdb])
        vals <- excr[HMDB %in% hmdbs, .(HMDB, beta_meta)]
        vals[, `:=`(class, testClass)]
        vals
    }), use.names = TRUE, fill = TRUE)
    plot_dt[, `:=`(class, factor(class, levels = res$class))]
    plot_dt <- merge(plot_dt, res[, .(class, padj)], by = "class", all.x = TRUE)
    ggplot2::ggplot(plot_dt, ggplot2::aes(class, beta_meta, fill = padj)) + ggplot2::geom_boxplot(outlier.shape = 16,
        outlier.size = 0.8, outlier.alpha = 0.7, colour = "#3A3A3A", linewidth = 0.35) + ggplot2::geom_hline(yintercept = 0,
        linetype = "dashed", colour = "#8C8C8C", linewidth = 0.4) + ggplot2::geom_hline(yintercept = global_median,
        linetype = "dashed", colour = "#B85C5C", linewidth = 0.4) + ggplot2::scale_fill_gradient(low = "#5A4E6D",
        high = "#D9D9D9", name = "FDR", guide = ggplot2::guide_colorbar(reverse = TRUE)) + ggplot2::labs(x = NULL,
        y = "Meta-analysis effect size") + ggplot2::coord_flip(clip = "off") + main_theme + ggplot2::theme(axis.text.y = ggplot2::element_text(size = 13),
        legend.position = "bottom", legend.direction = "horizontal", legend.justification = "center",
        legend.background = ggplot2::element_blank(), legend.key.width = grid::unit(1.5, "cm"), legend.key.height = grid::unit(0.3,
            "cm"), legend.title = ggplot2::element_text(size = 13), legend.text = ggplot2::element_text(size = 14))
}
p_mg_enrich <- performSecretionEnrichmentPlot(mg_file, hmdb_map_file, hmdb_class_file, min_class_size = 5,
    fdr_cutoff = sig_threshold)
mg_sig_scatter <- dplyr::mutate(dplyr::filter(mg_df, fdr_meta < sig_threshold, is.finite(hc3_coef.focus),
    is.finite(hc3_coef.metac)), agreement = ifelse(sign(hc3_coef.focus) == sign(hc3_coef.metac), "Concordant",
    "Discordant"), label = ifelse(!is.na(fullName) & fullName != "", fullName, metabolite))
mg_rho <- stats::cor(mg_sig_scatter$hc3_coef.focus, mg_sig_scatter$hc3_coef.metac, method = "spearman",
    use = "complete.obs")
mg_conc <- mean(mg_sig_scatter$agreement == "Concordant", na.rm = TRUE) * 100
mg_labels <- dplyr::slice_head(dplyr::arrange(mg_sig_scatter, dplyr::desc(abs(beta_meta))), n = 3)
p_mg_scatter <- ggplot2::ggplot(mg_sig_scatter, ggplot2::aes(hc3_coef.focus, hc3_coef.metac)) + ggplot2::geom_hline(yintercept = 0,
    linetype = "dashed", colour = "#9E9E9E", linewidth = 0.4) + ggplot2::geom_vline(xintercept = 0, linetype = "dashed",
    colour = "#9E9E9E", linewidth = 0.4) + ggplot2::geom_abline(slope = 1, intercept = 0, linewidth = 0.4,
    colour = "#BDBDBD") + ggplot2::geom_point(ggplot2::aes(colour = agreement), size = 1.7, alpha = 0.9) +
    ggrepel::geom_text_repel(data = mg_labels, ggplot2::aes(label = stringr::str_trunc(label, 30)), size = 5,
        max.overlaps = Inf, segment.colour = "#B0B0B0", show.legend = FALSE) + ggplot2::annotate("text",
    x = min(mg_sig_scatter$hc3_coef.focus, na.rm = TRUE), y = max(mg_sig_scatter$hc3_coef.metac, na.rm = TRUE) -
        0.15 * diff(range(mg_sig_scatter$hc3_coef.metac, na.rm = TRUE)), hjust = 0, vjust = 1, size = 5,
    label = paste0("Spearman ρ = ", round(mg_rho, 2), "\nConcordance = ", round(mg_conc, 1), "%")) +
    ggplot2::scale_colour_manual(values = pal_agree) + ggplot2::labs(x = "Kiel effect size", y = "MetaCardis effect size") +
    main_theme
top_row <- patchwork::free(p_taxa_kiel, side = "l", type = "space") + p_taxa_meta + patchwork::plot_layout(widths = c(1,
    1))
middle_row <- patchwork::free(p_anaerobe, side = "lb", type = "space") + p_react_enrich + p_react_scatter +
    patchwork::plot_layout(widths = c(1, 1.18, 1))
bottom_row <- p_mg_top30 + p_mg_enrich + p_mg_scatter + patchwork::plot_layout(widths = c(1.08, 1.15,
    1))
final_figure <- (top_row/middle_row/bottom_row + patchwork::plot_layout(heights = c(0.9, 1.08, 1.08)) +
    patchwork::plot_annotation(tag_levels = "A", theme = ggplot2::theme(plot.tag = ggplot2::element_text(family = "Arial",
        face = "bold", size = 24)))) & ggplot2::theme(legend.box.spacing = grid::unit(1, "mm"), legend.margin = ggplot2::margin(t = 0,
    r = 0, b = 0, l = 0), legend.spacing.y = grid::unit(0.5, "mm"), legend.spacing.x = grid::unit(1.5,
    "mm"), plot.margin = ggplot2::margin(t = 4, r = 4, b = 2, l = 4))
pdf_file <- file.path(OUT_DIR, "Cohort_Microbiome_Function.pdf")
png_file <- file.path(OUT_DIR, "Cohort_Microbiome_Function.png")
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
message("Figure 1 written to: ", OUT_DIR)
}))
