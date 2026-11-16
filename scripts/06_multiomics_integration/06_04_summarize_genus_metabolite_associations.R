# SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0
# Summarize associations
# Run through run.R so input and output paths resolve within the package.
# Inputs: see data/input_manifest.csv and docs/RUN_ORDER.md.
# Writes generated files under results/; bundled inputs are read-only.

suppressPackageStartupMessages({
    library(dplyr)
    library(tidyr)
    library(purrr)
    library(readr)
    library(stringr)
    library(tibble)
    library(ggplot2)
    library(ggrepel)
    library(writexl)
})
message("Loaded 06_04_summarize_genus_metabolite_associations.R")
script_path <- tryCatch(ppi_normalize_path(sys.frame(1)$ofile, winslash = "/", mustWork = FALSE), error = function(e) NA_character_)
script_dir <- if (!is.na(script_path) && nzchar(script_path)) dirname(script_path) else getwd()
config_path <- ppi_code("scripts/06_multiomics_integration/06_01_configure_multiomics.R")
if (!ppi_file_exists(config_path)) config_path <- "06_01_configure_multiomics.R"
ppi_source(config_path)
results_file <- file.path(association_root, "02_association_results", "PPI_multiomics_association_results.rds")
if (!ppi_file_exists(results_file)) {
    stop("Association result bundle not found: ", results_file, "\nRun PPI_multiomics_correlations_02_run_associations_FINAL.R first.")
}
association_bundle <- readRDS(ppi_input(results_file))
all_results <- association_bundle$all_results
significant_results <- association_bundle$significant_results
model_data <- association_bundle$family_model_data
settings <- association_bundle$settings
loo_results <- association_bundle$loo_results
if (!is.null(loo_results) && nrow(loo_results) > 0L) {
    all_results <- all_results %>% left_join(loo_results %>% select(family_id, predictor, outcome, proportion_same_direction, 
        minimum_LOO_beta, maximum_LOO_beta), by = c("family_id", "predictor", "outcome"))
}
if (!"proportion_same_direction" %in% names(all_results)) {
    all_results$proportion_same_direction <- NA_real_
}
if (!"proportion_same_direction" %in% names(significant_results)) {
    significant_results$proportion_same_direction <- NA_real_
}
summary_dir <- file.path(association_root, "03_summaries_and_figures")
heatmap_dir <- file.path(summary_dir, "heatmaps")
scatter_dir <- file.path(summary_dir, "top_scatterplots")
dir.create(ppi_output(heatmap_dir, directory = TRUE), recursive = TRUE, showWarnings = FALSE)
dir.create(ppi_output(scatter_dir, directory = TRUE), recursive = TRUE, showWarnings = FALSE)
if (nrow(significant_results) == 0L) {
    message("No associations passed the configured FDR threshold. Complete family tables and summaries will still be exported.")
}
safe_name <- function(x) {
    x <- str_replace_all(x, "[^A-Za-z0-9._-]+", "_")
    str_sub(x, 1, 180)
}
clean_label <- function(x) {
    x %>% str_replace("^HALLMARK_", "") %>% str_replace("^REACTOME_", "") %>% str_replace_all("_", " ") %>% 
        str_to_sentence()
}
subset_one_family <- function(tbl, requested_family_id, table_name = "results") {
    out <- tbl %>% filter(.data$family_id == .env$requested_family_id)
    if (nrow(out) == 0L) {
        stop("No rows found in ", table_name, " for family: ", requested_family_id)
    }
    observed_ids <- unique(as.character(out$family_id))
    if (length(observed_ids) != 1L || !identical(observed_ids[[1]], as.character(requested_family_id))) {
        stop("Family subsetting failed for ", requested_family_id, ". Observed IDs: ", paste(observed_ids, 
            collapse = "; "))
    }
    out
}
duplicate_significant_keys <- significant_results %>% count(family_id, predictor, outcome, name = "n") %>% 
    filter(.data$n > 1L)
if (nrow(duplicate_significant_keys) > 0L) {
    stop("significant_results contains duplicated family/predictor/outcome keys. ", "Rerun the corrected association script before plotting.")
}
old_heatmaps <- ppi_list_files(heatmap_dir, pattern = "\\.png$", full.names = TRUE)
old_scatterplots <- ppi_list_files(scatter_dir, pattern = "\\.png$", full.names = TRUE)
if (length(old_heatmaps) > 0L) unlink(old_heatmaps)
if (length(old_scatterplots) > 0L) unlink(old_scatterplots)
family_counts <- all_results %>% group_by(layer, family_id, predictor_dataset, outcome_dataset, organ, 
    library) %>% summarise(n_tests = n(), n_rawP0p05 = sum(p_value < 0.05, na.rm = TRUE), n_FDR0p1 = sum(FDR < 
    settings$association_fdr_cutoff, na.rm = TRUE), minimum_FDR = min(FDR, na.rm = TRUE), .groups = "drop") %>% 
    mutate(minimum_FDR = if_else(is.infinite(minimum_FDR), NA_real_, minimum_FDR))
layer_counts <- family_counts %>% group_by(layer) %>% summarise(n_families = n(), n_tests = sum(n_tests), 
    n_FDR0p1 = sum(n_FDR0p1), .groups = "drop")
p_counts <- ggplot(layer_counts, aes(x = reorder(layer, n_FDR0p1), y = n_FDR0p1)) + geom_col(width = 0.72) + 
    coord_flip() + labs(title = "FDR-significant multi-omics associations", subtitle = paste0("Primary family-level BH threshold: FDR < ", 
    settings$association_fdr_cutoff), x = NULL, y = "Number of associations") + theme_classic(base_size = 11)
ggsave(ppi_output(file.path(summary_dir, "significant_association_counts_by_layer.png"), name = "significant_association_counts_by_layer.png"), 
    p_counts, width = 7, height = 4.5, dpi = 450, bg = "white")
mg_all <- all_results %>% filter(layer == "metabolite_genus", outcome_dataset == "mouse_genus") %>% transmute(metabolomics_dataset = predictor_dataset, 
    metabolite = predictor, genus = outcome, beta_metabolite_genus = adjusted_rank_beta, rho_metabolite_genus = raw_spearman_rho, 
    P_metabolite_genus = p_value, FDR_metabolite_genus = FDR, LOO_metabolite_genus = proportion_same_direction)
mp_all <- all_results %>% filter(layer == "metabolite_pathway") %>% transmute(metabolomics_dataset = predictor_dataset, 
    organ, library, metabolite = predictor, pathway = outcome, pathway_GSEA_NES, pathway_GSEA_FDR, beta_metabolite_pathway = adjusted_rank_beta, 
    rho_metabolite_pathway = raw_spearman_rho, P_metabolite_pathway = p_value, FDR_metabolite_pathway = FDR, 
    LOO_metabolite_pathway = proportion_same_direction)
gp_all <- all_results %>% filter(layer == "genus_pathway") %>% transmute(organ, library, genus = predictor, 
    pathway = outcome, beta_genus_pathway = adjusted_rank_beta, rho_genus_pathway = raw_spearman_rho, 
    P_genus_pathway = p_value, FDR_genus_pathway = FDR, LOO_genus_pathway = proportion_same_direction)
triplets_two_edge <- mg_all %>% filter(FDR_metabolite_genus < settings$association_fdr_cutoff) %>% inner_join(mp_all %>% 
    filter(FDR_metabolite_pathway < settings$association_fdr_cutoff), by = c("metabolomics_dataset", 
    "metabolite")) %>% left_join(gp_all, by = c("organ", "library", "genus", "pathway")) %>% mutate(direct_genus_pathway_FDR0p1 = is.finite(FDR_genus_pathway) & 
    FDR_genus_pathway < settings$association_fdr_cutoff, three_edge_FDR0p1 = direct_genus_pathway_FDR0p1, 
    expected_indirect_sign = sign(beta_metabolite_genus * beta_metabolite_pathway), direct_sign = sign(beta_genus_pathway), 
    sign_consistent_three_layer = if_else(is.finite(beta_genus_pathway), expected_indirect_sign == direct_sign, 
        NA), maximum_two_edge_FDR = pmax(FDR_metabolite_genus, FDR_metabolite_pathway), maximum_three_edge_FDR = if_else(is.finite(FDR_genus_pathway), 
        pmax(FDR_metabolite_genus, FDR_metabolite_pathway, FDR_genus_pathway), NA_real_)) %>% distinct(metabolomics_dataset, 
    organ, library, metabolite, genus, pathway, .keep_all = TRUE) %>% arrange(desc(three_edge_FDR0p1), 
    maximum_two_edge_FDR, maximum_three_edge_FDR)
triplets_three_edge <- triplets_two_edge %>% filter(three_edge_FDR0p1)
write_csv(triplets_two_edge, ppi_output(file.path(summary_dir, "supported_triplets_two_metabolite_centered_edges_FDR0p1.csv.gz")))
write_csv(triplets_three_edge, ppi_output(file.path(summary_dir, "supported_triplets_all_three_edges_FDR0p1.csv.gz")))
network_edges <- significant_results %>% mutate(source_id = paste(predictor_type, predictor_dataset, 
    predictor, sep = "::"), target_id = paste(outcome_type, outcome_dataset, outcome, sep = "::"), source_label = predictor, 
    target_label = outcome, weight = abs(adjusted_rank_beta), sign = if_else(adjusted_rank_beta >= 0, 
        "positive", "negative")) %>% transmute(family_id, layer, organ, library, source_id, target_id, 
    source_label, target_label, source_type = predictor_type, target_type = outcome_type, adjusted_rank_beta, 
    raw_spearman_rho, p_value, FDR, weight, sign, n_samples, n_donors, proportion_same_direction)
network_nodes <- bind_rows(network_edges %>% transmute(node_id = source_id, label = source_label, node_type = source_type), 
    network_edges %>% transmute(node_id = target_id, label = target_label, node_type = target_type)) %>% 
    distinct(node_id, .keep_all = TRUE)
write_csv(network_edges, ppi_output(file.path(summary_dir, "network_edges_FDR0p1.csv.gz")))
write_csv(network_nodes, ppi_output(file.path(summary_dir, "network_nodes_FDR0p1.csv"), name = "network_nodes_fdr0p1.csv"))
plot_family_heatmap <- function(family_results, family_id) {
    sig <- family_results %>% filter(FDR < settings$association_fdr_cutoff)
    if (nrow(sig) == 0L) 
        return(NULL)
    predictor_order <- sig %>% group_by(predictor) %>% summarise(best_FDR = min(FDR), max_abs_beta = max(abs(adjusted_rank_beta)), 
        .groups = "drop") %>% arrange(best_FDR, desc(max_abs_beta)) %>% slice_head(n = settings$heatmap_max_predictors) %>% 
        pull(predictor)
    outcome_order <- sig %>% group_by(outcome) %>% summarise(best_FDR = min(FDR), max_abs_beta = max(abs(adjusted_rank_beta)), 
        .groups = "drop") %>% arrange(best_FDR, desc(max_abs_beta)) %>% slice_head(n = settings$heatmap_max_outcomes) %>% 
        pull(outcome)
    dat <- family_results %>% filter(predictor %in% predictor_order, outcome %in% outcome_order) %>% 
        mutate(predictor_label = factor(predictor, levels = rev(predictor_order)), outcome_label = factor(outcome, 
            levels = outcome_order), significance = case_when(FDR < 0.01 ~ "**", FDR < settings$association_fdr_cutoff ~ 
            "*", TRUE ~ ""))
    ggplot(dat, aes(outcome_label, predictor_label, fill = adjusted_rank_beta)) + geom_tile() + geom_text(aes(label = significance), 
        size = 2.5) + scale_fill_gradient2(low = "#2B6CB0", mid = "white", high = "#C53030", midpoint = 0, 
        name = "Adjusted\nrank beta") + scale_x_discrete(labels = clean_label) + scale_y_discrete(labels = clean_label) + 
        labs(title = family_id, subtitle = "* FDR < 0.1; ** FDR < 0.01", x = NULL, y = NULL) + theme_minimal(base_size = 9) + 
        theme(axis.text.x = element_text(angle = 55, hjust = 1, vjust = 1), panel.grid = element_blank(), 
            plot.title = element_text(face = "bold", size = 10))
}
significant_families <- unique(as.character(significant_results$family_id))
for (current_family_id in significant_families) {
    dat <- subset_one_family(all_results, requested_family_id = current_family_id, table_name = "all_results")
    p <- plot_family_heatmap(dat, current_family_id)
    if (!is.null(p)) {
        sig_dat <- dat %>% filter(.data$FDR < settings$association_fdr_cutoff)
        n_pred <- min(n_distinct(sig_dat$predictor), settings$heatmap_max_predictors)
        n_out <- min(n_distinct(sig_dat$outcome), settings$heatmap_max_outcomes)
        ggsave(ppi_output(file.path(heatmap_dir, paste0(safe_name(current_family_id), ".png"))), p, width = max(7, 
            0.28 * n_out + 4), height = max(5, 0.25 * n_pred + 3), dpi = 450, bg = "white", limitsize = FALSE)
    }
}
plot_top_association <- function(hit, family_data) {
    if (is.null(family_data)) 
        return(NULL)
    predictor <- hit$predictor[[1]]
    outcome <- hit$outcome[[1]]
    if (!predictor %in% colnames(family_data$X_original) || !outcome %in% colnames(family_data$Y_original)) 
        return(NULL)
    dat <- family_data$metadata %>% mutate(Predictor = family_data$X_original[, predictor], Outcome = family_data$Y_original[, 
        outcome], Donor = factor(Donor), Group = factor(Group, levels = c("ctrl", "ppi")))
    has_day <- "Day" %in% names(dat) && any(!is.na(dat$Day))
    if (has_day && !"DayF" %in% names(dat)) {
        dat <- dat %>% mutate(DayF = factor(.data$Day))
    }
    p <- ggplot(dat, aes(Predictor, Outcome, colour = Group, label = Donor)) + geom_point(size = 2.8, 
        alpha = 0.9) + geom_smooth(aes(group = Group), method = "lm", formula = y ~ x, se = FALSE, linewidth = 0.65) + 
        ggrepel::geom_text_repel(size = 2.3, show.legend = FALSE, max.overlaps = Inf) + scale_colour_manual(values = c(ctrl = "#2878B5", 
        ppi = "#C73E3A")) + labs(title = paste(clean_label(predictor), "versus", clean_label(outcome)), 
        subtitle = paste0("Adjusted rank beta = ", sprintf("%.2f", hit$adjusted_rank_beta), "; raw Spearman rho = ", 
            sprintf("%.2f", hit$raw_spearman_rho), "; FDR = ", format.pval(hit$FDR, digits = 2), "; n = ", 
            hit$n_samples), x = clean_label(predictor), y = clean_label(outcome), colour = "PPI group") + 
        theme_classic(base_size = 10) + theme(plot.title = element_text(face = "bold"))
    if (has_day) {
        p <- p + facet_wrap(~DayF)
    }
    p
}
for (current_family_id in significant_families) {
    top_hits <- subset_one_family(significant_results, requested_family_id = current_family_id, table_name = "significant_results") %>% 
        arrange(.data$FDR, desc(abs(.data$adjusted_rank_beta))) %>% slice_head(n = settings$scatterplot_top_n_per_family)
    current_model_data <- model_data[[current_family_id]]
    if (is.null(current_model_data)) {
        warning("No model data stored for family ", current_family_id, "; scatterplots for this family were skipped.")
        next
    }
    for (i in seq_len(nrow(top_hits))) {
        hit <- top_hits[i, ]
        p <- plot_top_association(hit, current_model_data)
        if (!is.null(p)) {
            filename <- paste0(safe_name(current_family_id), "__", i, "__", safe_name(hit$predictor), 
                "__", safe_name(hit$outcome), ".png")
            ggsave(ppi_output(file.path(scatter_dir, filename)), p, width = 6.8, height = if (isTRUE(current_model_data$include_day)) 
                6
            else 5.2, dpi = 450, bg = "white")
        }
    }
}
top_associations <- all_results %>% group_by(family_id) %>% arrange(FDR, desc(abs(adjusted_rank_beta)), 
    .by_group = TRUE) %>% slice_head(n = 50) %>% ungroup()
excel_row_limit <- 1000000L
significant_excel <- utils::head(significant_results, excel_row_limit)
three_edge_excel <- utils::head(triplets_three_edge, excel_row_limit)
two_edge_excel <- utils::head(triplets_two_edge, excel_row_limit)
network_nodes_excel <- utils::head(network_nodes, excel_row_limit)
writexl::write_xlsx(list(Layer_counts = layer_counts, Family_counts = family_counts, Significant_associations = significant_excel, 
    Top50_per_family = top_associations, Three_edge_triplets = three_edge_excel, Two_edge_triplets = two_edge_excel, 
    Network_nodes = network_nodes_excel), path = ppi_output(file.path(summary_dir, "PPI_multiomics_integrated_results_review.xlsx"), 
    name = "ppi_multiomics_integrated_results_review.xlsx"))
saveRDS(list(family_counts = family_counts, layer_counts = layer_counts, triplets_two_edge = triplets_two_edge, 
    triplets_three_edge = triplets_three_edge, network_edges = network_edges, network_nodes = network_nodes), 
    ppi_output(file.path(summary_dir, "PPI_multiomics_integrated_summaries.rds"), name = "ppi_multiomics_integrated_summaries.rds"), 
    compress = "xz")
message("Summaries and figures complete. Outputs written to: ", summary_dir)
