# SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0
# PCA/Procrustes statistics

source("scripts/helpers/setup.R")
invisible(run_analysis(
  section = "04_mouse_metabolomics",
  directory = "metabolomics",
  script = "04_11_additional_ordination_statistics.R",
  restricted = FALSE,
  code = {
suppressPackageStartupMessages({
    library(dplyr)
    library(tidyr)
    library(readr)
    library(stringr)
    library(tibble)
    library(purrr)
    library(ggplot2)
    library(ggrepel)
    library(lme4)
    library(lmerTest)
    library(broom.mixed)
    library(vegan)
    library(permute)
    library(writexl)
})
set.seed(20260714)
project_dir <- project_path("")
metabolomics_root <- file.path(project_dir, "metabolomics")
out_dir <- file.path(metabolomics_root, "Metabolomics_PCA")
pca_bundle_path <- file.path(out_dir, "PPI_metabolomics_PCA_bundle.rds")
analysis_out_dir <- file.path(out_dir, "PCA_statistical_analysis")
dir.create(project_output(analysis_out_dir, directory = TRUE), recursive = TRUE, showWarnings = FALSE)
if (!project_file_exists(pca_bundle_path)) {
    stop("PCA bundle not found: ", pca_bundle_path, "\nRun 04_10_prepare_metabolomics_ordination.R first.")
}
cumulative_variance_target <- 0.8
maximum_pcs_per_dataset <- 5L
minimum_pcs_per_dataset <- 2L
number_top_loadings_each_direction <- 15L
procrustes_permutations <- 9999L
bundle <- readRDS(project_input(pca_bundle_path))
pca_objects <- bundle$pca_objects
if (length(pca_objects) == 0L) stop("No PCA objects found in bundle.")
select_pcs <- function(obj) {
    v <- obj$variance
    first_at_target <- which(v$cumulative_variance >= cumulative_variance_target)[1]
    if (is.na(first_at_target))
        first_at_target <- nrow(v)
    n_keep <- max(minimum_pcs_per_dataset, min(maximum_pcs_per_dataset, first_at_target))
    v$PC[seq_len(min(n_keep, nrow(v)))]
}
safe_lmer <- function(formula, data) {
    warnings <- character()
    fit <- withCallingHandlers(tryCatch(lmerTest::lmer(formula, data = data, REML = TRUE, control = lme4::lmerControl(optimizer = "bobyqa",
        optCtrl = list(maxfun = 1e+05))), error = function(e) e), warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
    })
    if (inherits(fit, "error")) {
        return(list(fit = NULL, error = conditionMessage(fit), warnings = warnings))
    }
    list(fit = fit, error = NA_character_, warnings = warnings)
}
extract_fixed_coefficients <- function(fit, dataset_id, pc, model_id) {
    out <- broom.mixed::tidy(fit, effects = "fixed", conf.int = TRUE)
    for (nm in c("df", "statistic", "p.value", "conf.low", "conf.high")) {
        if (!nm %in% names(out))
            out[[nm]] <- NA_real_
    }
    out %>% mutate(dataset_id = dataset_id, PC = pc, model_id = model_id, singular = lme4::isSingular(fit,
        tol = 1e-04), n_observations = nobs(fit)) %>% select(dataset_id, PC, model_id, term, estimate,
        std.error, df, statistic, p.value, conf.low, conf.high, singular, n_observations)
}
extract_anova_terms <- function(fit, dataset_id, pc, model_id) {
    a <- as.data.frame(anova(fit, ddf = "Satterthwaite")) %>% rownames_to_column("term") %>% as_tibble()
    names(a) <- names(a) %>% str_replace_all(" ", "_") %>% str_replace_all("\\(", "") %>% str_replace_all("\\)",
        "") %>% str_replace_all(">", "gt") %>% str_replace_all("\\|", "")
    a %>% mutate(dataset_id = dataset_id, PC = pc, model_id = model_id, singular = lme4::isSingular(fit,
        tol = 1e-04), n_observations = nobs(fit))
}
fit_pc_models <- function(obj) {
    dataset_id <- obj$dataset_id
    pcs <- select_pcs(obj)
    dat <- obj$scores
    coefficient_results <- list()
    anova_results <- list()
    diagnostics <- list()
    for (pc in pcs) {
        dat_pc <- dat %>% mutate(PC_score = .data[[pc]])
        model_specs <- list()
        if (dataset_id %in% c("exvivo_nmr_raw", "exvivo_hilic_raw")) {
            model_specs$days1_4_average <- list(data = dat_pc %>% filter(Day_num %in% 1:4) %>% droplevels(),
                formula = PC_score ~ Group + DayF + (1 | Donor))
            model_specs$days0_4_interaction <- list(data = dat_pc %>% filter(Day_num %in% 0:4) %>% droplevels(),
                formula = PC_score ~ Group * DayF + (1 | Donor))
        }
        else if (dataset_id %in% c("exvivo_nmr_delta", "exvivo_hilic_delta")) {
            model_specs$baseline_centred_interaction <- list(data = dat_pc %>% filter(Day_num %in% 1:4) %>%
                droplevels(), formula = PC_score ~ Group * DayF + (1 | Donor))
        }
        else if (dataset_id == "intestinal_lcms_combined") {
            model_specs$paired_compartment <- list(data = dat_pc %>% droplevels(), formula = PC_score ~
                Group * Organ + (1 | Donor) + (1 | Pair_ID))
        }
        else if (dataset_id %in% c("intestinal_lcms_cecum", "intestinal_lcms_feces", "serum_lcms", "cecum_nmr")) {
            model_specs$recipient_clustered <- list(data = dat_pc %>% droplevels(), formula = PC_score ~
                Group + (1 | Donor))
        }
        else {
            next
        }
        for (model_id in names(model_specs)) {
            spec <- model_specs[[model_id]]
            result <- safe_lmer(spec$formula, spec$data)
            diagnostics[[paste(pc, model_id, sep = "__")]] <- tibble(dataset_id = dataset_id, PC = pc,
                model_id = model_id, model_formula = paste(deparse(spec$formula), collapse = " "), fit_failed = is.null(result$fit),
                error = result$error, warnings = paste(unique(result$warnings), collapse = " | "), singular = if (!is.null(result$fit))
                  lme4::isSingular(result$fit, tol = 1e-04)
                else NA, n_observations = nrow(spec$data), n_donors = n_distinct(spec$data$Donor))
            if (!is.null(result$fit)) {
                coefficient_results[[paste(pc, model_id, sep = "__")]] <- extract_fixed_coefficients(result$fit,
                  dataset_id, pc, model_id)
                anova_results[[paste(pc, model_id, sep = "__")]] <- extract_anova_terms(result$fit, dataset_id,
                  pc, model_id)
            }
        }
    }
    list(coefficients = bind_rows(coefficient_results), anova = bind_rows(anova_results), diagnostics = bind_rows(diagnostics),
        selected_pcs = pcs)
}
multivariate_partial_pseudo_F <- function(X, metadata, full_rhs, reduced_rhs) {
    X <- as.matrix(X)
    storage.mode(X) <- "double"
    if (nrow(X) != nrow(metadata)) {
        stop("X and metadata have different numbers of rows.")
    }
    full_formula <- stats::as.formula(paste("~", full_rhs))
    reduced_formula <- stats::as.formula(paste("~", reduced_rhs))
    M_full <- stats::model.matrix(full_formula, data = metadata)
    M_reduced <- stats::model.matrix(reduced_formula, data = metadata)
    rank_full <- qr(M_full)$rank
    rank_reduced <- qr(M_reduced)$rank
    df_effect <- rank_full - rank_reduced
    df_residual <- nrow(X) - rank_full
    if (df_effect <= 0L || df_residual <= 0L) {
        return(NA_real_)
    }
    residuals_full <- qr.resid(qr(M_full), X)
    residuals_reduced <- qr.resid(qr(M_reduced), X)
    rss_full <- sum(residuals_full^2, na.rm = TRUE)
    rss_reduced <- sum(residuals_reduced^2, na.rm = TRUE)
    ss_effect <- rss_reduced - rss_full
    if (is.finite(ss_effect) && ss_effect < 0 && abs(ss_effect) < 1e-10) {
        ss_effect <- 0
    }
    if (!is.finite(ss_effect) || !is.finite(rss_full) || rss_full <= 0) {
        return(NA_real_)
    }
    as.numeric((ss_effect/df_effect)/(rss_full/df_residual))
}
exact_donor_permanova <- function(obj, full_rhs, reduced_rhs, tested_effect, analysis_id) {
    X <- obj$X_processed
    meta <- obj$metadata %>% as.data.frame()
    donors <- sort(unique(as.character(meta$Donor)))
    donor_group <- meta %>% distinct(Donor, Group) %>% mutate(Donor = as.character(Donor), Group = as.character(Group))
    n_ppi <- sum(donor_group$Group == "ppi")
    if (n_ppi <= 0L || n_ppi >= length(donors)) {
        return(tibble(dataset_id = obj$dataset_id, analysis_id = analysis_id, full_model_rhs = full_rhs,
            reduced_model_rhs = reduced_rhs, tested_effect = tested_effect, pseudo_F = NA_real_, exact_p = NA_real_,
            n_assignments = NA_integer_, minimum_possible_p = NA_real_, note = "Invalid donor-group structure"))
    }
    assignments <- combn(donors, n_ppi, simplify = FALSE)
    observed_ppi <- sort(donor_group$Donor[donor_group$Group == "ppi"])
    get_F <- function(ppi_donors) {
        meta_perm <- meta %>% mutate(Group_perm = factor(if_else(as.character(Donor) %in% ppi_donors,
            "ppi", "ctrl"), levels = c("ctrl", "ppi")))
        multivariate_partial_pseudo_F(X = X, metadata = meta_perm, full_rhs = full_rhs, reduced_rhs = reduced_rhs)
    }
    all_F <- vapply(assignments, get_F, numeric(1))
    observed_index <- which(vapply(assignments, function(z) identical(sort(as.character(z)), observed_ppi),
        logical(1)))
    if (length(observed_index) != 1L) {
        stop("Could not uniquely identify the observed donor-group assignment.")
    }
    observed_F <- all_F[[observed_index]]
    valid_F <- is.finite(all_F)
    exact_p <- if (is.finite(observed_F) && any(valid_F)) {
        mean(all_F[valid_F] >= observed_F - 1e-12)
    }
    else {
        NA_real_
    }
    tibble(dataset_id = obj$dataset_id, analysis_id = analysis_id, full_model_rhs = full_rhs, reduced_model_rhs = reduced_rhs,
        tested_effect = tested_effect, pseudo_F = observed_F, exact_p = exact_p, n_assignments = length(assignments),
        n_valid_assignments = sum(valid_F), minimum_possible_p = if (any(valid_F))
            1/sum(valid_F)
        else NA_real_, note = paste0("Exact one-sided donor-label permutation; whole donor clusters retained. ",
            "Partial Euclidean pseudo-F obtained by comparing full and reduced ", "multivariate linear models."))
}
permanova_for_dataset <- function(obj) {
    id <- obj$dataset_id
    out <- list()
    if (id %in% c("exvivo_nmr_raw", "exvivo_hilic_raw")) {
        keep <- obj$metadata$Day_num %in% 1:4
        tmp <- obj
        tmp$X_processed <- obj$X_processed[keep, , drop = FALSE]
        tmp$metadata <- obj$metadata[keep, , drop = FALSE]
        out$group_adjusted_day <- exact_donor_permanova(tmp, full_rhs = "DayF + Group_perm", reduced_rhs = "DayF",
            tested_effect = "Group_perm", analysis_id = "days1_4_group_adjusted_for_day")
        out$interaction <- exact_donor_permanova(tmp, full_rhs = "DayF + Group_perm + DayF:Group_perm",
            reduced_rhs = "DayF + Group_perm", tested_effect = "DayF:Group_perm", analysis_id = "days1_4_group_by_day")
    }
    else if (id %in% c("exvivo_nmr_delta", "exvivo_hilic_delta")) {
        out$group_adjusted_day <- exact_donor_permanova(obj, full_rhs = "DayF + Group_perm", reduced_rhs = "DayF",
            tested_effect = "Group_perm", analysis_id = "baseline_centred_group_adjusted_for_day")
        out$interaction <- exact_donor_permanova(obj, full_rhs = "DayF + Group_perm + DayF:Group_perm",
            reduced_rhs = "DayF + Group_perm", tested_effect = "DayF:Group_perm", analysis_id = "baseline_centred_group_by_day")
    }
    else if (id == "intestinal_lcms_combined") {
        out$group_adjusted_organ <- exact_donor_permanova(obj, full_rhs = "Organ + Group_perm", reduced_rhs = "Organ",
            tested_effect = "Group_perm", analysis_id = "group_adjusted_for_compartment")
        out$interaction <- exact_donor_permanova(obj, full_rhs = "Organ + Group_perm + Organ:Group_perm",
            reduced_rhs = "Organ + Group_perm", tested_effect = "Organ:Group_perm", analysis_id = "group_by_compartment")
    }
    else if (id %in% c("intestinal_lcms_cecum", "intestinal_lcms_feces", "serum_lcms", "cecum_nmr")) {
        out$group <- exact_donor_permanova(obj, full_rhs = "Group_perm", reduced_rhs = "1", tested_effect = "Group_perm",
            analysis_id = "donor_label_group_effect")
    }
    bind_rows(out)
}
exact_two_group_difference <- function(values, donors, group) {
    dat <- tibble(value = values, Donor = donors, Group = group) %>% distinct(Donor, .keep_all = TRUE)
    donor_ids <- sort(dat$Donor)
    n_ppi <- sum(dat$Group == "ppi")
    assignments <- combn(donor_ids, n_ppi, simplify = FALSE)
    observed_ppi <- sort(dat$Donor[dat$Group == "ppi"])
    stat_fun <- function(ppi_ids) {
        mean(dat$value[dat$Donor %in% ppi_ids], na.rm = TRUE) - mean(dat$value[!dat$Donor %in% ppi_ids],
            na.rm = TRUE)
    }
    stats_all <- vapply(assignments, stat_fun, numeric(1))
    obs_idx <- which(vapply(assignments, function(z) identical(sort(z), observed_ppi), logical(1)))
    observed <- stats_all[obs_idx]
    tibble(difference_ppi_minus_ctrl = observed, exact_two_sided_p = mean(abs(stats_all) >= abs(observed) -
        1e-12), n_assignments = length(assignments), minimum_possible_two_sided_p = 2/length(assignments))
}
trajectory_metrics <- function(obj, n_pcs = 2L) {
    pcs <- paste0("PC", seq_len(min(n_pcs, ncol(obj$pca$x))))
    dat <- obj$scores %>% arrange(Donor, Day_num)
    metrics <- dat %>% group_by(Donor, Group) %>% group_modify(function(.x, .y) {
        .x <- .x %>% arrange(Day_num)
        baseline <- .x %>% filter(Day_num == min(Day_num))
        if (nrow(baseline) != 1L)
            return(tibble())
        score_mat <- as.matrix(.x[, pcs, drop = FALSE])
        baseline_vec <- as.numeric(baseline[1, pcs, drop = TRUE])
        distance_from_baseline <- sqrt(rowSums(sweep(score_mat, 2, baseline_vec, "-")^2))
        step_distance <- c(0, sqrt(rowSums((score_mat[-1, , drop = FALSE] - score_mat[-nrow(score_mat),
            , drop = FALSE])^2)))
        tibble(Day_num = .x$Day_num, distance_from_baseline = distance_from_baseline, step_distance = step_distance,
            total_trajectory_length = sum(step_distance), maximum_distance_from_baseline = max(distance_from_baseline),
            final_distance_from_baseline = distance_from_baseline[length(distance_from_baseline)], auc_distance = sum(diff(.x$Day_num) *
                (head(distance_from_baseline, -1) + tail(distance_from_baseline, -1))/2))
    }) %>% ungroup()
    donor_summary <- metrics %>% distinct(Donor, Group, total_trajectory_length, maximum_distance_from_baseline,
        final_distance_from_baseline, auc_distance)
    tests <- map_dfr(c("total_trajectory_length", "maximum_distance_from_baseline", "final_distance_from_baseline",
        "auc_distance"), function(metric) {
        test <- exact_two_group_difference(donor_summary[[metric]], donor_summary$Donor, as.character(donor_summary$Group))
        test %>% mutate(metric = metric, .before = 1)
    })
    list(long = metrics, donor_summary = donor_summary, tests = tests)
}
top_loading_table <- function(obj, pcs) {
    map_dfr(pcs, function(pc) {
        obj$loadings %>% select(Metabolite, loading = all_of(pc)) %>% mutate(dataset_id = obj$dataset_id,
            PC = pc, direction = if_else(loading >= 0, "positive", "negative"), abs_loading = abs(loading)) %>%
            group_by(direction) %>% slice_max(abs_loading, n = number_top_loadings_each_direction, with_ties = FALSE) %>%
            ungroup() %>% arrange(PC, direction, desc(abs_loading))
    })
}
run_exvivo_procrustes <- function(nmr_obj, hilic_obj) {
    nmr_scores <- nmr_obj$scores %>% select(Match_ID, Donor, Group, Day_num, starts_with("PC"))
    hilic_scores <- hilic_obj$scores %>% select(Match_ID, starts_with("PC"))
    matched <- inner_join(nmr_scores, hilic_scores, by = "Match_ID", suffix = c("_NMR", "_HILIC")) %>%
        arrange(Match_ID)
    if (nrow(matched) < 8L) {
        return(list(summary = tibble(note = "Insufficient matched samples")))
    }
    X <- as.matrix(matched[, c("PC1_NMR", "PC2_NMR")])
    Y <- as.matrix(matched[, c("PC1_HILIC", "PC2_HILIC")])
    proc <- vegan::procrustes(X, Y, symmetric = TRUE)
    ctrl <- permute::how(nperm = procrustes_permutations)
    permute::setBlocks(ctrl) <- factor(matched$Donor)
    test <- vegan::protest(X, Y, permutations = ctrl)
    rotated <- as.data.frame(proc$Yrot) %>% setNames(c("HILIC_rotated_1", "HILIC_rotated_2")) %>% bind_cols(matched %>%
        select(Match_ID, Donor, Group, Day_num), .)
    nmr_plot <- as.data.frame(proc$X) %>% setNames(c("NMR_1", "NMR_2")) %>% bind_cols(matched %>% select(Match_ID),
        .)
    plot_data <- left_join(rotated, nmr_plot, by = "Match_ID")
    summary <- tibble(n_matched_samples = nrow(matched), procrustes_correlation = as.numeric(test$t0),
        permutation_p = as.numeric(test$signif), permutations = procrustes_permutations, restriction = "Permutations restricted within donor")
    p <- ggplot(plot_data, aes(colour = Group)) + geom_segment(aes(x = NMR_1, y = NMR_2, xend = HILIC_rotated_1,
        yend = HILIC_rotated_2), linewidth = 0.45, alpha = 0.55, arrow = grid::arrow(length = grid::unit(0.06,
        "inches"))) + geom_point(aes(NMR_1, NMR_2), shape = 16, size = 2.6) + geom_point(aes(HILIC_rotated_1,
        HILIC_rotated_2), shape = 1, size = 2.6) + scale_colour_manual(values = c(ctrl = "#2878B5", ppi = "#C73E3A")) +
        labs(title = "Concordance of baseline-centred ex vivo PCA configurations", subtitle = paste0("Filled: NMR; open: rotated HILIC; Procrustes r = ",
            formatC(test$t0, digits = 3, format = "f"), ", permutation P = ", format.pval(test$signif,
                digits = 3)), x = "Procrustes axis 1", y = "Procrustes axis 2", colour = "Inoculum") +
        theme_classic(base_size = 11)
    list(summary = summary, plot_data = plot_data, plot = p)
}
pc_model_results <- imap(pca_objects, ~fit_pc_models(.x))
coefficients_all <- bind_rows(map(pc_model_results, "coefficients"))
anova_all <- bind_rows(map(pc_model_results, "anova"))
diagnostics_all <- bind_rows(map(pc_model_results, "diagnostics"))
selected_pc_table <- imap_dfr(pc_model_results, ~tibble(dataset_id = .y, PC = .x$selected_pcs))
if (nrow(coefficients_all) > 0L) {
    coefficients_all <- coefficients_all %>% group_by(dataset_id, model_id, term) %>% mutate(FDR_across_tested_PCs = p.adjust(p.value,
        method = "BH")) %>% ungroup()
}
if (nrow(anova_all) > 0L && "PrgtF" %in% names(anova_all)) {
    anova_all <- anova_all %>% group_by(dataset_id, model_id, term) %>% mutate(FDR_across_tested_PCs = p.adjust(PrgtF,
        method = "BH")) %>% ungroup()
}
permanova_all <- map_dfr(pca_objects, permanova_for_dataset)
trajectory_objects <- list()
for (id in c("exvivo_nmr_raw", "exvivo_hilic_raw")) {
    if (id %in% names(pca_objects)) {
        trajectory_objects[[id]] <- trajectory_metrics(pca_objects[[id]], n_pcs = 2L)
    }
}
trajectory_long <- imap_dfr(trajectory_objects, ~.x$long %>% mutate(dataset_id = .y, .before = 1))
trajectory_donor_summary <- imap_dfr(trajectory_objects, ~.x$donor_summary %>% mutate(dataset_id = .y,
    .before = 1))
trajectory_tests <- imap_dfr(trajectory_objects, ~.x$tests %>% mutate(dataset_id = .y, .before = 1))
if (nrow(trajectory_long) > 0L) {
    p_distance <- ggplot(trajectory_long, aes(Day_num, distance_from_baseline, colour = Group, group = Donor)) +
        geom_line(alpha = 0.65) + geom_point(size = 2.2) + facet_wrap(~dataset_id, scales = "free_y") +
        scale_colour_manual(values = c(ctrl = "#2878B5", ppi = "#C73E3A")) + scale_x_continuous(breaks = 0:4) +
        labs(title = "Ex vivo metabolic displacement from day 0", x = "Day", y = "Euclidean distance in PC1-PC2 space",
            colour = "Inoculum") + theme_classic(base_size = 11)
    ggsave(project_output(file.path(analysis_out_dir, "Exvivo_PCA_distance_from_day0.png"), name = "exvivo_pca_distance_from_day0.png"),
        p_distance, width = 8, height = 5.5, dpi = 450, bg = "white")
}
loadings_all <- imap_dfr(pca_objects, function(obj, id) {
    pcs <- pc_model_results[[id]]$selected_pcs
    top_loading_table(obj, pcs)
})
procrustes_result <- NULL
if (all(c("exvivo_nmr_delta", "exvivo_hilic_delta") %in% names(pca_objects))) {
    procrustes_result <- run_exvivo_procrustes(pca_objects$exvivo_nmr_delta, pca_objects$exvivo_hilic_delta)
    if (!is.null(procrustes_result$plot)) {
        ggsave(project_output(file.path(analysis_out_dir, "Exvivo_NMR_HILIC_Procrustes.png"), name = "exvivo_nmr_hilic_procrustes.png"),
            procrustes_result$plot, width = 7, height = 5.8, dpi = 450, bg = "white")
    }
}
write_csv(selected_pc_table, project_output(file.path(analysis_out_dir, "PCA_selected_PCs.csv"), name = "pca_selected_p_cs.csv"))
write_csv(coefficients_all, project_output(file.path(analysis_out_dir, "PCA_PC_mixed_model_coefficients.csv"),
    name = "pca_pc_mixed_model_coefficients.csv"))
write_csv(anova_all, project_output(file.path(analysis_out_dir, "PCA_PC_mixed_model_omnibus_tests.csv"),
    name = "pca_pc_mixed_model_omnibus_tests.csv"))
write_csv(diagnostics_all, project_output(file.path(analysis_out_dir, "PCA_PC_model_diagnostics.csv"), name = "pca_pc_model_diagnostics.csv"))
write_csv(permanova_all, project_output(file.path(analysis_out_dir, "PCA_exact_donor_permanova.csv"), name = "pca_exact_donor_permanova.csv"))
write_csv(loadings_all, project_output(file.path(analysis_out_dir, "PCA_top_loadings.csv"), name = "pca_top_loadings.csv"))
write_csv(trajectory_long, project_output(file.path(analysis_out_dir, "PCA_exvivo_trajectory_distances.csv"),
    name = "pca_exvivo_trajectory_distances.csv"))
write_csv(trajectory_donor_summary, project_output(file.path(analysis_out_dir, "PCA_exvivo_trajectory_donor_summary.csv"),
    name = "pca_exvivo_trajectory_donor_summary.csv"))
write_csv(trajectory_tests, project_output(file.path(analysis_out_dir, "PCA_exvivo_trajectory_exact_tests.csv"),
    name = "pca_exvivo_trajectory_exact_tests.csv"))
procrustes_summary <- if (!is.null(procrustes_result)) {
    procrustes_result$summary
} else {
    tibble(note = "Procrustes analysis not run because matched baseline-centred NMR/HILIC PCA objects were unavailable.")
}
write_csv(procrustes_summary, project_output(file.path(analysis_out_dir, "PCA_exvivo_NMR_HILIC_procrustes.csv"),
    name = "pca_exvivo_nmr_hilic_procrustes.csv"))
writexl::write_xlsx(list(Selected_PCs = selected_pc_table, PC_coefficients = coefficients_all, PC_omnibus_tests = anova_all,
    Model_diagnostics = diagnostics_all, Exact_PERMANOVA = permanova_all, Trajectory_tests = trajectory_tests,
    Trajectory_donor_values = trajectory_donor_summary, Procrustes = procrustes_summary, Top_loadings = loadings_all),
    path = project_output(file.path(analysis_out_dir, "PPI_metabolomics_PCA_statistical_results.xlsx"), name = "ppi_metabolomics_pca_statistical_results.xlsx"))
saveRDS(list(selected_pcs = selected_pc_table, pc_coefficients = coefficients_all, pc_anova = anova_all,
    model_diagnostics = diagnostics_all, permanova = permanova_all, trajectories = trajectory_objects,
    top_loadings = loadings_all, procrustes = procrustes_result), project_output(file.path(analysis_out_dir,
    "PPI_metabolomics_PCA_statistical_results.rds"), name = "ppi_metabolomics_pca_statistical_results.rds"))
message("PCA statistical analysis complete. Outputs written to: ", analysis_out_dir)
}))
