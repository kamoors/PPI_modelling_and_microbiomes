# SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0
# Blood-flux meta-analysis

source("scripts/helpers/setup.R")
invisible(run_analysis(
  section = "02_human_metabolism",
  directory = "human_cohorts/blood_flux/metacardis",
  script = "02_02_blood_flux_meta_analysis.R",
  restricted = FALSE,
  code = {
library(dplyr)
library(tidyr)
library(tibble)
library(stringr)
library(readr)
library(ggplot2)
kiel_rds_path <- project_path("human_cohorts/blood_flux/kiel/Blood_FBA_full_workflow_results.rds")
metac_rds_path <- project_path("human_cohorts/blood_flux/metacardis/Blood_FBA_Metacardis/resultFlux_Metacardis_v2/analysis_outputs/MetaCardis_Blood_FBA_v2_full_workflow_results.rds")
dict_path <- project_path("human_cohorts/kiel/recon-store-metabolites-1 (1).tsv")
comparison_out_dir <- project_path("human_cohorts/blood_flux/metacardis/Blood_FBA_Metacardis/resultFlux_Metacardis_v2/analysis_outputs/Kiel_Metacardis_comparison")
dir.create(project_output(comparison_out_dir, directory = TRUE), showWarnings = FALSE, recursive = TRUE)
result_fdr_cutoff <- 0.1
pick_object <- function(x, candidates, required = TRUE, object_description = "object") {
    for (candidate in candidates) {
        if (!is.null(x[[candidate]])) {
            message("Using ", object_description, ": ", candidate)
            return(x[[candidate]])
        }
    }
    if (required) {
        stop("Could not find ", object_description, ". Tried: ", paste(candidates, collapse = ", "))
    }
    NULL
}
extract_column <- function(df, candidates, required = TRUE, default = NA_real_) {
    found <- intersect(candidates, names(df))
    if (length(found) > 0) {
        return(df[[found[1]]])
    }
    if (required) {
        stop("Could not find required column. Tried: ", paste(candidates, collapse = ", "))
    }
    rep(default, nrow(df))
}
canonical_reaction_id <- function(x) {
    x %>% as.character() %>% trimws() %>% str_replace_all("\\s+", "") %>% str_replace("\\[BC\\]$", "[bc]")
}
extract_met_abbr <- function(reaction) {
    reaction %>% canonical_reaction_id() %>% str_remove("^DM_") %>% str_remove("\\[[^]]+\\]$")
}
safe_cor <- function(x, y, method = "pearson") {
    keep <- is.finite(x) & is.finite(y)
    if (sum(keep) < 3) {
        return(NA_real_)
    }
    cor(x[keep], y[keep], method = method)
}
first_nonmissing <- function(x) {
    x <- as.character(x)
    x <- x[!is.na(x) & nzchar(x)]
    if (length(x) == 0) {
        return(NA_character_)
    }
    x[1]
}
standardize_full_results <- function(df, cohort) {
    df <- as_tibble(df)
    if (!"reaction" %in% names(df)) {
        stop("Result table for ", cohort, " has no reaction column.")
    }
    reaction <- as.character(df$reaction)
    reaction_key <- canonical_reaction_id(reaction)
    beta <- extract_column(df, c("beta", "estimate", "Estimate"))
    SE <- extract_column(df, c("SE", "se", "Std.Error", "Std. Error"))
    pvalue <- extract_column(df, c("pvalue", "p", "pval", "P.value"))
    FDR <- extract_column(df, c("FDR", "padj", "qvalue"), required = FALSE)
    n <- extract_column(df, c("n", "N", "sample_size"), required = FALSE)
    if (all(is.na(FDR))) {
        FDR <- p.adjust(as.numeric(pvalue), method = "fdr")
    }
    met_abbr <- if ("met_abbr" %in% names(df)) {
        as.character(df$met_abbr)
    }
    else {
        extract_met_abbr(reaction)
    }
    metabolite_name <- if ("metabolite_name" %in% names(df)) {
        as.character(df$metabolite_name)
    }
    else if ("fullName" %in% names(df)) {
        as.character(df$fullName)
    }
    else {
        NA_character_
    }
    out <- tibble(reaction = reaction, reaction_key = reaction_key, met_abbr = met_abbr, metabolite_name = metabolite_name,
        beta = as.numeric(beta), SE = as.numeric(SE), pvalue = as.numeric(pvalue), FDR = as.numeric(FDR),
        n = as.numeric(n), cohort = cohort) %>% filter(!is.na(reaction_key), reaction_key != "")
    duplicate_reactions <- out %>% count(reaction_key) %>% filter(n > 1)
    if (nrow(duplicate_reactions) > 0) {
        warning(cohort, " has duplicated reaction keys. Keeping the row with lowest FDR/p-value.")
    }
    out %>% arrange(reaction_key, FDR, pvalue) %>% group_by(reaction_key) %>% slice(1) %>% ungroup()
}
two_cohort_fixed_meta <- function(beta_1, se_1, beta_2, se_2) {
    valid <- is.finite(beta_1) & is.finite(se_1) & se_1 > 0 & is.finite(beta_2) & is.finite(se_2) & se_2 >
        0
    w_1 <- rep(NA_real_, length(beta_1))
    w_2 <- rep(NA_real_, length(beta_1))
    w_1[valid] <- 1/se_1[valid]^2
    w_2[valid] <- 1/se_2[valid]^2
    beta_meta <- rep(NA_real_, length(beta_1))
    SE_meta <- rep(NA_real_, length(beta_1))
    beta_meta[valid] <- (w_1[valid] * beta_1[valid] + w_2[valid] * beta_2[valid])/(w_1[valid] + w_2[valid])
    SE_meta[valid] <- sqrt(1/(w_1[valid] + w_2[valid]))
    z_meta <- beta_meta/SE_meta
    p_meta <- 2 * pnorm(-abs(z_meta))
    lower_meta <- beta_meta - 1.96 * SE_meta
    upper_meta <- beta_meta + 1.96 * SE_meta
    Q <- rep(NA_real_, length(beta_1))
    Q[valid] <- w_1[valid] * (beta_1[valid] - beta_meta[valid])^2 + w_2[valid] * (beta_2[valid] - beta_meta[valid])^2
    Q_pvalue <- pchisq(Q, df = 1, lower.tail = FALSE)
    I2 <- rep(NA_real_, length(beta_1))
    I2[valid] <- ifelse(Q[valid] > 0, pmax(0, (Q[valid] - 1)/Q[valid]) * 100, 0)
    tibble(beta_meta = beta_meta, SE_meta = SE_meta, z_meta = z_meta, p_meta = p_meta, lower_meta = lower_meta,
        upper_meta = upper_meta, Q = Q, Q_pvalue = Q_pvalue, I2 = I2)
}
stopifnot(project_file_exists(kiel_rds_path), project_file_exists(metac_rds_path))
kiel_rds <- readRDS(project_input(kiel_rds_path))
metac_rds <- readRDS(project_input(metac_rds_path))
message("Kiel RDS objects:")
print(names(kiel_rds))
message("MetaCardis RDS objects:")
print(names(metac_rds))
kiel_full_raw <- pick_object(kiel_rds, candidates = c("blood_fba_results_annot", "blood_fba_results",
    "annotated_results", "full_results"), object_description = "Kiel full-cohort result table")
metac_full_raw <- pick_object(metac_rds, candidates = c("annotated_results", "full_results", "metac_blood_fba_results_annot",
    "metac_blood_fba_results", "blood_fba_results"), object_description = "MetaCardis full-cohort result table")
kiel_full <- standardize_full_results(kiel_full_raw, cohort = "Kiel")
metac_full <- standardize_full_results(metac_full_raw, cohort = "MetaCardis")
if (project_file_exists(dict_path)) {
    dict_tbl <- readr::read_tsv(project_input(dict_path), show_col_types = FALSE) %>% transmute(met_abbr = as.character(abbreviation),
        metabolite_name_dict = as.character(fullName)) %>% distinct(met_abbr, .keep_all = TRUE)
    kiel_full <- kiel_full %>% left_join(dict_tbl, by = "met_abbr") %>% mutate(metabolite_name = coalesce(metabolite_name,
        metabolite_name_dict, met_abbr)) %>% select(-metabolite_name_dict)
    metac_full <- metac_full %>% left_join(dict_tbl, by = "met_abbr") %>% mutate(metabolite_name = coalesce(metabolite_name,
        metabolite_name_dict, met_abbr)) %>% select(-metabolite_name_dict)
} else {
    kiel_full <- kiel_full %>% mutate(metabolite_name = coalesce(metabolite_name, met_abbr))
    metac_full <- metac_full %>% mutate(metabolite_name = coalesce(metabolite_name, met_abbr))
}
message("Kiel reactions: ", nrow(kiel_full))
message("MetaCardis reactions: ", nrow(metac_full))
exact_shared <- length(intersect(as.character(kiel_full_raw$reaction), as.character(metac_full_raw$reaction)))
canonical_shared <- length(intersect(kiel_full$reaction_key, metac_full$reaction_key))
message("Exact raw reaction-name overlap: ", exact_shared)
message("Canonical reaction-key overlap: ", canonical_shared)
reaction_coverage <- full_join(kiel_full %>% transmute(reaction_key, reaction_kiel = reaction, met_abbr_kiel = met_abbr,
    present_kiel = TRUE), metac_full %>% transmute(reaction_key, reaction_metacardis = reaction, met_abbr_metacardis = met_abbr,
    present_metacardis = TRUE), by = "reaction_key") %>% mutate(present_kiel = replace_na(present_kiel,
    FALSE), present_metacardis = replace_na(present_metacardis, FALSE), coverage = case_when(present_kiel &
    present_metacardis ~ "Both cohorts", present_kiel ~ "Kiel only", present_metacardis ~ "MetaCardis only",
    TRUE ~ "Neither"))
coverage_summary <- reaction_coverage %>% count(coverage, name = "n_reactions")
print(coverage_summary)
write.csv(reaction_coverage, project_output(file.path(comparison_out_dir, "Blood_FBA_reaction_coverage_diagnostics.csv"),
    name = "blood_fba_reaction_coverage_diagnostics.csv"), row.names = FALSE)
write.csv(coverage_summary, project_output(file.path(comparison_out_dir, "Blood_FBA_reaction_coverage_summary.csv"),
    name = "blood_fba_reaction_coverage_summary.csv"), row.names = FALSE)
full_comparison <- inner_join(kiel_full %>% select(reaction_key, reaction_kiel = reaction, met_abbr_kiel = met_abbr,
    metabolite_name_kiel = metabolite_name, beta_kiel = beta, SE_kiel = SE, pvalue_kiel = pvalue, FDR_kiel = FDR,
    n_kiel = n), metac_full %>% select(reaction_key, reaction_metacardis = reaction, met_abbr_metacardis = met_abbr,
    metabolite_name_metacardis = metabolite_name, beta_metacardis = beta, SE_metacardis = SE, pvalue_metacardis = pvalue,
    FDR_metacardis = FDR, n_metacardis = n), by = "reaction_key") %>% mutate(reaction = coalesce(reaction_kiel,
    reaction_metacardis, reaction_key), met_abbr = coalesce(met_abbr_kiel, met_abbr_metacardis), metabolite_name = coalesce(metabolite_name_kiel,
    metabolite_name_metacardis, met_abbr, reaction_key), same_direction = case_when(beta_kiel == 0 |
    beta_metacardis == 0 ~ NA, TRUE ~ sign(beta_kiel) == sign(beta_metacardis)), significance_group = case_when(FDR_kiel <
    result_fdr_cutoff & FDR_metacardis < result_fdr_cutoff ~ "FDR significant in both", FDR_kiel < result_fdr_cutoff ~
    "Kiel only", FDR_metacardis < result_fdr_cutoff ~ "MetaCardis only", TRUE ~ "Not FDR significant"),
    beta_difference = beta_metacardis - beta_kiel, absolute_beta_difference = abs(beta_difference))
message("Shared reactions used for comparison: ", nrow(full_comparison))
comparison_summary <- tibble(n_kiel_reactions = nrow(kiel_full), n_metacardis_reactions = nrow(metac_full),
    n_shared_reactions = nrow(full_comparison), pearson_beta_correlation = safe_cor(full_comparison$beta_kiel,
        full_comparison$beta_metacardis, method = "pearson"), spearman_beta_correlation = safe_cor(full_comparison$beta_kiel,
        full_comparison$beta_metacardis, method = "spearman"), pearson_beta_correlation_if_metacardis_flipped = safe_cor(full_comparison$beta_kiel,
        -full_comparison$beta_metacardis, method = "pearson"), same_direction_fraction = mean(full_comparison$same_direction,
        na.rm = TRUE), kiel_FDR_hits = sum(full_comparison$FDR_kiel < result_fdr_cutoff, na.rm = TRUE),
    metacardis_FDR_hits = sum(full_comparison$FDR_metacardis < result_fdr_cutoff, na.rm = TRUE), FDR_hits_in_both = sum(full_comparison$FDR_kiel <
        result_fdr_cutoff & full_comparison$FDR_metacardis < result_fdr_cutoff, na.rm = TRUE))
print(comparison_summary)
significance_overlap <- full_comparison %>% count(significance_group, name = "n_reactions") %>% arrange(desc(n_reactions))
print(significance_overlap)
replication_summary <- bind_rows(full_comparison %>% filter(FDR_kiel < result_fdr_cutoff) %>% summarise(discovery_cohort = "Kiel",
    n_discovery_FDR = n(), n_same_direction = sum(same_direction, na.rm = TRUE), n_nominal_replication = sum(same_direction &
        pvalue_metacardis < 0.05, na.rm = TRUE), n_FDR_replication = sum(same_direction & FDR_metacardis <
        result_fdr_cutoff, na.rm = TRUE)), full_comparison %>% filter(FDR_metacardis < result_fdr_cutoff) %>%
    summarise(discovery_cohort = "MetaCardis", n_discovery_FDR = n(), n_same_direction = sum(same_direction,
        na.rm = TRUE), n_nominal_replication = sum(same_direction & pvalue_kiel < 0.05, na.rm = TRUE),
        n_FDR_replication = sum(same_direction & FDR_kiel < result_fdr_cutoff, na.rm = TRUE))) %>% mutate(same_direction_fraction = n_same_direction/n_discovery_FDR,
    nominal_replication_fraction = n_nominal_replication/n_discovery_FDR, FDR_replication_fraction = n_FDR_replication/n_discovery_FDR)
print(replication_summary)
meta_stats <- two_cohort_fixed_meta(beta_1 = full_comparison$beta_kiel, se_1 = full_comparison$SE_kiel,
    beta_2 = full_comparison$beta_metacardis, se_2 = full_comparison$SE_metacardis)
full_meta_results <- bind_cols(full_comparison, meta_stats) %>% mutate(FDR_meta = p.adjust(p_meta, method = "fdr"),
    meta_significant = FDR_meta < result_fdr_cutoff, consistent_and_meta_significant = same_direction &
        meta_significant) %>% arrange(FDR_meta, p_meta)
meta_summary <- comparison_summary %>% mutate(meta_analysis_FDR_hits = sum(full_meta_results$FDR_meta <
    result_fdr_cutoff, na.rm = TRUE), meta_hits_same_direction = sum(full_meta_results$FDR_meta < result_fdr_cutoff &
    full_meta_results$same_direction, na.rm = TRUE), meta_hits_heterogeneity_p_lt_0.05 = sum(full_meta_results$FDR_meta <
    result_fdr_cutoff & full_meta_results$Q_pvalue < 0.05, na.rm = TRUE))
print(meta_summary)
message("Top fixed-effect meta-analysis results:")
print(full_meta_results %>% select(reaction, met_abbr, metabolite_name, beta_kiel, beta_metacardis, beta_meta,
    SE_meta, p_meta, FDR_meta, same_direction, Q_pvalue, I2) %>% head(30))
top_disagreements <- full_meta_results %>% arrange(desc(absolute_beta_difference)) %>% select(reaction,
    met_abbr, metabolite_name, beta_kiel, beta_metacardis, beta_difference, FDR_kiel, FDR_metacardis,
    same_direction, Q_pvalue, I2) %>% head(50)
print(top_disagreements)
p_scatter <- ggplot(full_meta_results, aes(x = beta_kiel, y = beta_metacardis, color = significance_group)) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.4) + geom_vline(xintercept = 0, linetype = "dashed",
    linewidth = 0.4) + geom_abline(slope = 1, intercept = 0, linetype = "dotted", linewidth = 0.5) +
    geom_point(alpha = 0.8, size = 2) + theme_bw() + labs(title = "Kiel vs MetaCardis Blood FBA PPI effects",
    subtitle = paste0("Pearson r = ", round(comparison_summary$pearson_beta_correlation, 3), "; Spearman rho = ",
        round(comparison_summary$spearman_beta_correlation, 3), "; same direction = ", round(100 * comparison_summary$same_direction_fraction,
            1), "%"), x = "Kiel PPI effect", y = "MetaCardis PPI effect", color = NULL)
print(p_scatter)
p_disagreement <- ggplot(full_meta_results, aes(x = beta_meta, y = I2)) + geom_point(alpha = 0.75) +
    theme_bw() + labs(title = "Kiel-MetaCardis heterogeneity", x = "Fixed-effect meta-analysis beta",
    y = "I2 (%)")
print(p_disagreement)
valid_primary <- full_meta_results %>% filter(same_direction, FDR_meta < result_fdr_cutoff, Q_pvalue >=
    0.05, I2 <= 25) %>% arrange(FDR_meta, p_meta)
valid_primary %>% select(reaction, met_abbr, metabolite_name, beta_kiel, beta_metacardis, beta_meta,
    FDR_kiel, FDR_metacardis, FDR_meta, Q_pvalue, I2)
library(data.table)
library(dplyr)
library(stringr)
library(ggplot2)
map_file <- project_path("figures/source_data/Cohort_Microbiome_Function/VMH_HMDB_map.csv")
hmdb_file <- project_path("figures/source_data/Cohort_Microbiome_Function/hmdb_class_annotations_collapsed.tsv")
hmdb_enrichment_min_class_size <- 5
hmdb_enrichment_min_foreground <- 2
hmdb_enrichment_fdr_cutoff <- 0.1
pad_hmdb <- function(x) {
    x <- trimws(as.character(x))
    x[x == ""] <- NA_character_
    out <- x
    need_pad <- !is.na(out) & grepl("^HMDB[0-9]+$", out)
    nums <- suppressWarnings(as.integer(sub("^HMDB", "", out[need_pad])))
    out[need_pad] <- paste0("HMDB", sprintf("%07d", nums))
    out
}
extract_blood_fba_met_abbr <- function(reaction) {
    reaction %>% as.character() %>% str_trim() %>% str_remove("^DM_") %>% str_remove("\\[[^]]+\\]$") %>%
        str_remove("\\([^)]*\\)$") %>% str_remove("_(bc|e|c)$")
}
read_vmh_hmdb_map <- function(map_file) {
    map <- data.table::fread(project_input(map_file), showProgress = FALSE)
    if (!all(c("met_id", "HMDB") %in% names(map))) {
        map <- data.table::fread(project_input(map_file), drop = 1, header = TRUE, showProgress = FALSE)
    }
    if (!all(c("met_id", "HMDB") %in% names(map))) {
        stop("Could not read VMH-HMDB map. Expected columns met_id and HMDB.")
    }
    map[, `:=`(met_id, trimws(as.character(met_id)))]
    map[, `:=`(HMDB, trimws(as.character(HMDB)))]
    map[HMDB == "", `:=`(HMDB, NA_character_)]
    map <- map[!is.na(HMDB) & HMDB != "HMDB0000216"]
    map[, `:=`(HMDB, pad_hmdb(HMDB))]
    unique(map[, .(met_abbr = met_id, HMDB)])
}
build_hmdb_class_membership <- function(hmdb_file, remove_root_terms = TRUE) {
    hmdb_annotation <- data.table::fread(project_input(hmdb_file), showProgress = FALSE)
    hmdb_col <- intersect(c("hmdb", "HMDB"), names(hmdb_annotation))[1]
    if (is.na(hmdb_col)) {
        stop("Could not identify HMDB column in HMDB annotation file.")
    }
    class_cols <- intersect(c("kingdom", "super_class", "superclass", "class", "sub_class", "subclass",
        "direct_parent", "class_label"), names(hmdb_annotation))
    class_cols <- setdiff(class_cols, hmdb_col)
    if (length(class_cols) == 0) {
        stop("Could not identify class columns in HMDB annotation file.")
    }
    hmdb_long_raw <- data.table::melt(hmdb_annotation, id.vars = hmdb_col, measure.vars = class_cols,
        variable.name = "source_column", value.name = "raw_class_label")
    data.table::setnames(hmdb_long_raw, hmdb_col, "HMDB")
    hmdb_long_raw <- hmdb_long_raw[!is.na(HMDB) & HMDB != "" & !is.na(raw_class_label) & raw_class_label !=
        ""]
    hmdb_long <- hmdb_long_raw[, {
        parts <- unlist(strsplit(as.character(raw_class_label), ";", fixed = TRUE))
        parts <- trimws(parts)
        parts <- parts[parts != ""]
        if (length(parts) == 0) {
            NULL
        }
        else {
            data.table::data.table(class_depth = seq_along(parts), class_label = parts)
        }
    }, by = .(HMDB, source_column)]
    hmdb_long[, `:=`(HMDB, pad_hmdb(HMDB))]
    hmdb_long[, `:=`(class_level, paste0("level_", class_depth))]
    if (remove_root_terms) {
        hmdb_long <- hmdb_long[!class_label %in% c("Organic compounds", "Inorganic compounds")]
    }
    unique(hmdb_long[, .(HMDB, class_level, class_label)])
}
perform_valid_primary_hmdb_ora <- function(full_meta_results, valid_primary, map_file, hmdb_file, min_class_size = 5,
    min_foreground = 2, fdr_cutoff = 0.1) {
    model <- data.table::as.data.table(full_meta_results)
    foreground <- data.table::as.data.table(valid_primary)
    if (!"reaction_key" %in% names(model)) {
        model[, `:=`(reaction_key, as.character(reaction))]
    }
    if (!"reaction_key" %in% names(foreground)) {
        foreground[, `:=`(reaction_key, as.character(reaction))]
    }
    if (!"met_abbr" %in% names(model)) {
        model[, `:=`(met_abbr, extract_blood_fba_met_abbr(reaction))]
    }
    model[, `:=`(met_abbr, trimws(as.character(met_abbr)))]
    model[, `:=`(reaction_key, trimws(as.character(reaction_key)))]
    foreground[, `:=`(reaction_key, trimws(as.character(reaction_key)))]
    map <- read_vmh_hmdb_map(map_file)
    hmdb_classes <- build_hmdb_class_membership(hmdb_file)
    model_hmdb <- merge(model[, .(reaction_key, reaction, met_abbr, metabolite_name, beta_kiel, beta_metacardis,
        beta_meta, FDR_meta, Q_pvalue, I2, same_direction)], map, by = "met_abbr", all.x = FALSE, allow.cartesian = TRUE)
    model_hmdb <- unique(model_hmdb, by = c("reaction_key", "HMDB"))
    model_class <- merge(model_hmdb, hmdb_classes, by = "HMDB", all.x = FALSE, allow.cartesian = TRUE)
    unit_class <- unique(model_class[, .(reaction_key, class_level, class_label)])
    background_units <- unique(unit_class$reaction_key)
    foreground_units <- intersect(unique(foreground$reaction_key), background_units)
    n_background <- length(background_units)
    n_foreground <- length(foreground_units)
    message("Model reactions in full_meta_results: ", data.table::uniqueN(model$reaction_key))
    message("Model reactions with HMDB mapping: ", data.table::uniqueN(model_hmdb$reaction_key))
    message("Background reactions with HMDB class annotation: ", n_background)
    message("Valid-primary foreground reactions: ", data.table::uniqueN(foreground$reaction_key))
    message("Foreground reactions with HMDB class annotation: ", n_foreground)
    class_counts <- unit_class[, .(n_background_class = data.table::uniqueN(reaction_key), n_foreground_class = data.table::uniqueN(reaction_key[reaction_key %in%
        foreground_units])), by = .(class_level, class_label)]
    class_counts <- class_counts[n_background_class >= min_class_size & n_foreground_class >= min_foreground]
    enrichment <- rbindlist(lapply(seq_len(nrow(class_counts)), function(i) {
        this_class_level <- class_counts$class_level[i]
        this_class <- class_counts$class_label[i]
        class_units <- unique(unit_class[class_level == this_class_level & class_label == this_class,
            reaction_key])
        a <- length(intersect(foreground_units, class_units))
        b <- n_foreground - a
        c <- length(setdiff(class_units, foreground_units))
        d <- n_background - n_foreground - c
        fisher_mat <- matrix(c(a, b, c, d), nrow = 2, byrow = TRUE, dimnames = list(foreground = c("in_class",
            "not_in_class"), background = c("in_class", "not_in_class")))
        fisher_res <- fisher.test(fisher_mat, alternative = "greater")
        data.table(class_level = this_class_level, class_label = this_class, n_foreground_class = a,
            n_foreground_total = n_foreground, n_background_class = length(class_units), n_background_total = n_background,
            fold_enrichment = (a/n_foreground)/(length(class_units)/n_background), odds_ratio = unname(fisher_res$estimate),
            pvalue = fisher_res$p.value)
    }))
    if (nrow(enrichment) == 0) {
        warning("No HMDB classes passed min_class_size / min_foreground filters.")
        return(list(enrichment_all = data.table(), enrichment_sig = data.table(), model_hmdb = model_hmdb,
            model_class = model_class, unit_class = unit_class))
    }
    enrichment[, `:=`(FDR, p.adjust(pvalue, method = "fdr"))]
    member_table <- model_class[reaction_key %in% foreground_units][, .(foreground_reactions = paste(unique(reaction),
        collapse = "; "), foreground_metabolites = paste(unique(metabolite_name), collapse = "; "), mean_foreground_beta_meta = mean(beta_meta,
        na.rm = TRUE)), by = .(class_level, class_label)]
    enrichment <- merge(enrichment, member_table, by = c("class_level", "class_label"), all.x = TRUE)
    data.table::setorder(enrichment, FDR, pvalue)
    enrichment_sig <- enrichment[FDR < fdr_cutoff]
    list(enrichment_all = enrichment, enrichment_sig = enrichment_sig, model_hmdb = model_hmdb, model_class = model_class,
        unit_class = unit_class)
}
hmdb_valid_primary_enrich <- perform_valid_primary_hmdb_ora(full_meta_results = full_meta_results, valid_primary = valid_primary,
    map_file = map_file, hmdb_file = hmdb_file, min_class_size = hmdb_enrichment_min_class_size, min_foreground = hmdb_enrichment_min_foreground,
    fdr_cutoff = hmdb_enrichment_fdr_cutoff)
hmdb_valid_primary_enrich$enrichment_all %>% arrange(FDR, pvalue) %>% dplyr::select(class_level, class_label,
    n_foreground_class, n_foreground_total, n_background_class, n_background_total, fold_enrichment,
    odds_ratio, pvalue, FDR, foreground_metabolites) %>% tibble::as_tibble() %>% print(n = 50)
hmdb_valid_primary_enrich$enrichment_sig %>% arrange(FDR, pvalue) %>% select(class_level, class_label,
    n_foreground_class, n_background_class, fold_enrichment, pvalue, FDR, foreground_metabolites) %>%
    tibble::as_tibble() %>% print(n = 50)
write.csv(kiel_full, project_output(file.path(comparison_out_dir, "Blood_FBA_Kiel_standardized_full_results.csv"),
    name = "blood_fba_kiel_standardized_full_results.csv"), row.names = FALSE)
write.csv(metac_full, project_output(file.path(comparison_out_dir, "Blood_FBA_MetaCardis_v2_standardized_full_results.csv"),
    name = "blood_fba_meta_cardis_standardized_full_results.csv"), row.names = FALSE)
write.csv(full_comparison, project_output(file.path(comparison_out_dir, "Blood_FBA_Kiel_vs_MetaCardis_v2_full_comparison.csv"),
    name = "blood_fba_kiel_vs_meta_cardis_full_comparison.csv"), row.names = FALSE)
write.csv(full_meta_results, project_output(file.path(comparison_out_dir, "Blood_FBA_Kiel_MetaCardis_v2_fixed_effect_meta_analysis.csv"),
    name = "blood_fba_kiel_meta_cardis_fixed_effect_meta_analysis.csv"), row.names = FALSE)
write.csv(comparison_summary, project_output(file.path(comparison_out_dir, "Blood_FBA_Kiel_MetaCardis_v2_comparison_summary.csv"),
    name = "blood_fba_kiel_meta_cardis_comparison_summary.csv"), row.names = FALSE)
write.csv(meta_summary, project_output(file.path(comparison_out_dir, "Blood_FBA_Kiel_MetaCardis_v2_meta_summary.csv"),
    name = "blood_fba_kiel_meta_cardis_meta_summary.csv"), row.names = FALSE)
write.csv(significance_overlap, project_output(file.path(comparison_out_dir, "Blood_FBA_Kiel_MetaCardis_v2_significance_overlap.csv"),
    name = "blood_fba_kiel_meta_cardis_significance_overlap.csv"), row.names = FALSE)
write.csv(replication_summary, project_output(file.path(comparison_out_dir, "Blood_FBA_Kiel_MetaCardis_v2_replication_summary.csv"),
    name = "blood_fba_kiel_meta_cardis_replication_summary.csv"), row.names = FALSE)
write.csv(top_disagreements, project_output(file.path(comparison_out_dir, "Blood_FBA_Kiel_MetaCardis_v2_top_disagreements.csv"),
    name = "blood_fba_kiel_meta_cardis_top_disagreements.csv"), row.names = FALSE)
write.csv(valid_primary, project_output(file.path(comparison_out_dir, "Blood_FBA_Kiel_MetaCardis_v2_filtered_meta_results.csv"),
    name = "blood_fba_kiel_meta_cardis_filtered_meta_results.csv"), row.names = FALSE)
ggsave(project_output(file.path(comparison_out_dir, "Blood_FBA_Kiel_vs_MetaCardis_v2_effect_scatter.pdf"),
    name = "blood_fba_kiel_vs_meta_cardis_effect_scatter.pdf"), p_scatter, width = 8, height = 7)
ggsave(project_output(file.path(comparison_out_dir, "Blood_FBA_Kiel_vs_MetaCardis_v2_heterogeneity.pdf"),
    name = "blood_fba_kiel_vs_meta_cardis_heterogeneity.pdf"), p_disagreement, width = 7, height = 6)
saveRDS(list(kiel_full = kiel_full, metac_full = metac_full, reaction_coverage = reaction_coverage, coverage_summary = coverage_summary,
    full_comparison = full_comparison, full_meta_results = full_meta_results, comparison_summary = comparison_summary,
    meta_summary = meta_summary, significance_overlap = significance_overlap, replication_summary = replication_summary,
    top_disagreements = top_disagreements), project_output(file.path(comparison_out_dir, "Blood_FBA_Kiel_vs_MetaCardis_v2_comparison_results.rds"),
    name = "blood_fba_kiel_vs_meta_cardis_comparison_results.rds"))
message("Comparison complete.")
message("Outputs saved to: ", comparison_out_dir)
}))
