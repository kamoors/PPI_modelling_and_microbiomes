# SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0
# Mixed-model helper called by existing serum loader/curation

run_mixedmeta_cross_sectional_single_compartment <- function(metadata, features, group_levels = c("ctrl",
    "ppi"), block_factor = NULL, log_transform = TRUE, tic_normalize = FALSE, allow_metabolite_duplicates = FALSE) {
    require(lmerTest)
    require(broom.mixed)
    require(dplyr)
    require(tidyr)
    require(stringr)
    require(tibble)
    require(purrr)
    meta <- metadata
    if (is.null(block_factor) || !block_factor %in% names(meta)) {
        meta$Donor <- str_extract(meta$Sample_Name, "(?:ctrl|ppi)[0-9]+")
        block_factor <- "Donor"
    }
    meta <- meta %>% mutate(Separation = factor(Separation, levels = group_levels), Label = factor(.data[[block_factor]]))
    feat_mat <- as.matrix(features)
    stopifnot(nrow(feat_mat) == nrow(meta))
    grp_vec <- factor(meta$Separation, levels = group_levels)
    ok_col <- apply(feat_mat, 2, function(x) {
        v_ok <- stats::sd(x, na.rm = TRUE) > 0
        n_per_grp <- tapply(x, grp_vec, function(v) sum(is.finite(v)))
        n_ok <- all(n_per_grp >= 2, na.rm = TRUE)
        isTRUE(v_ok) && isTRUE(n_ok)
    })
    feat_mat <- feat_mat[, ok_col, drop = FALSE]
    if (tic_normalize) {
        sample_totals <- rowSums(feat_mat, na.rm = TRUE)
        mean_total <- mean(sample_totals, na.rm = TRUE)
        feat_mat <- sweep(feat_mat, 1, sample_totals, "/") * mean_total
    }
    X <- feat_mat
    if (log_transform) {
        pc <- apply(feat_mat, 2, function(col) {
            m <- suppressWarnings(min(col[col > 0], na.rm = TRUE))
            if (!is.finite(m))
                m <- 1e-12
            m/2
        })
        X <- log2(sweep(feat_mat, 2, pc, "+"))
    }
    metabolite_name <- colnames(X)
    if (is.null(metabolite_name))
        metabolite_name <- paste0("feat", seq_len(ncol(X)))
    if (!allow_metabolite_duplicates && any(duplicated(metabolite_name))) {
        dups <- unique(metabolite_name[duplicated(metabolite_name)])
        stop("Duplicate metabolite names detected: ", paste(head(dups, 10), collapse = ", "), if (length(dups) >
            10)
            " ...", "\nSet allow_metabolite_duplicates = TRUE to aggregate duplicates (hybrid), ", "or provide unique metabolite names.")
    }
    feature_id <- make.unique(metabolite_name, sep = "__feat")
    colnames(X) <- feature_id
    feature_map <- tibble(feature_id = feature_id, metabolite_name = metabolite_name)
    X_df <- as.data.frame(X)
    X_df$Sample_Name <- meta$Sample_Name
    long_df <- X_df %>% relocate(Sample_Name) %>% pivot_longer(-Sample_Name, names_to = "feature_id",
        values_to = "log2int") %>% left_join(meta, by = "Sample_Name") %>% left_join(feature_map, by = "feature_id") %>%
        mutate(grp = factor(Separation, levels = group_levels), label = factor(Label))
    form <- as.formula("log2int ~ grp + (1|label)")
    grp_term <- paste0("grp", group_levels[2])
    fit_one_feature <- function(df_feat) {
        m <- suppressMessages(lmerTest::lmer(form, data = df_feat, REML = FALSE, control = lme4::lmerControl(check.conv.singular = "ignore")))
        singular_fit <- lme4::isSingular(m, tol = 1e-04)
        vc <- as.data.frame(lme4::VarCorr(m))
        donor_variance_values <- vc$vcov[vc$grp == "label"]
        donor_variance <- if (length(donor_variance_values) > 0) {
            donor_variance_values[[1]]
        }
        else {
            NA_real_
        }
        residual_variance <- stats::sigma(m)^2
        donor_icc <- if (is.finite(donor_variance) && is.finite(residual_variance) && donor_variance +
            residual_variance > 0) {
            donor_variance/(donor_variance + residual_variance)
        }
        else {
            NA_real_
        }
        tt <- broom.mixed::tidy(m, effects = "fixed", conf.int = TRUE, conf.method = "Wald")
        out <- tt %>% dplyr::filter(term == grp_term) %>% dplyr::transmute(estimate, std.error, df, statistic,
            p.value, conf.low, conf.high)
        if (nrow(out) == 0) {
            return(NULL)
        }
        out %>% dplyr::mutate(log2_FC = estimate, fold_change = 2^estimate, fc_ci_low = 2^conf.low, fc_ci_high = 2^conf.high,
            n_obs = stats::nobs(m), singular_fit = singular_fit, donor_variance = donor_variance, residual_variance = residual_variance,
            donor_icc = donor_icc)
    }
    safe_fit <- purrr::possibly(fit_one_feature, otherwise = NULL, quiet = TRUE)
    feature_results <- long_df %>% group_by(feature_id, metabolite_name) %>% group_modify(~safe_fit(.x)) %>%
        ungroup() %>% filter(!is.na(estimate)) %>% mutate(p_adj_feature = p.adjust(p.value, method = "BH")) %>%
        arrange(p_adj_feature)
    combine_metabolite <- function(df) {
        k <- nrow(df)
        w <- 1/(df$std.error^2)
        w[!is.finite(w)] <- NA
        if (all(is.na(w)))
            return(NULL)
        w[is.na(w)] <- 0
        wsum <- sum(w)
        beta_hat <- sum(w * df$estimate)/wsum
        se_hat <- sqrt(1/wsum)
        z_meta <- beta_hat/se_hat
        p_meta <- 2 * pnorm(abs(z_meta), lower.tail = FALSE)
        ci_low <- beta_hat - 1.96 * se_hat
        ci_high <- beta_hat + 1.96 * se_hat
        fisher_stat <- -2 * sum(log(df$p.value))
        fisher_df <- 2 * k
        p_fisher <- pchisq(fisher_stat, df = fisher_df, lower.tail = FALSE)
        signs <- sign(df$estimate)
        n_pos <- sum(signs > 0, na.rm = TRUE)
        n_neg <- sum(signs < 0, na.rm = TRUE)
        agree_prop <- if (k > 0)
            max(n_pos, n_neg)/k
        else NA_real_
        mixed_dir <- n_pos > 0 & n_neg > 0
        tibble(n_features = k, meta_log2FC = beta_hat, meta_SE = se_hat, meta_z = z_meta, meta_p = p_meta,
            meta_ci_low = ci_low, meta_ci_high = ci_high, meta_FC = 2^beta_hat, meta_FC_ci_low = 2^ci_low,
            meta_FC_ci_high = 2^ci_high, fisher_p = p_fisher, n_pos = n_pos, n_neg = n_neg, agree_prop = agree_prop,
            mixed_direction = mixed_dir)
    }
    metabolite_results <- feature_results %>% group_by(metabolite_name) %>% group_modify(~combine_metabolite(.x)) %>%
        ungroup() %>% mutate(p_adj_meta = p.adjust(meta_p, method = "BH"), p_adj_fisher = p.adjust(fisher_p,
        method = "BH")) %>% arrange(p_adj_meta)
    list(feature_level = feature_results, metabolite_level = metabolite_results)
}
