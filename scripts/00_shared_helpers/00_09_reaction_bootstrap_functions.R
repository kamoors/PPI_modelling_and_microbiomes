# SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0
# Additional bootstrap helpers
# Run through run.R so input and output paths resolve within the package.
# Inputs: see data/input_manifest.csv and docs/RUN_ORDER.md.
# Writes generated files under results/; bundled inputs are read-only.

hc3_boot_one <- function(test_df, formula, ppi_term, R = 1000) {
    n <- nrow(test_df)
    beta_star <- rep(NA_real_, R)
    for (b in seq_len(R)) {
        idx <- sample.int(n, replace = TRUE)
        dboot <- test_df[idx, , drop = FALSE]
        fit <- tryCatch(lm(formula, data = dboot), error = function(e) NULL)
        if (is.null(fit)) 
            next
        ct <- tryCatch(lmtest::coeftest(fit, vcov. = sandwich::vcovHC(fit, type = "HC3")), error = function(e) NULL)
        if (is.null(ct)) 
            next
        if (!(ppi_term %in% rownames(ct))) 
            next
        beta_star[b] <- ct[ppi_term, "Estimate"]
    }
    beta_star <- beta_star[is.finite(beta_star)]
    res <- data.frame(Estimate = NA_real_, p.value = NA_real_)
    rownames(res) <- ppi_term
    if (length(beta_star) == 0) {
        attr(res, "boot_beta") <- numeric(0)
        return(res)
    }
    res[ppi_term, "Estimate"] <- median(beta_star, na.rm = TRUE)
    res[ppi_term, "p.value"] <- 2 * min(mean(beta_star <= 0, na.rm = TRUE), mean(beta_star >= 0, na.rm = TRUE))
    attr(res, "boot_beta") <- beta_star
    res
}
bootstrap_subset <- function(data_in, test_df, fml_boot, ppi_term, R = 1000) {
    out <- vector("list", ncol(data_in))
    for (i in seq_len(ncol(data_in))) {
        test_df$met <- data_in[, i]
        out[[i]] <- hc3_boot_one(test_df, fml_boot, ppi_term = ppi_term, R = R)
    }
    names(out) <- colnames(data_in)
    out
}
build_tables <- function(react_ab_glm, sig_react_ab_boot, ppi_term) {
    hc3_tbl <- tibble(rxn = names(react_ab_glm), hc3_coef = sapply(react_ab_glm, function(ct) ct[ppi_term, 
        "Estimate"]), hc3_p = sapply(react_ab_glm, function(ct) ct[ppi_term, "Pr(>|t|)"])) %>% mutate(hc3_fdr = p.adjust(hc3_p, 
        method = "fdr"))
    boot_tbl <- tibble(rxn = names(sig_react_ab_boot), boot_coef_median = sapply(sig_react_ab_boot, function(x) x[ppi_term, 
        "Estimate"]), boot_p = sapply(sig_react_ab_boot, function(x) x[ppi_term, "p.value"]), boot_ci_low = sapply(sig_react_ab_boot, 
        function(x) {
            bb <- attr(x, "boot_beta")
            if (is.null(bb) || length(bb) == 0) 
                NA_real_
            else as.numeric(quantile(bb, 0.025, na.rm = TRUE))
        }), boot_ci_high = sapply(sig_react_ab_boot, function(x) {
        bb <- attr(x, "boot_beta")
        if (is.null(bb) || length(bb) == 0) 
            NA_real_
        else as.numeric(quantile(bb, 0.975, na.rm = TRUE))
    }))
    res_clean <- hc3_tbl %>% left_join(boot_tbl, by = "rxn") %>% mutate(hc3_sig_fdr_0.1 = hc3_fdr < 0.1, 
        boot_sig_p_0.05 = !is.na(boot_p) & boot_p < 0.05, sig_both = hc3_sig_fdr_0.1 & boot_sig_p_0.05) %>% 
        arrange(desc(sig_both), hc3_fdr, hc3_p)
    res_clean_report <- res_clean %>% dplyr::select(rxn, hc3_coef, hc3_p, hc3_fdr, boot_coef_median, 
        boot_ci_low, boot_ci_high, boot_p, sig_both)
    list(hc3_tbl = hc3_tbl, boot_tbl = boot_tbl, res_clean_report = res_clean_report)
}
