# SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0
# Additional FULL-bootstrap and BASE-model exports

hc3_p_full <- sapply(mg_sec_glm_full, function(ct) if (ppi_term %in% rownames(ct)) ct[ppi_term, "Pr(>|t|)"] else NA_real_)
hc3_fdr_full <- p.adjust(hc3_p_full, method = "fdr")
rxns_to_boot_full <- names(hc3_fdr_full)[is.finite(hc3_fdr_full) & hc3_fdr_full < 0.05]
data_in_full_sig <- data_in[, rxns_to_boot_full, drop = FALSE]
sig_boot_full <- bootstrap_subset(data_in_full_sig, test_df, fml_boot_full, ppi_term, R = 1000, transform_y = y_transform)
tables_full <- build_tables(mg_sec_glm_full, sig_boot_full, ppi_term, hc3_fdr_cut = 0.2, boot_p_cut = 0.05)
mg_net_sec_res_list_full <- list(Net_secretion_HC3_full = mg_sec_glm_full, Sig_net_secretion_boot_full_output = sig_boot_full,
    tables_full = tables_full)
saveRDS(mg_net_sec_res_list_full, project_output("MGPIPE_Net_secretion_HC3_FULL_sigGLMpadj_02_BOOT_results.rds",
    name = "kiel_metabolite_secretion_full_bootstrap.rds"))
write.csv(tables_full$res_clean_report, file = project_output(project_path("Results/mgpipe_net_secretion_HC3_FULL_clean_results.csv"),
    name = "kiel_metabolite_secretion_full_bootstrap_summary.csv"), row.names = FALSE)
write.csv(tables_full$res_clean_report %>% filter(sig_both), file = project_output(project_path("Results/mgpipe_net_secretion_HC3_FULL_sig_both.csv"),
    name = "kiel_metabolite_secretion_full_bootstrap_supported.csv"), row.names = FALSE)
mg_sec_glm_base <- fit_all_features_hc3(data_in, test_df, fml_base, transform_y = y_transform)
mgsec_results_preboot_base <- tibble(metabolite = names(mg_sec_glm_base), hc3_coef = sapply(mg_sec_glm_base,
    function(ct) {
        if (is.null(ct))
            return(NA_real_)
        if (ppi_term %in% rownames(ct))
            ct[ppi_term, "Estimate"]
        else NA_real_
    }), hc3_p = sapply(mg_sec_glm_base, function(ct) {
    if (is.null(ct))
        return(NA_real_)
    if (ppi_term %in% rownames(ct))
        ct[ppi_term, "Pr(>|t|)"]
    else NA_real_
})) %>% mutate(hc3_fdr = p.adjust(hc3_p, method = "fdr")) %>% arrange(hc3_fdr, hc3_p)
hc3_p_base <- sapply(mg_sec_glm_base, function(ct) if (ppi_term %in% rownames(ct)) ct[ppi_term, "Pr(>|t|)"] else NA_real_)
hc3_fdr_base <- p.adjust(hc3_p_base, method = "fdr")
rxns_to_boot_base <- names(hc3_fdr_base)[is.finite(hc3_fdr_base) & hc3_fdr_base < 0.1]
data_in_base_sig <- data_in[, rxns_to_boot_base, drop = FALSE]
sig_boot_base <- bootstrap_subset(data_in_base_sig, test_df, fml_boot_base, ppi_term, R = 1000, transform_y = y_transform)
tables_base <- build_tables(mg_sec_glm_base, sig_boot_base, ppi_term, hc3_fdr_cut = 0.2, boot_p_cut = 0.05)
mg_net_sec_res_list_base <- list(Net_secretion_HC3_base = mg_sec_glm_base, Sig_net_secretion_boot_base_output = sig_boot_base,
    tables_base = tables_base)
saveRDS(mg_net_sec_res_list_base, project_output("MGPIPE_Net_secretion_HC3_BASE_sigGLMpadj_02_BOOT_results.rds",
    name = "kiel_metabolite_secretion_basic_bootstrap.rds"))
write.csv(tables_base$res_clean_report, file = project_output(project_path("Results/mgpipe_net_secretion_HC3_BASE_clean_results.csv"),
    name = "kiel_metabolite_secretion_basic_bootstrap_summary.csv"), row.names = FALSE)
write.csv(tables_base$res_clean_report %>% filter(sig_both), file = project_output(project_path("Results/mgpipe_net_secretion_HC3_BASE_sig_both.csv"),
    name = "kiel_metabolite_secretion_basic_bootstrap_supported.csv"), row.names = FALSE)
message("FULL model: bootstrapped net secretion features = ", length(rxns_to_boot_full))
message("BASE model: bootstrapped net secretion features = ", length(rxns_to_boot_base))
message("FULL model: sig_both = ", sum(tables_full$res_clean_report$sig_both, na.rm = TRUE))
message("BASE model: sig_both = ", sum(tables_base$res_clean_report$sig_both, na.rm = TRUE))
