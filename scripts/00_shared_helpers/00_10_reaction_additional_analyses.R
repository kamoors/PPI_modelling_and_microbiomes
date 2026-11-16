# SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0
# Additional FULL-bootstrap and BASE-model exports
# Run through run.R so input and output paths resolve within the package.
# Inputs: see data/input_manifest.csv and docs/RUN_ORDER.md.
# Writes generated files under results/; bundled inputs are read-only.

hc3_p_full <- sapply(react_ab_glm_full, function(ct) ct[ppi_term, "Pr(>|t|)"])
hc3_fdr_full <- p.adjust(hc3_p_full, method = "fdr")
rxns_to_boot_full <- names(hc3_fdr_full)[hc3_fdr_full < 0.2]
data_in_full_sig <- data_in[, rxns_to_boot_full, drop = FALSE]
sig_react_ab_boot_full <- bootstrap_subset(data_in_full_sig, test_df, fml_boot_full, ppi_term, R = 1000)
tables_full <- build_tables(react_ab_glm_full, sig_react_ab_boot_full, ppi_term)
react_ab_res_list_full <- list(Reaction_abundance_HC3_full = react_ab_glm_full, Sig_react_ab_boot_full_output = sig_react_ab_boot_full, 
    tables_full = tables_full)
saveRDS(react_ab_res_list_full, ppi_output("MGPIPE_CLR_reaction_abundance_HC3_FULL_sigGLMpadj_02_BOOT_results.rds", 
    name = "kiel_reaction_abundance_full_bootstrap.rds"))
write.csv(tables_full$res_clean_report, file = ppi_output(ppi_path("Results/mgpipe_reaction_abundance_CLR_HC3_FULL_clean_results.csv"), 
    name = "kiel_reaction_abundance_full_bootstrap_summary.csv"), row.names = FALSE)
write.csv(tables_full$res_clean_report %>% filter(sig_both), file = ppi_output(ppi_path("Results/mgpipe_reaction_abundance_CLR_HC3_FULL_sig_both.csv"), 
    name = "kiel_reaction_abundance_full_bootstrap_supported.csv"), row.names = FALSE)
react_ab_glm_base <- fit_all_reactions_hc3(data_in, test_df, fml_base)
hc3_p_base <- sapply(react_ab_glm_base, function(ct) ct[ppi_term, "Pr(>|t|)"])
hc3_fdr_base <- p.adjust(hc3_p_base, method = "fdr")
rxns_to_boot_base <- names(hc3_fdr_base)[hc3_fdr_base < 0.2]
data_in_base_sig <- data_in[, rxns_to_boot_base, drop = FALSE]
sig_react_ab_boot_base <- bootstrap_subset(data_in_base_sig, test_df, fml_boot_base, ppi_term, R = 1000)
tables_base <- build_tables(react_ab_glm_base, sig_react_ab_boot_base, ppi_term)
react_ab_res_list_base <- list(Reaction_abundance_HC3_base = react_ab_glm_base, Sig_react_ab_boot_base_output = sig_react_ab_boot_base, 
    tables_base = tables_base)
saveRDS(react_ab_res_list_base, ppi_output("MGPIPE_CLR_reaction_abundance_HC3_BASE_sigGLMpadj_02_BOOT_results.rds", 
    name = "kiel_reaction_abundance_basic_bootstrap.rds"))
write.csv(tables_base$res_clean_report, file = ppi_output(ppi_path("Results/mgpipe_reaction_abundance_CLR_HC3_BASE_clean_results.csv"), 
    name = "kiel_reaction_abundance_basic_bootstrap_summary.csv"), row.names = FALSE)
write.csv(tables_base$res_clean_report %>% filter(sig_both), file = ppi_output(ppi_path("Results/mgpipe_reaction_abundance_CLR_HC3_BASE_sig_both.csv"), 
    name = "kiel_reaction_abundance_basic_bootstrap_supported.csv"), row.names = FALSE)
message("FULL model: reactions bootstrapped = ", length(rxns_to_boot_full))
message("BASE model: reactions bootstrapped = ", length(rxns_to_boot_base))
message("FULL model: sig_both = ", sum(tables_full$res_clean_report$sig_both, na.rm = TRUE))
message("BASE model: sig_both = ", sum(tables_base$res_clean_report$sig_both, na.rm = TRUE))
