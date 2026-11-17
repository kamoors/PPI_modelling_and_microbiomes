# SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0
# Configuration

project_dir <- project_path("")
metabolomics_root <- file.path(project_dir, "metabolomics")
rnaseq_root <- file.path(project_dir, "rnaseq")
association_root <- file.path(project_dir, "integration")
strict_input_paths <- TRUE
input_paths <- list(phyloseq = project_path("microbiome16s/16S_PPI_preprocessing_withRepeats.rds"),
    pathway_scores = file.path(rnaseq_root, "Results", "hallmark_reactome_metabolite_correlations_without_C02_G04",
        "Hallmark_Reactome_ssGSEA_scores_all_organs.csv"), rna_loader_script = project_code("scripts/05_mouse_transcriptomics/05_10_load_integration_deseq2_cache.R"),
    rna_metadata_table = file.path(rnaseq_root, "analysis_sample_metadata.csv"), serum_lcms = file.path(metabolomics_root, "Serum_metabolomics",
        "20260121_data", "Serum_cleaned_after_collaborator_review", "Serum_LCMS_limma_input_objects.rds"),
    cecum_feces_lcms = file.path(metabolomics_root, "20250714_MassSpecData", "Cecum_Feces_cleaned_after_collaborator_review",
        "Cecum_Feces_LCMS_limma_input_objects.rds"), cecum_nmr = file.path(metabolomics_root, "NMR_data",
        "NMR_cecum_annotated_limma", "NMR_cecum_annotated_limma_results.xlsx"), exvivo_nmr_raw = file.path(metabolomics_root,
        "NMR_data", "NMR_NEW_raw_data.csv"), exvivo_nmr_mapping = file.path(metabolomics_root, "NMR_data",
        "NMR_exvivo_annotated_limma_dupCor", "NMR_exvivo_cleaned_feature_mapping_v2.csv"), exvivo_hilic = file.path(metabolomics_root,
        "20250714_MassSpecData", "Ex_vivo_cleaned_one_per_metabolite", "HILIC_exvivo_cleaned_analysis_objects.rds"))
gsea_paths <- list(Colon_HALLMARK = file.path(rnaseq_root, "Results", "RNAseq_without_C02_G04", "GSEA_Colon_HALLMARK_without_C02_G04.csv"),
    Colon_REACTOME = file.path(rnaseq_root, "Results", "RNAseq_without_C02_G04", "GSEA_Colon_REACTOME_without_C02_G04.csv"),
    Liver_HALLMARK = file.path(rnaseq_root, "Results", "RNAseq_without_C02_G04", "GSEA_Liver_HALLMARK_without_C02_G04.csv"),
    Liver_REACTOME = file.path(rnaseq_root, "Results", "RNAseq_without_C02_G04", "GSEA_Liver_REACTOME_without_C02_G04.csv"),
    Brain_HALLMARK = file.path(rnaseq_root, "Results", "RNAseq_without_C02_G04", "GSEA_Brain_HALLMARK_without_C02_G04.csv"),
    Brain_REACTOME = file.path(rnaseq_root, "Results", "RNAseq_without_C02_G04", "GSEA_Brain_REACTOME_without_C02_G04.csv"))
settings <- list(gsea_fdr_cutoff = 0.1, association_fdr_cutoff = 0.1, exvivo_days = 1:4, taxonomic_rank = "Genus",
    genus_prevalence_fraction = 0.25, genus_mean_relative_abundance = 0.001, genus_minimum_donors = 2L,
    genus_pseudocount = 0.5, metabolite_minimum_positive_samples = 2L, metabolite_minimum_positive_fraction = 0.25,
    metabolite_maximum_missing_fraction = 0.25, minimum_samples_mouse = 8L, minimum_samples_exvivo = 12L,
    minimum_donors = 4L, rank_standardize = TRUE, robust_ebayes = TRUE, trend_ebayes = FALSE, heatmap_max_predictors = 30L,
    heatmap_max_outcomes = 30L, scatterplot_top_n_per_family = 3L, expected_mouse_16s_samples = 16L,
    expected_culture_16s_samples_all_days = 30L, expected_culture_16s_samples_days1_4 = 24L, expected_original_donor_16s_samples = 6L,
    expected_unmatched_16s_samples = 0L, fail_on_audit_mismatch = TRUE, run_leave_one_donor_out = TRUE)
run_layers <- list(metabolite_pathway = TRUE, metabolite_genus = TRUE, genus_pathway = TRUE, supported_triplets = TRUE)
mouse_metabolomics_datasets <- c("serum_lcms", "cecum_lcms", "feces_lcms", "cecum_nmr")
exvivo_metabolomics_datasets <- c("exvivo_nmr", "exvivo_hilic")
