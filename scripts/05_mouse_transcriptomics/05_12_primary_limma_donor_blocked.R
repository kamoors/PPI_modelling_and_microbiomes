# SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0
# Primary donor-blocked limma-voom contrasts excluding C02/G04.
# Run through run.R; inputs and exports use package-local routing.
# Requires the original/all-sample cache from 05_01; removes C02/G04 within the recovered limma function.
# Fitting helper bodies and settings are retained. Saved source output matches all 38,694 Supplementary Data 3 rows.
# The additional GSEA and original-versus-outlier comparison branches are not run here.

suppressPackageStartupMessages({
    library(DESeq2)
    library(SummarizedExperiment)
    library(edgeR)
    library(limma)
    library(dplyr)
    library(tibble)
})
load_rnaseq_script <- file.path(ppi_path("Validation/Mouse_RNASeq/Load_all_deseq2_files.R"))
ppi_source(load_rnaseq_script)
required_objects <- c("dds_list")
missing_objects <- required_objects[!vapply(required_objects, exists, logical(1))]
if (length(missing_objects) > 0) {
    stop("The following required objects are missing after source(load_rnaseq_script): ", paste(missing_objects, 
        collapse = ", "))
}
out_dir <- ppi_path("rnaseq/Results/RNAseq_limma_model_sensitivity_tests")
dir.create(ppi_output(out_dir, directory = TRUE), recursive = TRUE, showWarnings = FALSE)
exclude_samples <- c("25Jul1866-C02", "25Jul1866-G04")
DONOR_COLUMN <- "donor"
MOUSE_COLUMN <- "name"
GROUP_COLUMN <- "group"
CTRL_LEVEL <- "Ctrl"
PPI_LEVEL <- "PPI"
alpha <- 0.1
guess_column <- function(df, requested = NULL, candidates, label, required = TRUE) {
    if (!is.null(requested)) {
        if (!requested %in% colnames(df)) {
            stop(label, " column requested but not found: ", requested, "\nAvailable columns are:\n", 
                paste(colnames(df), collapse = ", "))
        }
        return(requested)
    }
    found <- candidates[candidates %in% colnames(df)]
    if (length(found) > 0) {
        message("Using ", label, " column: ", found[1])
        return(found[1])
    }
    if (required) {
        stop("Could not identify ", label, " column.\nAvailable columns are:\n", paste(colnames(df), 
            collapse = ", "), "\nSet the correct column manually at the top of the script.")
    }
    NULL
}
merge_external_meta <- function(meta_all, external_meta = NULL) {
    if (is.null(external_meta)) {
        return(meta_all)
    }
    ext <- as.data.frame(external_meta)
    ext$.meta_rownames <- rownames(ext)
    sample_candidates <- c("sampleID", "SampleID", "sample_id", "Sample_ID", "sample", "Sample", "ID", 
        "id", "sampleID_original", "Run", "run")
    sample_col <- sample_candidates[sample_candidates %in% colnames(ext)][1]
    if (is.na(sample_col)) {
        if (any(meta_all$sampleID_original %in% ext$.meta_rownames)) {
            ext$sampleID_original <- ext$.meta_rownames
            sample_col <- "sampleID_original"
        }
        else {
            warning("External meta object exists, but I could not identify a sample ID column. ", "Proceeding with colData(dds) only.")
            return(meta_all)
        }
    }
    ext <- ext %>% dplyr::rename(sampleID_original = all_of(sample_col)) %>% dplyr::distinct(sampleID_original, 
        .keep_all = TRUE)
    overlapping_cols <- intersect(colnames(meta_all), colnames(ext))
    overlapping_cols <- setdiff(overlapping_cols, "sampleID_original")
    if (length(overlapping_cols) > 0) {
        ext <- ext %>% dplyr::rename_with(~paste0("meta_", .x), all_of(overlapping_cols))
    }
    meta_all %>% dplyr::left_join(ext, by = "sampleID_original")
}
build_combined_rnaseq_object <- function(dds_list, exclude_samples = character(0), external_meta = NULL, 
    donor_col = NULL, mouse_col = NULL, group_col = NULL, ctrl_level = "Ctrl", ppi_level = "PPI") {
    if (is.null(names(dds_list)) || any(names(dds_list) == "")) {
        stop("dds_list must be a named list, with names corresponding to organs.")
    }
    common_genes <- Reduce(intersect, lapply(dds_list, rownames))
    if (length(common_genes) == 0) {
        stop("No common genes found across dds_list.")
    }
    message("Number of genes common to all organs: ", length(common_genes))
    counts_all <- list()
    meta_all <- list()
    for (org in names(dds_list)) {
        dds <- dds_list[[org]]
        keep_samples <- !colnames(dds) %in% exclude_samples
        dds <- dds[, keep_samples]
        cts <- DESeq2::counts(dds)[common_genes, , drop = FALSE]
        original_sample_id <- colnames(cts)
        combined_sample_id <- paste(org, original_sample_id, sep = "__")
        colnames(cts) <- combined_sample_id
        md <- as.data.frame(SummarizedExperiment::colData(dds))
        md$sampleID_original <- original_sample_id
        md$sampleID_combined <- combined_sample_id
        md$organ <- org
        rownames(md) <- combined_sample_id
        counts_all[[org]] <- cts
        meta_all[[org]] <- md
    }
    counts_all <- do.call(cbind, counts_all)
    meta_all <- dplyr::bind_rows(meta_all)
    meta_all <- meta_all[match(colnames(counts_all), meta_all$sampleID_combined), ]
    if (!all(meta_all$sampleID_combined == colnames(counts_all))) {
        stop("Metadata and count matrix are not aligned.")
    }
    meta_all <- merge_external_meta(meta_all = meta_all, external_meta = external_meta)
    group_candidates <- c("group", "Group", "treatment", "Treatment", "condition", "Condition", "PPI", 
        "ppi", "meta_group", "meta_Group", "meta_treatment", "meta_Treatment", "meta_condition", "meta_Condition")
    donor_candidates <- c("donor_id", "DonorID", "donorID", "donor", "Donor", "fecal_donor", "FecalDonor", 
        "Fecal_Donor", "HumanDonor", "human_donor", "meta_donor_id", "meta_DonorID", "meta_donorID", 
        "meta_donor", "meta_Donor", "meta_fecal_donor", "meta_FecalDonor")
    mouse_candidates <- c("mouse_id", "MouseID", "mouseID", "mouse", "Mouse", "animal_id", "AnimalID", 
        "animalID", "animal", "Animal", "recipient", "Recipient", "recipient_mouse", "RecipientMouse", 
        "meta_mouse_id", "meta_MouseID", "meta_mouseID", "meta_mouse", "meta_Mouse", "meta_animal_id", 
        "meta_AnimalID", "meta_animalID")
    group_col_use <- guess_column(df = meta_all, requested = group_col, candidates = group_candidates, 
        label = "group", required = TRUE)
    donor_col_use <- guess_column(df = meta_all, requested = donor_col, candidates = donor_candidates, 
        label = "donor", required = TRUE)
    mouse_col_use <- guess_column(df = meta_all, requested = mouse_col, candidates = mouse_candidates, 
        label = "mouse", required = FALSE)
    meta_all$group <- factor(meta_all[[group_col_use]])
    meta_all$group <- droplevels(meta_all$group)
    if (!all(c(ctrl_level, ppi_level) %in% levels(meta_all$group))) {
        stop("The group column does not contain both expected levels: ", ctrl_level, " and ", ppi_level, 
            "\nObserved group levels are: ", paste(levels(meta_all$group), collapse = ", "))
    }
    meta_all$group <- relevel(meta_all$group, ref = ctrl_level)
    meta_all$donor_id <- factor(meta_all[[donor_col_use]])
    if (!is.null(mouse_col_use)) {
        meta_all$mouse_id <- factor(meta_all[[mouse_col_use]])
    }
    else {
        meta_all$mouse_id <- NA
    }
    meta_all$organ <- factor(meta_all$organ, levels = names(dds_list))
    message("Number of samples after exclusion: ", ncol(counts_all))
    message("Number of donors: ", nlevels(meta_all$donor_id))
    message("Organs: ", paste(levels(meta_all$organ), collapse = ", "))
    message("Groups: ", paste(levels(meta_all$group), collapse = ", "))
    donor_group_check <- meta_all %>% dplyr::distinct(donor_id, group) %>% dplyr::count(donor_id, name = "n_groups_per_donor") %>% 
        dplyr::filter(n_groups_per_donor > 1)
    if (nrow(donor_group_check) > 0) {
        warning("Some donors appear in more than one treatment group. ", "Check whether donor_id and group are correctly defined.")
    }
    list(counts = counts_all, meta = meta_all, group_col = group_col_use, donor_col = donor_col_use, 
        mouse_col = mouse_col_use)
}
run_limma_voom_dupCor_donor_block <- function(dds_list, exclude_samples = character(0), external_meta = NULL, 
    donor_col = NULL, mouse_col = NULL, group_col = NULL, out_dir, label, ctrl_level = "Ctrl", ppi_level = "PPI", 
    min_count = 10, min_prop_samples = 0.2) {
    message("============================================================")
    message("Running limma-voom duplicateCorrelation: ", label)
    message("============================================================")
    obj <- build_combined_rnaseq_object(dds_list = dds_list, exclude_samples = exclude_samples, external_meta = external_meta, 
        donor_col = donor_col, mouse_col = mouse_col, group_col = group_col, ctrl_level = ctrl_level, 
        ppi_level = ppi_level)
    counts <- obj$counts
    meta_df <- obj$meta
    meta_df$organ_group <- factor(paste(meta_df$organ, meta_df$group, sep = "."))
    design <- model.matrix(~0 + organ_group, data = meta_df)
    colnames(design) <- make.names(levels(meta_df$organ_group))
    y <- edgeR::DGEList(counts = counts)
    keep_genes <- rowSums(y$counts >= min_count) >= ceiling(min_prop_samples * ncol(y))
    y <- y[keep_genes, , keep.lib.sizes = FALSE]
    message("Genes retained after filtering: ", nrow(y))
    y <- edgeR::calcNormFactors(y)
    pdf(ppi_output(file.path(out_dir, paste0("voom_mean_variance_", label, "_first_pass.pdf"))), width = 6, 
        height = 5)
    v0 <- limma::voom(y, design = design, plot = TRUE)
    dev.off()
    corfit <- limma::duplicateCorrelation(v0, design = design, block = meta_df$donor_id)
    message("Consensus donor correlation for ", label, ": ", round(corfit$consensus.correlation, 3))
    pdf(ppi_output(file.path(out_dir, paste0("voom_mean_variance_", label, "_dupCor.pdf"))), width = 6, 
        height = 5)
    v <- limma::voom(y, design = design, block = meta_df$donor_id, correlation = corfit$consensus.correlation, 
        plot = TRUE)
    dev.off()
    fit <- limma::lmFit(v, design = design, block = meta_df$donor_id, correlation = corfit$consensus.correlation)
    contrast_strings <- c()
    for (org in levels(meta_df$organ)) {
        ppi_col <- make.names(paste(org, ppi_level, sep = "."))
        ctrl_col <- make.names(paste(org, ctrl_level, sep = "."))
        if (ppi_col %in% colnames(design) && ctrl_col %in% colnames(design)) {
            contrast_strings[paste0(org, "_", ppi_level, "_vs_", ctrl_level)] <- paste0(ppi_col, " - ", 
                ctrl_col)
        }
        else {
            warning("Skipping contrast for ", org, ": missing ", ppi_col, " or ", ctrl_col)
        }
    }
    if (length(contrast_strings) == 0) {
        stop("No valid organ-specific PPI vs Ctrl contrasts could be constructed.")
    }
    contrast_matrix <- limma::makeContrasts(contrasts = contrast_strings, levels = design)
    fit2 <- limma::contrasts.fit(fit, contrast_matrix)
    fit2 <- limma::eBayes(fit2)
    res_by_organ <- list()
    for (coef_name in colnames(contrast_matrix)) {
        org <- sub(paste0("_", ppi_level, "_vs_", ctrl_level, "$"), "", coef_name)
        res_df <- limma::topTable(fit2, coef = coef_name, number = Inf, sort.by = "P") %>% tibble::rownames_to_column("gene") %>% 
            dplyr::transmute(gene = gene, log2FoldChange = logFC, AveExpr = AveExpr, stat = t, pvalue = P.Value, 
                padj = adj.P.Val, B = B, organ = org, sig_0.1 = !is.na(padj) & padj < 0.1, sig_0.05 = !is.na(padj) & 
                  padj < 0.05) %>% dplyr::arrange(padj)
        res_by_organ[[org]] <- list(res_df = res_df)
        write.csv(res_df, ppi_output(file.path(out_dir, paste0("limma_voom_dupCor_", org, "_", label, 
            ".csv"))), row.names = FALSE)
    }
    res_all <- dplyr::bind_rows(lapply(res_by_organ, `[[`, "res_df"))
    write.csv(res_all, ppi_output(file.path(out_dir, paste0("limma_voom_dupCor_all_organs_", label, ".csv"))), 
        row.names = FALSE)
    saveRDS(list(res_by_organ = res_by_organ, res_all = res_all, fit = fit, fit2 = fit2, voom = v, corfit = corfit, 
        design = design, contrast_matrix = contrast_matrix, meta = meta_df, exclude_samples = exclude_samples, 
        group_col = obj$group_col, donor_col = obj$donor_col, mouse_col = obj$mouse_col), ppi_output(file.path(out_dir, 
        paste0("limma_voom_dupCor_full_objects_", label, ".rds"))))
    list(res_by_organ = res_by_organ, res_all = res_all, fit = fit, fit2 = fit2, voom = v, corfit = corfit, 
        design = design, contrast_matrix = contrast_matrix, meta = meta_df, exclude_samples = exclude_samples, 
        group_col = obj$group_col, donor_col = obj$donor_col, mouse_col = obj$mouse_col)
}
external_meta <- if (exists("meta")) meta else NULL
limma_no_outliers <- run_limma_voom_dupCor_donor_block(dds_list = dds_list, exclude_samples = exclude_samples, 
    external_meta = external_meta, donor_col = DONOR_COLUMN, mouse_col = MOUSE_COLUMN, group_col = GROUP_COLUMN, 
    out_dir = out_dir, label = "without_C02_G04", ctrl_level = CTRL_LEVEL, ppi_level = PPI_LEVEL)
write.csv(limma_no_outliers$res_all, ppi_output(file.path(out_dir, "limma_all_organs_M1_dupCor_noCov.csv"), 
    name = "rna_limma_primary_all_organs.csv"), row.names = FALSE)
