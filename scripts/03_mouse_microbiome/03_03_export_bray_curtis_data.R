# SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0
# Export the Figure 3 Bray-Curtis data without plotting dependencies.
# Calculations and row matching are copied unchanged from 03_02.
library(phyloseq)
library(dplyr)
ps <- readRDS(ppi_input(ppi_path("microbiome16s/processed_phyloseq.rds")))
donor_cols <- c(`Ctrl 1`="#FF6E00",`Ctrl 2`="#8A3324",`Ctrl 3`="#FFAD00",
                `PPI 1`="#08306B",`PPI 2`="#2171B5",`PPI 3`="#6BAED6")
phyloseq::sample_data(ps)$time <- "D0"
for (d in 1:4) {
    hit <- grepl(paste0("D", d), phyloseq::sample_data(ps)$sample_name)
    phyloseq::sample_data(ps)$time[hit] <- paste0("D", d)
}
phyloseq::sample_data(ps)$time <- factor(phyloseq::sample_data(ps)$time, levels = c("D0", "D1", "D2", 
    "D3", "D4"))
sample_meta <- as.data.frame(phyloseq::sample_data(ps))
sample_meta$subject_name <- as.character(sample_meta$subject_name)
sample_meta$experiment <- as.character(sample_meta$experiment)
rownames(sample_meta) <- sample_meta$barcode
donor_ids_manual <- c("250200001321", "250200001322", "250200001323", "250200001324", "250200001325", 
    "250200001326")
donor_groups <- c("PPI_1", "PPI_2", "PPI_3", "Ctrl_1", "Ctrl_2", "Ctrl_3")
names(donor_groups) <- donor_ids_manual
rownames(sample_meta) <- as.character(sample_meta$barcode)
present_donors <- intersect(donor_ids_manual, rownames(sample_meta))
sample_meta[present_donors, "subject_name"] <- donor_groups[present_donors]
sample_meta <- as.data.frame(phyloseq::sample_data(ps))
sample_meta$subject_name <- as.character(sample_meta$subject_name)
sample_meta$experiment <- as.character(sample_meta$experiment)
rownames(sample_meta) <- phyloseq::sample_names(ps)
donor_ids <- rownames(sample_meta)[sample_meta$experiment == "OG_sample"]
ps_norm <- phyloseq::transform_sample_counts(ps, function(x) x/sum(x))
bray_mat <- as.matrix(phyloseq::distance(ps_norm, method = "bray", type = "samples"))
rownames(bray_mat) <- as.character(rownames(bray_mat))
colnames(bray_mat) <- as.character(colnames(bray_mat))
distance_rows <- list()
k <- 1L
for (donor_id in donor_ids) {
    donor_group <- unname(donor_groups[donor_id])
    if (is.na(donor_group)) 
        next
    mouse_rec <- as.character(sample_meta$barcode[sample_meta$subject_name == donor_group & sample_meta$experiment == 
        "mouse"])
    culture_rec <- as.character(sample_meta$barcode[sample_meta$subject_name == donor_group & sample_meta$experiment == 
        "culture" & sample_meta$time != "D0"])
    mouse_non <- as.character(sample_meta$barcode[sample_meta$subject_name != donor_group & sample_meta$experiment == 
        "mouse"])
    culture_non <- as.character(sample_meta$barcode[sample_meta$subject_name != donor_group & sample_meta$experiment == 
        "culture" & sample_meta$time != "D0"])
    add_rows <- function(ids, sample_type, group) {
        ids <- intersect(ids, colnames(bray_mat))
        if (!length(ids)) 
            return(NULL)
        data.frame(donor_group = donor_group, donor_id = donor_id, sample_type = sample_type, group = group, 
            distance = as.numeric(bray_mat[donor_id, ids]), stringsAsFactors = FALSE)
    }
    for (z in list(add_rows(mouse_rec, "Mouse", "recipient"), add_rows(culture_rec, "Ex vivo", "recipient"), 
        add_rows(mouse_non, "Mouse", "nonrecipient"), add_rows(culture_non, "Ex vivo", "nonrecipient"))) {
        if (!is.null(z)) {
            distance_rows[[k]] <- z
            k <- k + 1L
        }
    }
}
results_df <- dplyr::mutate(dplyr::bind_rows(distance_rows), donor_label = gsub("_", " ", donor_group), 
    donor_label = factor(donor_label, levels = names(donor_cols)), group = factor(group, levels = c("nonrecipient", 
        "recipient")))
write.csv(results_df,ppi_output(ppi_path("figures/source_data/Experimental_16S/Figure3_BrayCurtis_plot_data.csv"),
    name="figure3_bray_curtis_plot_data.csv"),row.names=FALSE)
message("Figure 3 Bray-Curtis table exported from processed counts.")
