# SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0
# Figure 3: experimental microbiome ordination and Bray-Curtis distances.
# Run through run.R so input and output paths resolve within the package.
# Inputs: see data/input_manifest.csv and docs/RUN_ORDER.md.
# Writes generated files under results/; bundled inputs are read-only.

required_packages <- c("cowplot", "phyloseq", "dplyr", "ggplot2")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) {
    stop("Install required packages first: ", paste(missing_packages, collapse = ", "))
}
ppi_check_plot_versions(required_packages)
FIG_ROOT <- ppi_path("figures")
DATA_DIR <- file.path(FIG_ROOT, "source_data", "Experimental_16S")
OUT_DIR <- file.path(FIG_ROOT, "final")
dir.create(ppi_output(DATA_DIR, directory = TRUE), recursive = TRUE, showWarnings = FALSE)
dir.create(ppi_output(OUT_DIR, directory = TRUE), recursive = TRUE, showWarnings = FALSE)
REFRESH_SOURCE_DATA <- FALSE
PPI_ROOT <- ppi_path("")
stage_input <- function(filename, candidates) {
    target <- file.path(DATA_DIR, filename)
    if (REFRESH_SOURCE_DATA || !ppi_file_exists(target)) {
        hits <- candidates[ppi_file_exists(candidates)]
        if (!length(hits)) {
            stop("Could not stage ", filename, ".\nChecked:\n", paste(candidates, collapse = "\n"))
        }
        if (!file.copy(ppi_input(hits[[1]]), ppi_output(target), overwrite = TRUE)) 
            stop("Failed to copy ", hits[[1]])
        message("Staged: ", filename)
    }
    target
}
phyloseq_file <- stage_input("16S_PPI_preprocessing_withRepeats.RData", c(file.path(PPI_ROOT, "Validation", 
    "16S_results", "PPI_rawdata", "16S_PPI_preprocessing_withRepeats.RData")))
readr_manifest <- data.frame(file = basename(phyloseq_file), centralized_path = phyloseq_file)
utils::write.csv(readr_manifest, ppi_output(file.path(DATA_DIR, "Experimental_16S_source_manifest.csv"), 
    name = "experimental_16_s_source_manifest.csv"), row.names = FALSE)
load(ppi_input(phyloseq_file))
if (!exists("ps")) stop("The staged RData file does not contain a phyloseq object named 'ps'.")
BASE_SIZE <- 18
FIG_W_IN <- 18
FIG_H_IN <- 13.5
DPI <- 600
donor_cols <- c(`Ctrl 1` = "#FF6E00", `Ctrl 2` = "#8A3324", `Ctrl 3` = "#FFAD00", `PPI 1` = "#08306B", 
    `PPI 2` = "#2171B5", `PPI 3` = "#6BAED6")
group_cols <- c(nonrecipient = "white", recipient = "grey70")
base_theme <- ggplot2::theme_minimal(base_size = BASE_SIZE, base_family = "Arial") + ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
    panel.grid.major = ggplot2::element_line(colour = "#E6E6E6", linewidth = 0.25),
    axis.text = ggplot2::element_text(size = 17, colour = "#2F2F2F"), axis.title = ggplot2::element_text(size = 18, 
        colour = "#2F2F2F"), strip.text = ggplot2::element_text(size = 16), legend.title = ggplot2::element_text(size = 16), 
    legend.text = ggplot2::element_text(size = 15), plot.margin = ggplot2::margin(7, 7, 7, 7))
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
make_distance_panel <- function(sample_type) {
    # Select the function argument explicitly; the data column has the same name.
    ggplot2::ggplot(dplyr::filter(results_df, .data$sample_type == .env$sample_type), ggplot2::aes(group,
        distance, fill = group)) + ggplot2::geom_boxplot(width = 0.62, outlier.shape = NA, linewidth = 0.4) +
        ggplot2::geom_point(position = ggplot2::position_jitter(width = 0.08, height = 0, seed = 1),
            size = 1.2, alpha = 0.65, colour = "#333333") +
        ggplot2::facet_wrap(~donor_label, nrow = 2, scales = "free_y") + ggplot2::scale_fill_manual(values = group_cols,
        guide = "none") + ggplot2::scale_x_discrete(labels = c(nonrecipient = "Other", recipient = "Matched")) +
        ggplot2::labs(x = NULL, y = "Bray–Curtis dissimilarity") + base_theme + ggplot2::theme(legend.position = "none",
        strip.background = ggplot2::element_rect(fill = "#EEEEEE", colour = "#444444", linewidth = 0.35), 
        axis.text.x = ggplot2::element_text(size = 15))
}
p_culture <- make_distance_panel("Ex vivo")
p_mouse <- make_distance_panel("Mouse")
ps.culture.norm <- phyloseq::subset_samples(ps, experiment == "culture" | experiment == "OG_sample")
ps.culture.norm <- phyloseq::transform_sample_counts(ps.culture.norm, function(x) x/sum(x))
ord_culture <- phyloseq::ordinate(ps.culture.norm, method = "NMDS", distance = "bray")
ps.mouse.norm <- phyloseq::subset_samples(ps, experiment == "mouse" | experiment == "OG_sample")
ps.mouse.norm <- phyloseq::transform_sample_counts(ps.mouse.norm, function(x) x/sum(x))
ord_mouse <- phyloseq::ordinate(ps.mouse.norm, method = "NMDS", distance = "bray")
attach_donor_display <- function(ps_obj) {
    sd <- data.frame(phyloseq::sample_data(ps_obj))
    sd$subject_name <- as.character(sd$subject_name)
    ids <- rownames(sd)
    for (id in intersect(donor_ids_manual, ids)) {
        sd[id, "subject_name"] <- donor_groups[[id]]
    }
    sd$donor_label <- gsub("_", " ", sd$subject_name)
    sd$donor_label <- factor(sd$donor_label, levels = names(donor_cols))
    phyloseq::sample_data(ps_obj) <- phyloseq::sample_data(sd)
    ps_obj
}
ps.culture.norm <- attach_donor_display(ps.culture.norm)
ps.mouse.norm <- attach_donor_display(ps.mouse.norm)
# Construct one point layer explicitly: recipient samples are circles and the
# original donor samples are triangles, as in the manuscript caption.
culture_scores <- phyloseq::plot_ordination(ps.culture.norm, ord_culture, justDF = TRUE)
mouse_scores <- phyloseq::plot_ordination(ps.mouse.norm, ord_mouse, justDF = TRUE)
culture_nmds <- ggplot2::ggplot(culture_scores,
    ggplot2::aes(NMDS1, NMDS2, colour = donor_label, shape = experiment)) +
    ggplot2::geom_point(ggplot2::aes(size = time), alpha = 0.9) + ggplot2::scale_colour_manual(values = donor_cols, 
    drop = FALSE, name = "Donor") + ggplot2::scale_size_discrete(name = "Time", range = c(3, 5), na.translate = FALSE) + 
    ggplot2::scale_shape_manual(name = "Sample type", values = c(culture = 16, mouse = 16, OG_sample = 17),
        breaks = c("culture", "OG_sample"), labels = c("Culture/mouse", "Donor")) + ggplot2::labs(x = "NMDS1",
    y = "NMDS2", title = "Ex vivo") + base_theme + ggplot2::theme(legend.position = "bottom")
mouse_nmds <- ggplot2::ggplot(mouse_scores,
    ggplot2::aes(NMDS1, NMDS2, colour = donor_label, shape = experiment)) +
    ggplot2::geom_point(size = 4.2, alpha = 0.9) + ggplot2::scale_colour_manual(values = donor_cols, 
    drop = FALSE, name = "Donor") + ggplot2::scale_shape_manual(name = "Sample type",
    values = c(culture = 16, mouse = 16, OG_sample = 17), breaks = c("mouse", "OG_sample"),
    labels = c("Culture/mouse", "Donor")) + ggplot2::labs(x = "NMDS1", y = "NMDS2", title = "In vivo") +
    base_theme + ggplot2::theme(legend.position = "bottom")
column_title <- ggplot2::theme(plot.title = ggplot2::element_text(family = "Arial", face = "bold", size = 20, 
    hjust = 0.5, margin = ggplot2::margin(b = 8)))
culture_nmds <- culture_nmds + column_title
mouse_nmds <- mouse_nmds + column_title
# Extract one common legend instead of collecting context-dependent copies.
# This keeps the donor, time and sample-type keys inside the figure canvas.
culture_nmds <- culture_nmds + ggplot2::guides(
    size = ggplot2::guide_legend(order = 1),
    colour = ggplot2::guide_legend(order = 2, nrow = 1, override.aes = list(shape = 16, size = 3)),
    shape = ggplot2::guide_legend(order = 3, override.aes = list(size = 3)))
shared_legend <- cowplot::get_legend(culture_nmds + ggplot2::theme(
    legend.box = "horizontal", legend.position = "bottom"))
top_panels <- cowplot::plot_grid(
    culture_nmds + ggplot2::theme(legend.position = "none"),
    mouse_nmds + ggplot2::theme(legend.position = "none"),
    ncol = 2, labels = c("A", "B"), label_size = 24)
top_row <- cowplot::plot_grid(top_panels, shared_legend, ncol = 1, rel_heights = c(1, 0.12))
bottom_row <- cowplot::plot_grid(p_culture, p_mouse, ncol = 2, labels = c("C", "D"), label_size = 24)
final_figure <- cowplot::plot_grid(top_row, bottom_row, ncol = 1, rel_heights = c(1, 0.85))
print(final_figure)
pdf_file <- file.path(OUT_DIR, "Experimental_16S.pdf")
png_file <- file.path(OUT_DIR, "Experimental_16S.png")
ggplot2::ggsave(ppi_output(pdf_file), final_figure, width = FIG_W_IN, height = FIG_H_IN, units = "in", 
    device = grDevices::cairo_pdf, limitsize = FALSE)
if (requireNamespace("ragg", quietly = TRUE)) {
    ragg::agg_png(ppi_output(png_file), width = FIG_W_IN, height = FIG_H_IN, units = "in", res = DPI, 
        background = "white")
    print(final_figure)
    grDevices::dev.off()
} else {
    ggplot2::ggsave(ppi_output(png_file), final_figure, width = FIG_W_IN, height = FIG_H_IN, units = "in", 
        dpi = DPI, bg = "white", limitsize = FALSE)
}
readr_like <- results_df
utils::write.csv(readr_like, ppi_output(file.path(DATA_DIR, "Figure3_BrayCurtis_plot_data.csv"), name = "figure3_bray_curtis_plot_data.csv"), 
    row.names = FALSE)
message("Figure 3 written to: ", OUT_DIR)
