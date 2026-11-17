# SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0
# Cross-cohort taxonomic comparison statistics
# Human participant-level input matrices and covariates are not distributed.

source("scripts/helpers/setup.R")
invisible(run_analysis(
  section = "01_human_microbiome",
  directory = "human_cohorts/metacardis",
  script = "01_03_taxonomy_concordance.R",
  restricted = TRUE,
  code = {
maas <- read_tsv(project_input(project_path("human_cohorts/results/Maaslin2_FULL/all_results.tsv")))
maas_clean <- maas %>% dplyr::filter(metadata == "Proton.pump.inhibitor", value == "Yes") %>% dplyr::select(feature,
    maas_coef = coef, maas_se = stderr, maas_p = pval, maas_q = qval)
ald <- read.csv(project_input(project_path("Results/Aldex2/Focus_bacterial_abundance_aldex2_newPPI.csv")))
ald_clean <- ald %>% dplyr::select(feature = X, aldex_coef = PPI1.Est, aldex_se = PPI1.SE, aldex_p = PPI1.pval,
    aldex_q = PPI1.pval.padj)
extract_genus <- function(x) {
    gsub("[._].*", "", x)
}
ald_clean$genus <- extract_genus(ald_clean$feature)
maas_clean$genus <- extract_genus(maas_clean$feature)
ald_genus <- ald_clean %>% group_by(genus) %>% summarise(aldex_coef = mean(aldex_coef, na.rm = TRUE),
    aldex_q = min(aldex_q, na.rm = TRUE))
maas_genus <- maas_clean %>% group_by(genus) %>% summarise(maas_coef = mean(maas_coef, na.rm = TRUE),
    maas_q = min(maas_q, na.rm = TRUE))
merged <- inner_join(ald_genus, maas_genus, by = "genus")
merged$direction_agree <- sign(merged$aldex_coef) == sign(merged$maas_coef)
mean(merged$direction_agree)
cor.test(merged$aldex_coef, merged$maas_coef, method = "spearman")
sig_ald <- merged$aldex_q < 0.1
sig_maas <- merged$maas_q < 0.1
table(sig_ald, sig_maas)
clean_species <- function(x) {
    x <- gsub("\\.", "_", x)
    sapply(strsplit(x, "_"), function(parts) {
        if (length(parts) >= 2) {
            paste(parts[1], parts[2], sep = "_")
        }
        else {
            NA
        }
    })
}
ald_clean$species <- clean_species(ald_clean$feature)
maas_clean$species <- clean_species(maas_clean$feature)
merged_species <- inner_join(ald_clean, maas_clean, by = "species")
mean(sign(merged_species$aldex_coef) == sign(merged_species$maas_coef))
cor.test(merged_species$aldex_coef, merged_species$maas_coef, method = "spearman")
table(merged_species$aldex_q < 0.1, merged_species$maas_q < 0.1)
pal_sig <- c(`Not significant` = "#B8B8B8", `Kiel only` = "#5B8FA8", `Metacardis only` = "#C27D61", `Significant in both` = "#5A4E6D")
pal_direction <- c(`TRUE` = "#5A4E6D", `FALSE` = "#B8B8B8")
pal_cohort <- c(Kiel = "#5B8FA8", Metacardis = "#C27D61", `Meta-analysis` = "#2F2F2F")
merged$group <- "Not significant"
merged$group[merged$aldex_q < 0.1 & merged$maas_q < 0.1] <- "Significant in both"
merged$group[merged$aldex_q < 0.1 & merged$maas_q >= 0.1] <- "Kiel only"
merged$group[merged$aldex_q >= 0.1 & merged$maas_q < 0.1] <- "Metacardis only"
merged$group <- factor(merged$group, levels = c("Not significant", "Kiel only", "Metacardis only", "Significant in both"))
library(ggrepel)
sig_to_label <- subset(merged, group != "Not significant")
mb_plot <- ggplot(merged, aes(x = aldex_coef, y = maas_coef, color = group)) + geom_point(alpha = 0.8,
    size = 2) + geom_hline(yintercept = 0, linetype = "dashed", color = "black") + geom_vline(xintercept = 0,
    linetype = "dashed", color = "black") + geom_text_repel(data = sig_to_label, aes(label = genus),
    size = 3, show.legend = FALSE, max.overlaps = Inf) + scale_color_manual(values = pal_sig) + theme_classic() +
    labs(x = "ALDEx2 effect size (Kiel)", y = "MaAsLin2 effect size (MetaCardis)", color = "Association")
ggsave(filename = project_output(project_path("ancillary_figures/mb_plot.png"), name = "mb_plot.png"), plot = mb_plot,
    dpi = 600, width = 8, height = 8)
library(dplyr)
library(stringr)
clean_species <- function(x) {
    x <- gsub("\\.", "_", x)
    sapply(strsplit(x, "_"), function(parts) {
        if (length(parts) >= 2)
            paste(parts[1], parts[2], sep = "_")
        else NA_character_
    })
}
pretty_species <- function(x) {
    x <- gsub("_", " ", x)
    stringr::str_replace(x, "^([A-Za-z]+) ([a-z]+)$", "*\\1 \\2*")
}
ald_clean$species <- clean_species(ald_clean$feature)
maas_clean$species <- clean_species(maas_clean$feature)
ald_species <- ald_clean %>% filter(!is.na(species)) %>% group_by(species) %>% summarise(aldex_coef = mean(aldex_coef,
    na.rm = TRUE), aldex_se = mean(aldex_se, na.rm = TRUE), aldex_q = min(aldex_q, na.rm = TRUE), n_aldex = n(),
    genus = first(genus), .groups = "drop")
maas_species <- maas_clean %>% filter(!is.na(species)) %>% group_by(species) %>% summarise(maas_coef = mean(maas_coef,
    na.rm = TRUE), maas_se = mean(maas_se, na.rm = TRUE), maas_q = min(maas_q, na.rm = TRUE), n_maas = n(),
    genus = first(genus), .groups = "drop")
merged_species_unique <- inner_join(ald_species, maas_species, by = "species", suffix = c("_ald", "_maas"))
merged_species_unique$group <- "Not significant"
merged_species_unique$group[merged_species_unique$aldex_q < 0.1 & merged_species_unique$maas_q < 0.1] <- "Significant in both"
merged_species_unique$group[merged_species_unique$aldex_q < 0.1 & merged_species_unique$maas_q >= 0.1] <- "Kiel only"
merged_species_unique$group[merged_species_unique$aldex_q >= 0.1 & merged_species_unique$maas_q < 0.1] <- "Metacardis only"
merged_species_unique$group <- factor(merged_species_unique$group, levels = c("Not significant", "Kiel only",
    "Metacardis only", "Significant in both"))
merged_species_unique$direction_agree <- sign(merged_species_unique$aldex_coef) == sign(merged_species_unique$maas_coef)
merged_species_unique$species_label <- gsub("_", " ", merged_species_unique$species)
species_for_plot <- merged_species_unique %>% filter(aldex_q < 0.1 | maas_q < 0.1) %>% mutate(species_label = gsub("_",
    " ", species))
species_long <- species_for_plot %>% dplyr::select(species, species_label, aldex_coef, aldex_q, maas_coef,
    maas_q) %>% tidyr::pivot_longer(cols = c(aldex_coef, aldex_q, maas_coef, maas_q), names_to = c("cohort",
    ".value"), names_pattern = "(aldex|maas)_(coef|q)") %>% dplyr::mutate(cohort = dplyr::case_when(cohort ==
    "aldex" ~ "Kiel", cohort == "maas" ~ "Metacardis"), sig_level = dplyr::case_when(q < 0.1 ~ "FDR < 0.1",
    TRUE ~ "Not significant"))
species_order <- species_for_plot %>% mutate(mean_effect = (aldex_coef + maas_coef)/2) %>% arrange(mean_effect) %>%
    pull(species_label)
species_long$species_label <- factor(species_long$species_label, levels = species_order)
species_for_plot$species_label <- factor(species_for_plot$species_label, levels = species_order)
species_long$sig_level <- factor(species_long$sig_level, levels = c("Not significant", "FDR < 0.1"))
species_pairs <- ggplot() + geom_segment(data = species_for_plot, aes(x = aldex_coef, xend = maas_coef,
    y = species_label, yend = species_label), color = "#B8B8B8", linewidth = 0.6) + geom_point(data = species_long,
    aes(x = effect, y = species_label, color = cohort), size = 3) + geom_vline(xintercept = 0, linetype = "dashed",
    color = "black") + scale_color_manual(values = pal_cohort) + theme_classic(base_size = 12) + labs(x = "Adjusted effect size for PPI exposure",
    y = NULL, color = NULL)
sig_shapes <- c(`Not significant` = 1, `FDR < 0.1` = 16)
species_pairs <- ggplot() + geom_segment(data = species_for_plot, aes(x = aldex_coef, xend = maas_coef,
    y = species_label, yend = species_label), color = "#B8B8B8", linewidth = 0.6) + geom_point(data = species_long,
    aes(x = coef, y = species_label, color = cohort, shape = sig_level), size = 3) + geom_vline(xintercept = 0,
    linetype = "dashed", color = "black") + scale_color_manual(values = pal_cohort) + scale_shape_manual(values = sig_shapes) +
    theme_classic(base_size = 12) + labs(x = "Adjusted effect size for PPI exposure", y = NULL, color = NULL,
    shape = NULL)
ggsave(filename = project_output(project_path("ancillary_figures/species_pairs.png"), name = "species_pairs.png"),
    plot = species_pairs, dpi = 600, width = 8, height = 5.5)
ggsave(filename = project_output(project_path("ancillary_figures/species_pairs.png"), name = "species_pairs.png"),
    plot = species_pairs, dpi = 600, width = 8, height = 5.5)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(patchwork)
theme_set(theme_classic(base_family = "Arial", base_size = 30) + theme(text = element_text(family = "Arial"),
    panel.background = element_rect(fill = "white", color = NA), plot.background = element_rect(fill = "white",
        color = NA), legend.background = element_rect(fill = "white", color = NA), legend.box.background = element_rect(fill = "white",
        color = NA), axis.title = element_text(size = 22), axis.text = element_text(size = 18), legend.text = element_text(size = 18),
    legend.title = element_text(size = 20), plot.title = element_text(size = 24, face = "bold"), plot.margin = margin(12,
        12, 12, 12)))
mb_plot <- ggplot(merged, aes(x = aldex_coef, y = maas_coef, color = group)) + geom_point(alpha = 0.85,
    size = 3.5) + geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.7) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.7) + geom_text_repel(data = merged %>%
    filter(group != "Not significant") %>% slice_max(order_by = abs(aldex_coef) + abs(maas_coef), n = 10),
    aes(label = genus), size = 12, fontface = "bold", family = "Arial", show.legend = FALSE, max.overlaps = 10,
    box.padding = 0.6, point.padding = 0.4, segment.color = "grey50") + scale_color_manual(values = pal_sig) +
    labs(title = "Genus-level concordance", x = "ALDEx2 effect size (Kiel)", y = "MaAsLin2 effect size (MetaCardis)",
        color = "Association") + theme(legend.position = "bottom", axis.text.y = element_text(size = 30,
    family = "Arial"), axis.text.x = element_text(size = 30, family = "Arial"))
species_pairs <- ggplot() + geom_segment(data = species_for_plot, aes(x = aldex_coef, xend = maas_coef,
    y = species_label, yend = species_label), color = "#B8B8B8", linewidth = 0.8) + geom_point(data = species_long,
    aes(x = coef, y = species_label, color = cohort, shape = sig_level), size = 4.5) + geom_vline(xintercept = 0,
    linetype = "dashed", color = "black", linewidth = 0.7) + scale_color_manual(values = pal_cohort) +
    scale_shape_manual(values = c(`Not significant` = 1, `FDR < 0.1` = 16)) + labs(title = "Species-level effects",
    x = "Adjusted effect size for PPI exposure", y = NULL, color = NULL, shape = NULL) + theme(legend.position = "bottom",
    axis.text.y = element_text(size = 30, family = "Arial", face = "italic"))
ggsave(filename = project_output(project_path("ancillary_figures/genus_concordance.png"), name = "genus_concordance.png"),
    plot = mb_plot, width = 18, height = 8, units = "in", dpi = 600, bg = "white")
ggsave(filename = project_output(project_path("ancillary_figures/species_concordance.png"), name = "species_concordance.png"),
    plot = species_pairs, width = 16, height = 8, units = "in", bg = "white")
poster_combined <- (mb_plot | species_pairs) + plot_layout(widths = c(1.15, 1.35), guides = "collect") +
    theme(legend.position = "bottom", legend.box = "horizontal")
ggsave(filename = project_output(project_path("ancillary_figures/poster_taxa_combined_final.png"), name = "poster_taxa_combined.png"),
    plot = poster_combined, width = 16, height = 9, units = "in", dpi = 600, bg = "white")
ggsave(filename = project_output(project_path("ancillary_figures/poster_taxa_combined_final.pdf"), name = "poster_taxa_combined.pdf"),
    plot = poster_combined, width = 16, height = 9, units = "in", bg = "white")
}))
