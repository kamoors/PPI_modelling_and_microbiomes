# Figure reproduction

The seven plotting scripts were traced to `Manuscript/Figures/scripts/` and run from the public package inputs. All six main figures and the supplementary transcriptomics figure generated PDF and PNG files. The manuscript source scripts have not changed since the earlier tracing pass.

## Input coverage

All 50 source-data assets are represented:

- 40 files are byte-identical copies.
- Six gene-result files retain every field of the exact 3-HPP rows used by the cross-layer figure. Their retained rows match the full originals exactly.
- The experimental phyloseq workspace contains only `ps`; unused birth-date and sex metadata were removed. Counts, taxonomy, sample order and retained metadata match.
- The saved oxygen plot was replaced by its 6,357-row aggregate table. Every column matches; captured plotting environments are excluded.
- Two large RNA caches were replaced by compact plotting tables. The three PCA coordinate tables and variance attributes match the original figure function. All 14,070 Liver, 15,987 Colon and 16,041 Brain gene-symbol assignments match its annotation function. Saved supplementary GO/expression tables are included unchanged.

Source and coverage records are in `provenance/figure_script_sources.csv`, `figure_input_coverage.csv` and `transformed_input_checks.csv`.

## Experimental 16S plotting correction

The current source script used `filter(results_df, .data$sample_type == sample_type)` inside a function whose argument was also named `sample_type`. Both sides resolved to the data column, so both distance panels selected all contexts. The package now uses `.env$sample_type` to select the function argument explicitly. This restores the separate comparisons in the saved manuscript figure.

The culture panel contains 144 donor/sample comparisons (four recipients and 20 nonrecipients per donor); the mouse panel contains 96 (two or three recipients per donor, with the remaining mice as nonrecipients). The donor-comparison Wilcoxon values were checked separately. The current manuscript image has no P-value brackets; the drawing now follows that version. No ASV counts, abundance transformations, Bray-Curtis values or fitted model results changed. A regression check is provided:

```sh
Rscript tests/check_microbiome_panels.R
```

The donor facets use two rows, as in the saved figure. One shared legend replaces clipped duplicate legends. The original manuscript scripts and figures were not overwritten.

## Rendering environment and visual limits

Successful rendering used R 4.5.3, ggplot2 4.0.1, cowplot 1.2.0, patchwork 1.3.2, forcats 1.0.1 and Rcpp 1.1.2. The cohort plot used clusterProfiler 4.18.4 and DOSE 4.4.0, with enrichplot 1.24.4 and ggtree 3.12.0 satisfying its declared dependencies. Dependency versions are recorded in `environment/figure_render_environment.csv`. The isolated test libraries are outside the repository and are not redistributed.

Rendered figures were visually compared with the saved final figures. The WBM PNG in the source folder is a cropped panel label, so its PDF was used for reference. Figures are not promised to be byte-identical: the original scripts leave jitter, label repulsion and NMDS initialization unseeded, and typography/layout can vary with graphics-package versions. Numerical inputs and the corrected 16S comparisons were checked separately.

Regenerating figures from saved results does not independently revalidate the unavailable human cohort fits or the full RNA/gene-metabolite enrichment workflows. Those limits remain in `VALIDATION.md`.

The 12 November follow-up directly inspected the images embedded in the closed manuscript and corrected Figure 3 presentation and Figure 4 legend spacing. See FIGURE_LAYOUT_CORRECTIONS.md.
