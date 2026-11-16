# Code provenance and packaging changes

`original_scripts.json` preserves the source tracing and checks from the earlier curated package. Its portable hashes describe that earlier package, not this release. `script_manifest.csv` records this release's script hashes and links to those earlier sources. `packaging_changes.json` records the initial copying pass; its initial destination hashes are superseded by `script_manifest.csv`.

Changes in this release:

- Noncommercial licence notices, package-local input mapping, relative working paths and descriptive output names. No original project root is searched automatically.
- Removal of top-level interactive inspection calls and directory-guessing code; explicit loading of prerequisites previously supplied by an RStudio session.
- Compact data tables replace serialized plotting environments, fitted PCA objects and fitted RNA caches. RNA caches are regenerated locally from tximport. Frozen gene annotations replace live annotation lookup for the primary cache's labels; the fitting expressions are unchanged.
- The primary RNA contrast script retains the original fitting/filtering code. Its additional plotting/GO/GSEA branch is retained under `scripts/optional/`. The independent integration cache remains separate.
- The primary Procrustes script retains the original test function. Broader historical PCA statistics are retained as optional code. They are not added to the manuscript supplementary workbooks.
- Figure 5 uses exact exported PCA coordinate tables. Supplementary transcriptomics plotting uses the published GO and normalized expression tables. The original refitting branch is retained separately.
- 03_03 extracts the unchanged Figure 3 Bray-Curtis calculations into a plotting-independent table exporter.
- Supplementary assembly uses the final five-workbook table snapshots. It is a reporting operation and does not refresh statistics after refitting.
- Unused birth-date and sex fields were removed from public experimental 16S metadata and matching copies. ASV counts and taxonomy were verified unchanged. Cohort participant matrices and covariates were never selected for inclusion.
- Compatible plot-assembler versions are required with ggplot2 4 to prevent invalid rendering. Fonts, labels, numerical transformations, model formulas, contrasts and thresholds were not intentionally changed.

See docs/VALIDATION.md for what was tested and the remaining limitations. Original manuscript and analysis files were not modified by this repository build.

Follow-up figure validation corrected the experimental 16S context filter to select the function argument explicitly, restored two-row donor facets and one unclipped shared legend, and added a regression check. These drawing changes restore the separate saved mouse/culture comparisons; the distance table and model results are unchanged. See docs/FIGURE_VALIDATION.md.

Version 0.4.2 follows the actual manuscript-embedded Figure 3 presentation and corrects metabolomics legend clipping/colour-bar spacing. These are drawing changes; input and statistical tables are unchanged. See docs/FIGURE_LAYOUT_CORRECTIONS.md.

Version 0.5.0 records newly recovered sources and validated MetaCardis/meta-analysis reconstructions in source_updates_20261113.json. The earlier original_scripts.json remains a historical record of the original curated sources. New MetaCardis primary outputs and both meta-analyses match the saved downstream numerical tables. Nine recovered helper bodies were checked after removing filesystem wrappers only. The MetaCardis blood script omits unreported sex-stratified and interaction branches; retained expressions were verified identical. The current manuscript and five supplementary workbooks were checked for this scope decision. Further full enrichment and integration reruns were outside the requested update. The external Persephone pipeline is referenced, with only its study-specific assembly wrapper bundled. See docs/SOURCE_UPDATES_20261113.md.
