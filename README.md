# PPI modelling and microbiomes

Code and processed data for *Proton pump inhibitor-associated microbiomes shape systemic host responses*, by Karlis A. Moors et al. Link to Preprint: https://www.biorxiv.org/content/10.64898/2026.11.10.765563v1

The numbered folders follow the order of the manuscript. Inputs are in `data/` and generated tables, figures and R objects go to `results/`.

## Getting started

Download this repository or clone it:

```sh
git clone https://github.com/kamoors/PPI_modelling_and_microbiomes.git
cd PPI_modelling_and_microbiomes
```

Use R 4.5. Open R or RStudio and set the working directory to this folder. Each analysis can then be run with `source()`. For example, this refits the serum LC-MS analysis from the prepared measurements:

```r
source("scripts/04_mouse_metabolomics/04_03_analyse_serum_lcms.R")
```

The same script can be run from a terminal in the repository folder:

```sh
Rscript scripts/04_mouse_metabolomics/04_03_analyse_serum_lcms.R
```

Start with a fresh R session for each analysis sequence. The scripts load their shared file-handling functions automatically; there is no installation step for the repository itself.

## R packages

For the metabolomics analyses:

```r
install.packages(c("BiocManager", "statmod", "dplyr", "tidyr", "purrr",
                   "stringr", "tibble", "readr", "readxl", "writexl"))
BiocManager::install("limma")
```

For experimental 16S analysis, add `phyloseq`, `vegan` and `permute`. For RNA-seq, add `DESeq2`, `SummarizedExperiment`, `edgeR` and `limma` through BiocManager. Enrichment scripts also require the annotation and gene-set packages named in those scripts, including `org.Mm.eg.db`, `clusterProfiler`, `fgsea`, `msigdbr` and `GSVA`. Other dependencies are listed in the individual scripts.

Figure scripts use `ggplot2`, `ggrepel`, `cowplot`, `patchwork`, `scales`, `forcats` and `ggsignif`, alongside the packages for the corresponding analysis. With ggplot2 4, use cowplot >= 1.2.0 and patchwork >= 1.3.2. The figures use Arial; Cairo support and this font are needed to reproduce the typography. Supplementary workbook assembly uses `openxlsx`.

## Analyses

Run the scripts below using `source("scripts/<folder>/<filename>")`. Filenames are shown in execution order. The public starting points use the supplied processed measurements or aggregate statistics.

| Folder | Scripts to run | Results |
| --- | --- | --- |
| `01_human_microbiome` | `01_11_meta_analyse_reaction_abundance.R`, then `01_12_meta_analyse_predicted_secretion.R` | Fixed-effect meta-analyses from the complete Kiel and MetaCardis cohort-level tables |
| `02_human_metabolism` | `02_02_blood_flux_meta_analysis.R`; `02_04_integrate_blood_flux_and_serum.R` | Blood-flux meta-analysis and comparison with saved serum results |
| `03_mouse_microbiome` | `03_02_plot_mouse_microbiome.R` | Experimental 16S ordination, Bray-Curtis comparisons and Figure 3; `03_03_export_bray_curtis_data.R` exports the distances without plotting |
| `04_mouse_metabolomics` | `04_03_analyse_serum_lcms.R`, `04_05_analyse_intestinal_lcms.R`, `04_06_analyse_caecum_nmr.R`, `04_07_analyse_ex_vivo_nmr.R`, `04_09_analyse_ex_vivo_hilic.R` | Primary model tables from prepared LC-MS/HILIC matrices and NMR intensities |
| `04_mouse_metabolomics` | `04_10_prepare_metabolomics_ordination.R`, then `04_11_test_metabolomics_ordination.R` | PCA configurations and donor-restricted Procrustes comparison; 04_11 can also use the supplied PCA configurations |
| `05_mouse_transcriptomics` | `05_01_prepare_primary_deseq2_cache.R`, then `05_04_primary_deseq2_and_gsea.R` | Primary DESeq2 results from tximport and mouse metadata, excluding C02 and G04 |
| `05_mouse_transcriptomics` | `05_01_prepare_primary_deseq2_cache.R`, then `05_12_primary_limma_donor_blocked.R` | Donor-blocked limma-voom results, excluding C02 and G04 |
| `06_multiomics_integration` | `06_02_prepare_multiomics_inputs.R`, `06_03_genus_metabolite_associations.R`, then `06_04_summarize_genus_metabolite_associations.R` | Genus-metabolite associations and leave-one-donor-out summaries; 06_03 can also start from the supplied prepared input object |

To repeat LC-MS feature import and curation, run 04_01 then 04_02 before 04_03; run 04_04 before 04_05; and run 04_08 before 04_09. The initial serum loading also runs historical mixed models and needs `lmerTest`.

Additional RNA analyses are in 05_03 (sample diagnostics), 05_05 (covariate sensitivity) and 05_06 (liver GO comparison using saved differential-expression tables). For gene-metabolite models, first run 05_09 to create the independent integration RNA cache, then 06_02 and 06_05. Script 05_11 computes pathway scores and pathway-metabolite associations from the integration cache and GSEA outputs. These analyses can take longer and require the enrichment dependencies. Script 06_06 exports the focused integration figure tables after the full integration fits.

Shared analysis functions are in `scripts/00_shared_helpers/`. Historical exploratory analyses and additional enrichment scripts are in `scripts/optional/`; they are not needed for the primary analyses above.

## Figures and supplementary tables

The figure scripts can run directly using the saved plotting tables included in `data/figure_inputs/`:

| Figure | Script |
| --- | --- |
| 1 | `scripts/01_human_microbiome/01_09_plot_human_microbiome.R` |
| 2 | `scripts/02_human_metabolism/02_05_plot_human_metabolism.R` |
| 3 | `scripts/03_mouse_microbiome/03_02_plot_mouse_microbiome.R` |
| 4 | `scripts/04_mouse_metabolomics/04_13_plot_mouse_metabolomics.R` |
| 5 | `scripts/05_mouse_transcriptomics/05_07_plot_mouse_transcriptomics.R` |
| 6 | `scripts/06_multiomics_integration/06_07_plot_multiomics_integration.R` |
| Supplementary transcriptomics figure | `scripts/05_mouse_transcriptomics/05_08_plot_supplementary_transcriptomics.R` |

For example:

```r
source("scripts/04_mouse_metabolomics/04_13_plot_mouse_metabolomics.R")
```

PDF and PNG files appear under the corresponding `results/<analysis>/figures/` folder. `results/files.csv` lists generated files and their producing scripts.

Figure generation uses saved tables where provided; it does not refit every underlying model. The reaction and secretion meta-analyses (01_11 and 01_12) update their Figure 1 inputs when run before plotting. Other figures may retain saved selection tables, PCA coordinates, annotations or pathway results. Newly generated outputs are used where their input references match; this is not an automatic rebuild of every manuscript result. See [data/README.md](data/README.md) for the saved inputs.

Run `scripts/07_supplementary_data/07_01_build_supplementary_data.R` to assemble the five Supplementary Data workbooks from the final result snapshots in `data/reference/`. This step formats those tables; it does not replace them with new fits.

## Data that are not included

Participant-level human cohort data and clinical covariates cannot be shared here. The Kiel and MetaCardis taxonomy, reaction-abundance, secretion, blood-flux and serum model scripts are included, but their individual-level models cannot be rerun from the public repository. Aggregate statistics support the meta-analyses and figures. The MetaCardis primary scripts retain exclusion of reported non-users with metabolomics-detected PPI.

Raw sequencing reads, instrument acquisition files and metabolic model reconstructions are not included. DADA2 preprocessing (03_01) requires reads and references. Oxygen-reaction extraction (01_07) requires AGORA2 models; set `PPI_AGORA2_MODELS` to their directory. The saved reaction statistics support threshold sensitivity analysis (01_08) without models. Model construction used Persephone from the Thiele lab; the MATLAB file in `02_human_metabolism/` records the study-specific host/microbiome assembly step and requires that external software and models.

Authorised users with the original private cohort directory structure can set `PPI_PRIVATE_DATA_ROOT` to that directory before running the restricted scripts. Keep these data outside the repository. Public analyses use only the files in this repository. Enrichment results can vary with annotation database versions, and stochastic plots can differ in layout.

## Licence and citation

The authors' code is licensed for noncommercial use under [PolyForm Noncommercial 1.0.0](LICENSE); see [NOTICE](NOTICE). Original processed data and aggregate tables are under CC BY-NC 4.0. Third-party resources retain their own terms; see [DATA_LICENSE.md](DATA_LICENSE.md).

Please cite the accompanying manuscript and this repository. Citation metadata are in [CITATION.cff](CITATION.cff).
