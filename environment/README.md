# Software environment

R 4.5.3 was used for the local packaging checks. `packages.csv` records versions recovered from the analysis environment; it is not a tested lockfile. `validation_versions.csv` contains earlier validation records. `run.R --check` reports installed versions and missing namespaces without installing anything.

For the primary metabolomics analyses, install:

```r
install.packages(c("BiocManager", "statmod", "dplyr", "tidyr", "purrr", "stringr", "tibble", "readr", "readxl", "writexl"))
BiocManager::install("limma")
```

For primary RNA refitting add `DESeq2` and `SummarizedExperiment` through BiocManager. The tximport object is already provided, so refitting does not require re-running Salmon or loading the tximport package. For 16S add phyloseq and vegan; for Procrustes add vegan and permute. Plotting additionally uses ggplot2, ggrepel, cowplot, patchwork, scales, forcats and ggsignif. Supplementary assembly uses openxlsx.

Full enrichment, sensitivity, upstream and optional analyses require their additional packages listed at the top of each script, including org.Mm.eg.db, clusterProfiler, fgsea, msigdbr, GSVA and lmerTest where applicable. Dependencies and third-party databases retain their own licences. Annotation/database updates can change gene labels, memberships and enrichment results.

Figure scripts retain the manuscript's Arial font. Headless plotting uses a Cairo device for font-aware layout; install Arial or an equivalent system font for matching typography. This does not affect numerical model results.

The recorded plotting combination is ggplot2 4.0.3, cowplot 1.2.0 and patchwork 1.3.2. When using ggplot2 4, use at least cowplot 1.2.0 and patchwork 1.3.2. Older assemblers were present locally during validation and caused a patchwork operator failure or missing legends. Figure scripts check these minimum versions before drawing. See the [patchwork changelog](https://patchwork.data-imaginist.com/news/index.html) for compatibility updates.

The follow-up figure check successfully rendered all seven figures in an isolated compatible library. Actual versions, including the enrichment visualization dependencies used for the cohort plot, are in `figure_render_environment.csv`; see `docs/FIGURE_VALIDATION.md`.

The recovered primary limma branch additionally needs edgeR, limma, DESeq2 and SummarizedExperiment. It loads the original/all-sample cache created by 05_01. The private Kiel blood script needs R.matlab (3.7.0 was used for the raw-file validation); its MAT import temporarily uses the C character locale on Windows and restores the prior locale afterwards. MATLAB and the separately installed Persephone/COBRA Toolbox are needed only for the optional external model-construction wrapper. No third-party pipeline source is redistributed.
