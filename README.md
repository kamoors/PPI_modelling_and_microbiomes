# PPI modelling and microbiomes

Analysis code and processed data accompanying **Proton pump inhibitor-associated microbiomes shape systemic host responses** by Karlis A. Moors and colleagues.

The scripts follow manuscript order: human microbiome analyses (01), human metabolic analyses (02), experimental 16S (03), metabolomics (04), mouse transcriptomics (05), cross-layer integration (06), and supplementary table assembly (07). Shared functions are in `scripts/00_shared_helpers/`.

## Included data

The repository contains processed experimental 16S counts and taxonomy, metabolomics feature intensities and annotations, gene-level tximport matrices and mouse metadata, and aggregate cohort result tables needed for the figures. Each input is listed in [`data/input_manifest.csv`](data/input_manifest.csv), with its source, purpose and SHA-256 checksum.

Participant-level human cohort data cannot be provided in this repository. This includes secretion and reaction-abundance matrices and their demographic or clinical covariates. Their analysis scripts are retained, with the restriction explained in [`docs/RESTRICTED_DATA.md`](docs/RESTRICTED_DATA.md). Raw sequencing reads, raw instrument files and metabolic model reconstructions are not bundled.

## Running the code

Use R 4.5 and install the packages for the workflow you want to run; see [`environment/README.md`](environment/README.md). Open a terminal in the repository folder:

```sh
Rscript run.R --list
Rscript run.R --check
Rscript tests/check_inputs.R
Rscript run.R --workflow cohort-meta
Rscript run.R --workflow metabolomics-refit
Rscript run.R --workflow rna-primary
Rscript run.R --workflow integration
Rscript run.R --workflow figures
Rscript tests/compare_generated_tables.R results
```

Individual scripts can be selected by their numbered filename:

```sh
Rscript run.R --script 04_11_test_metabolomics_ordination.R
```

The repository can be moved to another folder. Inputs default to `data/` and generated files go to `results/`. No original project directory is required. `results/output_manifest.csv` records the producer and path of each export. Use `--results-root /path/to/results` for another destination. See [`docs/RUN_ORDER.md`](docs/RUN_ORDER.md) for prerequisites and supported starting points.

## Reproduction scope

Refitting from processed experimental measurements, redrawing figures from saved tables, and assembling supplementary workbooks are separate operations. `data/reference/` contains published result snapshots for comparison and workbook assembly; these are not substituted for measurements in the model-fitting scripts. The supplementary assembly step does not automatically incorporate new model fits.

The primary RNA workflow preserves the original filtering order, including exclusions of C02 and G04. The independent integration cache uses its original, different filtering order. Optional enrichment and historical additional analyses are retained in `scripts/optional/`; these may require annotation packages, databases and locally generated caches.

The MetaCardis reaction and secretion scripts reconstruct the primary models that match the downstream tables. Both retain exclusion of recorded non-users with metabolomics-detected PPI. Their provenance and validation are described in [`docs/SOURCE_UPDATES_20261113.md`](docs/SOURCE_UPDATES_20261113.md). The `cohort-meta` workflow reruns both meta-analyses using bundled aggregate cohort statistics and registers refreshed Figure 1 inputs.

Model construction used the external Persephone pipeline from the Thiele lab. Its source is referenced rather than redistributed; the study-specific assembly wrapper and documented settings are included. See [`docs/MODEL_CONSTRUCTION.md`](docs/MODEL_CONSTRUCTION.md). The reviewed manuscript and Supplementary Data 1–5 do not report sex-stratified association or PPI-by-sex interaction results. Those exploratory branches are outside this release's required workflow. See [`docs/VALIDATION.md`](docs/VALIDATION.md). Figure-specific checks and the corrected 16S panel filter are described in [`docs/FIGURE_VALIDATION.md`](docs/FIGURE_VALIDATION.md).

## Licence and citation

The authors' code is available under the **PolyForm Noncommercial License 1.0.0**; see [`LICENSE`](LICENSE) and [`NOTICE`](NOTICE). Commercial use is not licensed by these terms. This is a source-available repository, not an OSI open-source release. Original experimental data are provided under **CC BY-NC 4.0**; third-party annotations retain their own terms, as described in [`DATA_LICENSE.md`](DATA_LICENSE.md).

Please cite the accompanying manuscript and this repository: https://github.com/kamoors/PPI_modelling_and_microbiomes. A manuscript DOI can be added to `CITATION.cff` when available.
