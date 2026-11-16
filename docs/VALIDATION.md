# Validation and reproduction limits

Initial local checks were performed on 11 November 2026 with R 4.5.3 on Windows. Follow-up source recovery, cohort-model validation and package updates were performed on 13 November; see SOURCE_UPDATES_20261113.md. Numerical comparisons use an absolute tolerance of `1e-12 + 1e-8 * abs(reference)` and check row coverage and missing values. `provenance/validation_table_summary.csv` records comparisons with the saved manuscript/project tables.

## Completed checks

- In the initial release all 61 analysis/helper scripts parsed and all 113 input aliases resolved within the repository's `data/` folder, with the recorded file sizes and SHA-256 hashes. Registered workflows and script references exist.
- Serum LC-MS, intestinal LC-MS, caecum NMR, ex vivo NMR and combined-mode HILIC analyses ran from bundled inputs. The compared primary numerical results match the saved tables. Intestinal and HILIC curation workflows also ran from the included feature workbooks.
- Both stages of the primary RNA workflow ran from tximport and mouse metadata. The saved initial and outlier-excluded result tables matched, including the Liver, Colon and Brain primary contrasts. The original two-stage filtering order was preserved.
- Experimental Bray–Curtis values regenerated from the processed 16S counts match the saved Figure 3 table. Counts and taxonomy were also verified unchanged after removal of unused metadata fields.
- PCA preparation ran from the processed metabolomics inputs. The donor-restricted Procrustes test matched the saved result: 24 samples, correlation **0.8572254868783488**, permutation **P = 0.0001**, and 9,999 permutations.
- Integration input preparation, genus-metabolite/pathway association fits, leave-one-donor-out analyses and summaries ran from bundled experimental matrices and saved pathway tables. The compared statistical tables match. The RNA metadata bridge excludes C02 and G04, as in the original integration cache.
- Blood-flux meta-analysis ran from aggregate cohort statistics, without participant-level data. Its compared result tables matched.
- The supplementary expression preparation preserves all 123 source-table rows, including genes without a mapped symbol. All compared columns match.
- The five-workbook reporting workflow ran. All **41** assembled sheets from Supplementary Data 1–4 matched their reference snapshots, including **3,357,486 numeric cells**. Supplementary Data 5 was copied byte-for-byte.
- A complete repository copy was moved to another directory containing spaces. The input check, serum refit, Procrustes test and Bray–Curtis export all ran successfully there.
- Restricted cohort scripts stop with a clear data-unavailability message when no separately authorised private data directory is supplied. Bundled R objects were inspected for excluded metadata and fitted/model classes; serialized workspaces contain only the processed phyloseq object.

There were **407 complete table/worksheet matches** and **eight reporting differences** in the saved-project comparison. These counts include related exports and repeated curation/refit outputs, not 407 independent analyses. The reporting differences are:

- Three caecum-NMR workbook sheets: the current original script enables an additional adjusted sensitivity analysis (38 rows), whereas the saved reference's sensitivity tab is empty; the associated summary and input-path settings also differ. Primary results match. The original script's flag was retained.
- The ex vivo NMR settings sheet records a portable path.
- Two figure-source manifests record the new paths and compact RNA plotting input.
- Two integration input-audit exports record new paths, file sizes and timestamps.

These reporting differences are retained in the comparison summary; they are not labelled as complete output matches.

## Rendering and validation scope

All six main figures and the supplementary transcriptomics figure have now rendered successfully from the packaged inputs using a compatible isolated environment. The current manuscript 16S source script contained a panel-selection bug; the package corrects it and restores the separate published mouse/culture comparisons. Input checks, environment details, the plotting correction and visual limits are documented in `FIGURE_VALIDATION.md`. This figure check does not constitute independent revalidation of every manuscript claim. The manuscript figures themselves were not replaced.

The full annotation-dependent GO/GSEA/ssGSEA and gene-metabolite enrichment workflows were not independently rerun during packaging. Their recovered code and saved inputs are retained; further full enrichment/integration reruns are outside the agreed validation scope. Database and annotation versions can affect results, and saved reference tables are not independent recomputations. Historical additional mixed models are optional. Model construction is documented as an external Persephone dependency. No sex-stratified association or PPI-by-sex interaction results were found in the reviewed manuscript or Supplementary Data 1–5, so their exploratory blood-analysis branches are not required release checks.

Participant-level human cohort analyses cannot be rerun from this public package because their matrices and covariates cannot be provided. The matching primary MetaCardis model settings have now been reconstructed and validated locally; their aggregate meta-analysis inputs are bundled. Primary Kiel blood fits also matched the saved results from the original MAT inputs. These checks do not extend to a complete raw-to-paper model-construction rebuild. The recovered primary RNA limma output matches Supplementary Data 3; its fitting branch is included without claiming a fresh packaging refit.

## Reviewer checks

```sh
Rscript tests/check_inputs.R
Rscript run.R --workflow metabolomics-refit
Rscript tests/compare_generated_tables.R results
```

The second test compares generated CSV result tables that have an unambiguous bundled reference. It checks numeric values, missingness, row coverage and schema; it does not validate every workbook, figure or database-dependent enrichment. The local broader validation used the original saved project files as additional references.
