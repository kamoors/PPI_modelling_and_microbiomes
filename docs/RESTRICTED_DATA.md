# Human cohort data

Participant-level data from the Kiel and MetaCardis human cohorts cannot be provided in this repository. The excluded inputs include:

- Microbial taxonomic and reaction-abundance matrices.
- Predicted metabolite secretion and participant-level blood-flux matrices.
- Participant-level serum metabolomics.
- Demographic, clinical, medication, dietary and other adjustment covariates.

The available aggregate tables contain cohort-level statistics and annotations. Their presence does not permit reconstruction or refitting of participant-level models.

The cohort analysis scripts are included for transparency. The runner stops before these scripts execute unless the user explicitly supplies separately authorised inputs with `--private-data-root /path/to/private/project`. No request or transfer of restricted data is performed by the code. `config/script_paths.csv`, the script input statements and the input manifest document expected filenames. Private data must remain outside the repository.

The experimental 16S and ex vivo datasets are separate from the restricted cohort datasets and are included in processed form.

The included primary MetaCardis generators reconstruct the model settings validated against the final results: reaction abundance omits systolic blood pressure and includes log total reaction abundance; secretion retains systolic blood pressure. Both apply the same metabolomics-based non-user exclusion. The matching serum generator is retained unchanged. `cohort-meta` uses only aggregate statistics, so it can run without private data. The alternate same-named standalone secretion file from Results/Kiel_Metacardis is deliberately not used. See SOURCE_UPDATES_20261113.md.
