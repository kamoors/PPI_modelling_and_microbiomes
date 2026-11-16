# Model construction with Persephone

Model construction used the external Persephone pipeline from the Thiele lab. Its complete implementation is referenced rather than redistributed. See the [official tutorial](https://opencobra.github.io/cobratoolbox/stable/tutorials/analysis/persephone/TutorialFour_HumanMicrobiomeModels.html), [COBRA Toolbox source](https://github.com/opencobra/cobratoolbox), and Nap et al. (2026), *Personalization and evaluation of human whole-body metabolic models and community microbiome models using Persephone*, [Nature Protocols](https://doi.org/10.1038/s41596-026-01437-8).

The study-specific MATLAB wrapper is `scripts/02_human_metabolism/s02_00_construct_metacardis_wbm.m`. With the relevant software on the MATLAB path:

```matlab
addpath('scripts/02_human_metabolism');
s02_00_construct_metacardis_wbm(hostModelDir, microbiomeModelDir, metadataFile, outputDir);
```

The four arguments name separately supplied host-model files, microbiome models, participant metadata and an output folder. Expected host files are Harvey_1_03d.mat (variable male) and Harvetta_1_03d.mat (variable female). The wrapper extracts sample IDs from microbiota_model_diet_<ID>.mat, matches metadata ID/Gender, selects the corresponding host model and calls `createMWBM(modelM, modelH, 'EUAverageDiet')`. It saves modelHM under a sex-prefixed filename, skips existing outputs and retains the received handling of missing metadata and failed assembly calls. Only filesystem paths were changed into function arguments.

The received source has no additional explicit diet-changing call after assembly. The EUAverageDiet argument remains part of assembly; no assumption is made about undocumented historical defaults. Host-model filenames are recorded as supplied, not independently established release identifiers. The exact study Persephone/COBRA commit and MATLAB/solver versions remain unrecorded. No construction rerun or model-equivalence check is claimed. Study-specific flux objectives and bounds are separate from this assembly wrapper. External software retains its original licence.
