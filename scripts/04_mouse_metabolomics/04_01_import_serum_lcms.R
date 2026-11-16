# SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0
# Serum raw-data loading
# Run through run.R so input and output paths resolve within the package.
# Inputs: see data/input_manifest.csv and docs/RUN_ORDER.md.
# Writes generated files under results/; bundled inputs are read-only.

ppi_setwd(ppi_path("metabolomics/Serum_metabolomics/20260121_data"))
ppi_source(ppi_code("scripts/00_shared_helpers/00_02_import_serum_metabolite_sheet.R"))
ppi_source(ppi_code("scripts/00_shared_helpers/00_03_fit_serum_mixed_models.R"))
files <- ppi_list_files(ppi_path("metabolomics/Serum_metabolomics/20260121_data"), 
    pattern = ".xlsx")
library(readxl)
f <- sapply(files, ppi_reader(read_excel))
pos_serum <- process_metabolite_sheet(data = f$`2026-01-21_20251205_AW_PPI_KIEL_STUDYII_MS_C8_RP_BA_pos_data.xlsx`, 
    metid = f$`2026-01-21_20251205_AW_PPI_KIEL_STUDYII_MS_C8_RP_BA_pos_metid.xlsx`, meta_cols = 7)
neg_serum <- process_metabolite_sheet(data = f$`2026-01-21_20251202_AW_PPI_KIEL_STUDYII_MS_C8_RP_BA_neg_data.xlsx`, 
    metid = f$`2026-01-21_20251202_AW_PPI_KIEL_STUDYII_MS_C8_RP_BA_neg_metid.xlsx`, meta_cols = 7)
