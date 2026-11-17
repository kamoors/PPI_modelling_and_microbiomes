# SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0
# Serum raw-data loading

source("scripts/helpers/setup.R")
invisible(run_analysis(
  section = "04_mouse_metabolomics",
  directory = "metabolomics/Serum_metabolomics",
  script = "04_01_import_serum_lcms.R",
  restricted = FALSE,
  code = {
project_setwd(project_path("metabolomics/Serum_metabolomics/20260121_data"))
project_source(project_code("scripts/00_shared_helpers/00_02_import_serum_metabolite_sheet.R"))
project_source(project_code("scripts/00_shared_helpers/00_03_fit_serum_mixed_models.R"))
files <- project_list_files(project_path("metabolomics/Serum_metabolomics/20260121_data"),
    pattern = ".xlsx")
library(readxl)
f <- sapply(files, project_reader(read_excel))
pos_serum <- process_metabolite_sheet(data = f$`2026-01-21_20251205_AW_PPI_KIEL_STUDYII_MS_C8_RP_BA_pos_data.xlsx`,
    metid = f$`2026-01-21_20251205_AW_PPI_KIEL_STUDYII_MS_C8_RP_BA_pos_metid.xlsx`, meta_cols = 7)
neg_serum <- process_metabolite_sheet(data = f$`2026-01-21_20251202_AW_PPI_KIEL_STUDYII_MS_C8_RP_BA_neg_data.xlsx`,
    metid = f$`2026-01-21_20251202_AW_PPI_KIEL_STUDYII_MS_C8_RP_BA_neg_metid.xlsx`, meta_cols = 7)
}))
