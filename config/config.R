# SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0
# Bundled processed inputs are the default. No original project tree is required.
data_root <- Sys.getenv("PPI_DATA_ROOT", unset = file.path(code_root, "data"))
results_root <- Sys.getenv("PPI_RESULTS_ROOT", unset = file.path(code_root, "results"))
private_data_root <- Sys.getenv("PPI_PRIVATE_DATA_ROOT", unset = "")
external <- list(
  agora2_models = Sys.getenv("PPI_AGORA2_MODELS", unset = "")
)
# Optional private overrides; this file is deliberately absent from the archive.
local_config <- file.path(code_root, "config", "local.R")
if (file.exists(local_config)) source(local_config, local = TRUE)
