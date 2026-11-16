# SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0
# Run from any directory: Rscript /path/to/PPI_modelling_and_microbiomes/run.R --list
# Original metadata include non-ASCII identifiers; use a UTF-8 character locale.
if (.Platform$OS.type == "windows" && !isTRUE(l10n_info()[["UTF-8"]]))
  invisible(Sys.setlocale("LC_CTYPE", "English_United States.utf8"))
file_argument <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (length(file_argument) != 1L) stop("Run this entry point using Rscript.")
code_root <- dirname(normalizePath(sub("^--file=", "", file_argument), winslash = "/"))
args <- commandArgs(TRUE)
# Headless layout measurements need the same font-aware device as final PDFs.
# This affects graphics only; model fitting and exported numbers are unchanged.
if (!interactive() && capabilities("cairo"))
  options(device = function(...) grDevices::cairo_pdf(filename = tempfile(fileext = ".pdf"), ...))
value_after <- function(flag) {
  at <- match(flag, args)
  if (is.na(at)) return(NULL)
  if (at == length(args)) stop("Missing value after ", flag)
  args[at + 1L]
}
if (!is.null(value_after("--data-root"))) Sys.setenv(PPI_DATA_ROOT = value_after("--data-root"))
if (!is.null(value_after("--results-root"))) Sys.setenv(PPI_RESULTS_ROOT = value_after("--results-root"))
if (!is.null(value_after("--private-data-root"))) Sys.setenv(PPI_PRIVATE_DATA_ROOT = value_after("--private-data-root"))
source(file.path(code_root, "config", "config.R"))
source(file.path(code_root, "R", "paths.R"))
source(file.path(code_root, "R", "output_names.R"))
source(file.path(code_root, "R", "bundled_inputs.R"))
dir.create(data_root, recursive = TRUE, showWarnings = FALSE)
ppi_initialize(code_root, data_root, results_root, external)
index <- read.csv(file.path(code_root, "config", "entries.csv"), stringsAsFactors = FALSE)
if ("--list" %in% args) {
  print(index[, c("script", "role", "availability")], row.names = FALSE)

} else if ("--check" %in% args) {
  source(file.path(code_root, "R", "check_environment.R"))
} else {
  # Several original model scripts use %>% after requireNamespace-only checks.
  # Make that existing operator available without relying on an RStudio session.
  if (!exists("%>%", mode = "function") && requireNamespace("magrittr", quietly = TRUE))
    `%>%` <- magrittr::`%>%`
  workflow <- value_after("--workflow")
  if (!is.null(workflow)) {
    workflows <- read.csv(file.path(code_root, "config", "workflows.csv"), stringsAsFactors = FALSE)
    planned <- workflows$script[workflows$workflow == workflow]
    if (!length(planned)) stop("Unknown workflow. See docs/RUN_ORDER.md.")
    .ppi$plan <- index[match(planned, basename(index$script)), , drop = FALSE]
  } else {
    script <- value_after("--script")
    if (is.null(script)) stop("Use --list, --check, --workflow <name>, or --script <script>.")
    selected <- which(index$script == script | basename(index$script) == script)
    if (length(selected) != 1L) stop("Script is missing or ambiguous; use a path shown by --list.")
    .ppi$plan <- index[selected, , drop = FALSE]
  }
  stopifnot(!anyNA(.ppi$plan$script))
  for (.ppi_step in seq_len(nrow(.ppi$plan))) {
    if (.ppi$plan$availability[.ppi_step] == "restricted_human" && !nzchar(.ppi$private_data))
      stop("This human cohort analysis requires participant-level matrices and covariates that cannot be provided. ",
           "The scripts are included for transparency. See docs/RESTRICTED_DATA.md; ",
           "use --private-data-root only with separately authorised inputs.")
    if (.ppi$plan$availability[.ppi_step] == "optional_upstream")
      message("Optional upstream script: raw reads or metabolic models are not bundled. See docs/RUN_ORDER.md.")
    ppi_setwd(ppi_path(.ppi$plan$initial_directory[.ppi_step]))
    .ppi$active_group <- .ppi$plan$output_group[.ppi_step]
    .ppi$current_script <- basename(.ppi$plan$script[.ppi_step])
    message("Running ", .ppi$plan$script[.ppi_step])
    message("Reading bundled data from: ", .ppi$data)
    message("Writing results under: ", .ppi$results)
    source(file.path(.ppi$code, .ppi$plan$script[.ppi_step]), local = .GlobalEnv)
  }
  writeLines(capture.output(sessionInfo()), file.path(.ppi$results, "last_sessionInfo.txt"))
}
