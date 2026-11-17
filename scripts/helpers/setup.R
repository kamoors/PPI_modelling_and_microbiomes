# SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0
# Shared file handling for the numbered analysis scripts.

project_slashes <- function(x) gsub("\\\\", "/", x)

project_under <- function(path, root) {
    path <- sub("/+$", "", project_slashes(path))
    root <- sub("/+$", "", project_slashes(root))
    tolower(path) == tolower(root) | startsWith(tolower(path), paste0(tolower(root), "/"))
}

project_absolute <- function(x) {
    x <- project_slashes(path.expand(x))
    ifelse(grepl("^/|^[A-Za-z]:/", x), x, file.path(getwd(), x))
}

project_relative <- function(path) {
    path <- project_absolute(path)
    for (root in c(.project$work, .project$data)) {
        if (project_under(path, root)) {
            relative <- sub("^/+", "", substring(path, nchar(root) + 1L))
            return(project_canonical_path(gsub("/{2,}", "/", relative)))
        }
    }
    NA_character_
}

project_path <- function(relative = "") file.path(.project$work, project_canonical_path(relative))

project_external <- function(key) {
    value <- .project$external[[key]]
    if (is.null(value) || !nzchar(value))
        stop("Configure external resource '", key, "' using PPI_AGORA2_MODELS; see README.md.")
    project_absolute(value)
}

project_script <- function(path, required = FALSE) {
    if (length(path) != 1L || !is.character(path))
        return(NULL)
    path <- project_slashes(path)
    if (project_under(project_absolute(path), .project$code) && file.exists(project_absolute(path)))
        return(project_absolute(path))
    relative <- project_relative(path)
    keys <- tolower(.project$map$legacy_path)
    hit <- if (!is.na(relative))
        which(keys == tolower(relative))
    else integer()
    if (!length(hit)) {
        candidates <- which(tolower(basename(.project$map$legacy_path)) == tolower(basename(path)))
        if (length(unique(.project$map$package_path[candidates])) == 1L)
            hit <- candidates[1L]
    }
    if (length(hit))
        return(file.path(.project$code, .project$map$package_path[hit[1L]]))
    if (required)
        stop("No selected code file maps to source path: ", path)
    NULL
}

project_input <- function(path) {
    if (!is.character(path))
        return(path)
    out <- vapply(path, function(value) {
        if (is.na(value) || !nzchar(value) || grepl("^[a-z]+://", value))
            return(value)
        relative <- project_relative(value)
        if (is.na(relative))
            return(project_absolute(value))
        generated <- if (exists("project_find_output", mode = "function"))
            project_find_output(relative)
        else NULL
        if (!is.null(generated))
            return(generated)
        work <- project_path(relative)
        if (file.exists(work))
            return(work)
        hit <- which(tolower(.project$inputs$logical_path) == tolower(relative))
        if (length(hit))
            return(file.path(.project$data, .project$inputs$package_path[hit[1L]]))
        bundled <- file.path(.project$data, relative)
        if (file.exists(bundled))
            return(bundled)
        private <- project_private_input(relative)
        if (!is.null(private))
            return(private)
        bundled
    }, character(1), USE.NAMES = FALSE)
    names(out) <- names(path)
    out
}

project_reader <- function(fun) {
    force(fun)
    function(file, ...) fun(project_input(file), ...)
}

project_output <- function(path, name = NULL, directory = FALSE) {
    if (!is.character(path))
        return(path)
    result <- vapply(path, function(value) {
        if (is.na(value) || !nzchar(value))
            return(value)
        absolute <- project_absolute(value)
        relative <- project_relative(absolute)
        if (!is.na(relative)) {
            if (isTRUE(directory))
                dir.create(project_path(relative), recursive = TRUE, showWarnings = FALSE)
            target <- if (exists("project_named_output", mode = "function")) {
                project_named_output(relative, name = name, directory = directory)
            }
            else project_path(relative)
        }
        else if (project_under(absolute, tempdir())) {
            target <- absolute
        }
        else {
            target <- file.path(.project$results, "ancillary", basename(absolute))
        }
        dir.create(dirname(target), recursive = TRUE, showWarnings = FALSE)
        target
    }, character(1), USE.NAMES = FALSE)
    names(result) <- names(path)
    result
}

project_setwd <- function(dir) {
    relative <- project_relative(dir)
    if (is.na(relative))
        stop("Working directory is outside the configured project: ", dir)
    target <- project_path(relative)
    dir.create(target, recursive = TRUE, showWarnings = FALSE)
    setwd(target)
}

project_file_exists <- function(...) {
    paths <- unlist(list(...), use.names = FALSE)
    vapply(paths, function(path) {
        script <- if (grepl("\\.[Rr]$", path))
            project_script(path)
        else NULL
        if (!is.null(script))
            return(file.exists(script))
        file.exists(project_input(path))
    }, logical(1), USE.NAMES = FALSE)
}

project_dir_exists <- function(paths) vapply(paths, function(path) {
    if (dir.exists(project_input(path)))
        return(TRUE)
    relative <- project_relative(path)
    !is.na(relative) && any(startsWith(.project$inputs$logical_path, paste0(sub("/+$", "", relative), "/")))
}, logical(1), USE.NAMES = FALSE)

project_normalize_path <- function(path, winslash = "\\", mustWork = NA) {
    vapply(path, function(value) {
        relative <- project_relative(value)
        if (!is.na(relative)) {
            if (isTRUE(mustWork) && !project_file_exists(value))
                stop("Path does not exist: ", value)
            return(project_path(relative))
        }
        normalizePath(value, winslash = winslash, mustWork = mustWork)
    }, character(1), USE.NAMES = FALSE)
}

project_list_files <- function(path = ".", pattern = NULL, all.files = FALSE, full.names = FALSE, recursive = FALSE,
    ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE) {
    result <- character()
    for (one in path) {
        relative <- project_relative(one)
        locations <- if (is.na(relative))
            project_input(one)
        else c(project_path(relative), file.path(.project$data, relative))
        found <- unique(unlist(lapply(locations, function(location) list.files(location, pattern = pattern,
            all.files = all.files, full.names = FALSE, recursive = recursive, ignore.case = ignore.case,
            include.dirs = include.dirs, no.. = no..)), use.names = FALSE))
        if (!is.na(relative)) {
            prefix <- if (nzchar(relative))
                paste0(sub("/+$", "", relative), "/")
            else ""
            aliases <- .project$inputs$logical_path[startsWith(.project$inputs$logical_path, prefix)]
            suffix <- substring(aliases, nchar(prefix) + 1L)
            if (!recursive)
                suffix <- suffix[!grepl("/", suffix, fixed = TRUE)]
            if (!is.null(pattern))
                suffix <- suffix[grepl(pattern, basename(suffix), ignore.case = ignore.case)]
            if (!all.files)
                suffix <- suffix[!startsWith(basename(suffix), ".")]
            found <- c(found, suffix, project_registered_listing(relative, recursive, pattern, ignore.case,
                all.files))
        }
        found <- sort(unique(found))
        if (full.names)
            found <- file.path(if (is.na(relative))
                one
            else project_path(relative), found)
        result <- c(result, found)
    }
    unique(result)
}

project_source <- function(file, local = FALSE, ...) {
    environment <- if (isTRUE(local))
        parent.frame()
    else local
    base::source(project_script(file, required = TRUE), local = environment, ...)
}

project_sys_source <- function(file, envir = baseenv(), ...) {
    base::sys.source(project_script(file, required = TRUE), envir = envir, ...)
}

project_initialize_output_names <- function() {
    names_file <- file.path(.project$code, "scripts", "helpers", "output_names.csv")
    .project$named_outputs <- file.exists(names_file)
    if (!.project$named_outputs)
        return(invisible(NULL))
    .project$output_names <- read.csv(names_file, stringsAsFactors = FALSE)
    .project$registry_file <- file.path(.project$results, "files.csv")
    .project$output_registry <- if (file.exists(.project$registry_file)) {
        read.csv(.project$registry_file, stringsAsFactors = FALSE)
    }
    else data.frame(logical_path = character(), output_path = character(), producer = character(), section = character(),
        stringsAsFactors = FALSE)
    .project$active_group <- "00_shared_helpers"
    .project$current_script <- "shared_helpers"
}

project_code <- function(relative) file.path(.project$code, relative)

project_snake_case <- function(value) {
    value <- gsub("([A-Z]+)([A-Z][a-z])", "\\1_\\2", value)
    value <- gsub("([a-z0-9])([A-Z])", "\\1_\\2", value)
    value <- tolower(gsub("[^A-Za-z0-9]+", "_", value))
    value <- gsub("(^|_)(final|newppi|with_new_ppi|new_ppi|v[0-9]+)(_|$)", "_", value)
    value <- gsub("focus", "kiel", value, fixed = TRUE)
    value <- gsub("cecum", "caecum", value, fixed = TRUE)
    value <- gsub("feces", "faeces", value, fixed = TRUE)
    value <- gsub("dup_cor", "duplicate_correlation", value, fixed = TRUE)
    value <- gsub("_{2,}", "_", value)
    gsub("^_|_$", "", value)
}

project_reviewer_filename <- function(filename) {
    override <- match(filename, .project$output_names$original_filename)
    if (!is.na(override))
        return(.project$output_names$reviewer_filename[override])
    if (grepl("\\.(fastq|fq|fasta|fa)(\\.gz)?$", filename, ignore.case = TRUE))
        return(filename)
    extension <- if (grepl("\\.(csv|tsv|txt)\\.gz$", filename, ignore.case = TRUE)) {
        sub("^.*\\.((csv|tsv|txt)\\.gz)$", "\\1", filename, ignore.case = TRUE)
    }
    else tools::file_ext(filename)
    stem <- if (nzchar(extension))
        substr(filename, 1L, nchar(filename) - nchar(extension) - 1L)
    else filename
    stem <- project_snake_case(stem)
    if (!nzchar(stem))
        stem <- "analysis_output"
    if (grepl("reaction_abundance", .project$current_script, fixed = TRUE))
        stem <- gsub("_clr(_|$)", "_", stem)
    stem <- sub("_$", "", stem)
    extension <- if (tolower(extension) == "rdata")
        "RData"
    else tolower(extension)
    paste0(stem, if (nzchar(extension))
        paste0(".", extension)
    else "")
}

project_find_output <- function(relative) {
    if (!isTRUE(.project$named_outputs))
        return(NULL)
    hit <- match(tolower(relative), tolower(.project$output_registry$logical_path))
    if (is.na(hit))
        return(NULL)
    target <- file.path(.project$results, .project$output_registry$output_path[hit])
    if (file.exists(target))
        target
    else NULL
}

project_named_output <- function(relative, name = NULL, directory = FALSE) {
    if (!isTRUE(.project$named_outputs))
        return(project_path(relative))
    group <- .project$active_group
    if (isTRUE(directory))
        return(file.path(.project$results, group))
    existing <- match(tolower(relative), tolower(.project$output_registry$logical_path))
    if (!is.na(existing))
        return(file.path(.project$results, .project$output_registry$output_path[existing]))
    filename <- if (is.null(name))
        project_reviewer_filename(basename(relative))
    else name
    extension <- tolower(tools::file_ext(filename))
    category <- if (grepl("(^|/)source_data/", relative))
        "figure_inputs"
    else if (extension %in% c("png", "pdf", "tiff", "tif", "svg", "jpeg", "jpg"))
        "figures"
    else if (extension %in% c("rds", "rdata"))
        "objects"
    else if (grepl("\\.(fastq|fq|fasta|fa)(\\.gz)?$", filename, ignore.case = TRUE))
        "reads"
    else if (extension %in% c("csv", "tsv", "xlsx", "xls", "txt") || grepl("\\.(csv|tsv|txt)\\.gz$",
        filename, ignore.case = TRUE))
        "tables"
    else "other"
    proposed <- file.path(group, category, filename)
    if (proposed %in% .project$output_registry$output_path) {
        parent_label <- project_snake_case(basename(dirname(relative)))
        proposed <- file.path(group, category, parent_label, filename)
        if (proposed %in% .project$output_registry$output_path)
            stop("Output-name collision for ", relative, "; add an explicit distinct name.")
    }
    record <- data.frame(logical_path = relative, output_path = project_slashes(proposed), producer = .project$current_script,
        section = group, stringsAsFactors = FALSE)
    .project$output_registry <- rbind(.project$output_registry, record)
    write.csv(.project$output_registry, .project$registry_file, row.names = FALSE)
    file.path(.project$results, proposed)
}

project_registered_listing <- function(relative, recursive, pattern, ignore.case, all.files) {
    if (!isTRUE(.project$named_outputs) || !nrow(.project$output_registry))
        return(character())
    prefix <- if (nzchar(relative))
        paste0(sub("/+$", "", relative), "/")
    else ""
    paths <- .project$output_registry$logical_path
    keep <- startsWith(paths, prefix)
    suffix <- substring(paths, nchar(prefix) + 1L)
    if (!recursive)
        keep <- keep & !grepl("/", suffix, fixed = TRUE)
    if (!is.null(pattern))
        keep <- keep & grepl(pattern, basename(suffix), ignore.case = ignore.case)
    if (!all.files)
        keep <- keep & !startsWith(basename(suffix), ".")
    keep <- keep & file.exists(file.path(.project$results, .project$output_registry$output_path))
    suffix[keep]
}

project_check_plot_versions <- function(packages) {
    if (utils::packageVersion("ggplot2") < "4.0.0")
        return(invisible(NULL))
    minimum <- c(cowplot = "1.2.0", patchwork = "1.3.2")
    for (package in intersect(packages, names(minimum))) if (utils::packageVersion(package) < minimum[[package]])
        stop("With ggplot2 4, install ", package, " >= ", minimum[[package]], " before plotting. See README.md.")
    invisible(NULL)
}

project_canonical_path <- function(path) {
    path <- gsub("\\\\", "/", path)
    aliases <- .project$aliases
    if (is.null(aliases) || !nrow(aliases))
        return(path)
    for (i in order(nchar(aliases$legacy_prefix), decreasing = TRUE)) {
        old <- aliases$legacy_prefix[i]
        if (identical(path, old) || startsWith(path, paste0(old, "/")))
            return(paste0(aliases$package_prefix[i], substring(path, nchar(old) + 1L)))
    }
    path
}

project_private_input <- function(relative) {
    if (!nzchar(.project$private_data))
        return(NULL)
    candidates <- file.path(.project$private_data, relative)
    for (i in seq_len(nrow(.project$aliases))) {
        prefix <- .project$aliases$package_prefix[i]
        if (identical(relative, prefix) || startsWith(relative, paste0(prefix, "/")))
            candidates <- c(candidates, file.path(.project$private_data, paste0(.project$aliases$legacy_prefix[i],
                substring(relative, nchar(prefix) + 1L))))
    }
    found <- candidates[file.exists(candidates)]
    if (length(found))
        found[1L]
    else NULL
}


# Resolve the stored input names and keep generated files under results/.
project_initialize <- function(root) {
  .project <<- new.env(parent = emptyenv())
  .project$code <- normalizePath(root, winslash = "/", mustWork = TRUE)
  data_root <- Sys.getenv("PPI_DATA_ROOT", file.path(root, "data"))
  results_root <- Sys.getenv("PPI_RESULTS_ROOT", file.path(root, "results"))
  .project$data <- normalizePath(data_root, winslash = "/", mustWork = TRUE)
  dir.create(results_root, recursive = TRUE, showWarnings = FALSE)
  .project$results <- normalizePath(results_root, winslash = "/", mustWork = TRUE)
  if (identical(tolower(.project$data), tolower(.project$results)))
    stop("Data and results directories must be different.")
  .project$work <- file.path(.project$results, ".work")
  dir.create(.project$work, recursive = TRUE, showWarnings = FALSE)
  .project$external <- list(agora2_models = Sys.getenv("PPI_AGORA2_MODELS", ""))
  .project$aliases <- data.frame(legacy_prefix = c("Validation/16S_results/PPI_rawdata/PPI_rawdata", "Validation/16S_results/PPI_rawdata", "Validation/Munich Metabolomics", "Validation/Mouse_RNASeq", "Validation/CrossLayer_Correlations", "Manuscript/Figures", "Results/Kiel_Metacardis", "Results/Blood_FBA", "Bram/transfer_130444_files_9d5afaf0/mgpipRes/totalresults", "Metacardis/data/MetaCardis", "Metacardis/data/WBM", "data/P2N_2022-006_Karlis_Moors", "Metacardis/Results", "Metacardis/data/WBM/Blood_FBA_Kiel"),
    package_prefix = c("microbiome16s", "microbiome16s", "metabolomics", "rnaseq", "integration", "figures", "human_cohorts/results", "human_cohorts/blood_flux/kiel", "human_cohorts/kiel", "human_cohorts/metacardis", "human_cohorts/blood_flux/metacardis", "human_cohorts/kiel", "human_cohorts/metacardis_results", "human_cohorts/blood_flux/kiel_models"), stringsAsFactors = FALSE)
  .project$inputs <- read.csv(file.path(.project$data, "files.csv"), stringsAsFactors = FALSE)
  .project$inputs$logical_path <- vapply(.project$inputs$logical_path, project_canonical_path, character(1))
  if (any(grepl("^/|^[A-Za-z]:|(^|/)\\.\\.(/|$)", .project$inputs$package_path)))
    stop("An input path points outside data/.")
  .project$map <- read.csv(file.path(root, "scripts/helpers/source_files.csv"), stringsAsFactors = FALSE)
  .project$map$legacy_path <- vapply(.project$map$legacy_path, project_canonical_path, character(1))
  .project$private_data <- Sys.getenv("PPI_PRIVATE_DATA_ROOT", "")
  if (exists("private_data_root", inherits = TRUE) && nzchar(get("private_data_root", inherits = TRUE)))
    .project$private_data <- get("private_data_root", inherits = TRUE)
  project_initialize_output_names()
  invisible(.project)
}

# Evaluate scripts in the usual source() environment; restore the user's folder
# on completion or error. Model calculations are in the numbered scripts.
run_analysis <- function(section, directory, script, code, restricted = FALSE) {
  old_wd <- getwd()
  if (!file.exists(file.path(old_wd, "data/files.csv")) ||
      !file.exists(file.path(old_wd, "scripts/helpers/setup.R")))
    stop("Set the working directory to the repository folder first. See README.md.")
  on.exit(setwd(old_wd), add = TRUE)
  if (.Platform$OS.type == "windows" && !isTRUE(l10n_info()[["UTF-8"]]))
    invisible(Sys.setlocale("LC_CTYPE", "English_United States.utf8"))
  if (!interactive() && capabilities("cairo"))
    options(device = function(...) grDevices::cairo_pdf(filename = tempfile(fileext = ".pdf"), ...))
  project_initialize(old_wd)
  .project$active_group <- section
  .project$current_script <- script
  if (restricted && !nzchar(.project$private_data))
    stop("This analysis needs participant-level human cohort data and covariates, which cannot be shared. ",
         "See README.md. Aggregate tables and figure inputs are provided separately.")
  if (!exists("%>%", mode = "function") && requireNamespace("magrittr", quietly = TRUE))
    assign("%>%", magrittr::`%>%`, envir = .GlobalEnv)
  project_setwd(project_path(directory))
  message("Running ", script, "; results: ", file.path(.project$results, section))
  eval(substitute(code), envir = .GlobalEnv)
  invisible(NULL)
}
