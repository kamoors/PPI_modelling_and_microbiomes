# SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0
# Portable filesystem routing. Statistical expressions live in scripts/ unchanged.
ppi_initialize <- function(code_root, data_root, results_root, external = list()) {
  clean <- function(x) gsub("\\\\", "/", x)
  .ppi <<- new.env(parent = emptyenv())
  .ppi$code <- clean(normalizePath(code_root, winslash = "/", mustWork = TRUE))
  .ppi$data <- clean(normalizePath(data_root, winslash = "/", mustWork = TRUE))
  dir.create(results_root, recursive = TRUE, showWarnings = FALSE)
  .ppi$results <- clean(normalizePath(results_root, winslash = "/", mustWork = TRUE))
  if (identical(tolower(.ppi$data), tolower(.ppi$results)))
    stop("Data and results directories must be different.")
  .ppi$work <- file.path(.ppi$results, ".work")
  .ppi$external <- external
  .ppi$map <- read.csv(file.path(.ppi$code, "config", "script_paths.csv"),
                       stringsAsFactors = FALSE, check.names = FALSE)
  dir.create(.ppi$work, recursive = TRUE, showWarnings = FALSE)
  if (exists("ppi_initialize_output_names", mode = "function")) ppi_initialize_output_names()
  invisible(.ppi)
}

ppi_slashes <- function(x) gsub("\\\\", "/", x)
ppi_under <- function(path, root) {
  path <- sub("/+$", "", ppi_slashes(path))
  root <- sub("/+$", "", ppi_slashes(root))
  tolower(path) == tolower(root) | startsWith(tolower(path), paste0(tolower(root), "/"))
}
ppi_absolute <- function(x) {
  x <- ppi_slashes(path.expand(x))
  ifelse(grepl("^/|^[A-Za-z]:/", x), x, file.path(getwd(), x))
}
ppi_relative <- function(path) {
  path <- ppi_absolute(path)
  for (root in c(.ppi$work, .ppi$data)) {
    if (ppi_under(path, root)) {
      relative <- sub("^/+", "", substring(path, nchar(root) + 1L))
      return(gsub("/{2,}", "/", relative))
    }
  }
  NA_character_
}
ppi_path <- function(relative = "") file.path(.ppi$work, relative)
ppi_external <- function(key) {
  value <- .ppi$external[[key]]
  if (is.null(value) || !nzchar(value))
    stop("Configure external resource '", key, "' in config/local.R; see data/README.md.")
  ppi_absolute(value)
}
ppi_script <- function(path, required = FALSE) {
  if (length(path) != 1L || !is.character(path)) return(NULL)
  path <- ppi_slashes(path)
  if (ppi_under(ppi_absolute(path), .ppi$code) && file.exists(ppi_absolute(path)))
    return(ppi_absolute(path))
  relative <- ppi_relative(path)
  keys <- tolower(.ppi$map$legacy_path)
  hit <- if (!is.na(relative)) which(keys == tolower(relative)) else integer()
  if (!length(hit)) {
    candidates <- which(tolower(basename(.ppi$map$legacy_path)) == tolower(basename(path)))
    if (length(unique(.ppi$map$package_path[candidates])) == 1L) hit <- candidates[1L]
  }
  if (length(hit)) return(file.path(.ppi$code, .ppi$map$package_path[hit[1L]]))
  if (required) stop("No selected code file maps to source path: ", path)
  NULL
}
ppi_input <- function(path) {
  if (!is.character(path)) return(path)
  result <- vapply(path, function(value) {
    if (is.na(value) || !nzchar(value) || grepl("^[a-z]+://", value)) return(value)
    relative <- ppi_relative(value)
    if (!is.na(relative)) {
      renamed <- if (exists("ppi_find_output", mode = "function")) ppi_find_output(relative) else NULL
      if (!is.null(renamed)) return(renamed)
      generated <- ppi_path(relative)
      if (file.exists(generated)) return(generated)
      return(file.path(.ppi$data, relative))
    }
    ppi_absolute(value)
  }, character(1), USE.NAMES = FALSE)
  names(result) <- names(path)
  result
}
ppi_reader <- function(fun) {
  force(fun)
  function(file, ...) fun(ppi_input(file), ...)
}
ppi_output <- function(path, name = NULL, directory = FALSE) {
  if (!is.character(path)) return(path)
  result <- vapply(path, function(value) {
    if (is.na(value) || !nzchar(value)) return(value)
    absolute <- ppi_absolute(value)
    relative <- ppi_relative(absolute)
    if (!is.na(relative)) {
      if (isTRUE(directory)) dir.create(ppi_path(relative), recursive = TRUE, showWarnings = FALSE)
      target <- if (exists("ppi_named_output", mode = "function")) {
        ppi_named_output(relative, name = name, directory = directory)
      } else ppi_path(relative)
    } else if (ppi_under(absolute, tempdir())) {
      target <- absolute
    } else {
      # External ancillary figure exports cannot write into the input tree.
      target <- file.path(.ppi$results, "ancillary", basename(absolute))
    }
    dir.create(dirname(target), recursive = TRUE, showWarnings = FALSE)
    target
  }, character(1), USE.NAMES = FALSE)
  names(result) <- names(path)
  result
}
ppi_setwd <- function(dir) {
  relative <- ppi_relative(dir)
  if (is.na(relative)) stop("Working directory is outside the configured project: ", dir)
  target <- ppi_path(relative)
  dir.create(target, recursive = TRUE, showWarnings = FALSE)
  setwd(target)
}
ppi_file_exists <- function(...) {
  paths <- unlist(list(...), use.names = FALSE)
  vapply(paths, function(path) {
    script <- if (grepl("\\.[Rr]$", path)) ppi_script(path) else NULL
    if (!is.null(script)) return(file.exists(script))
    file.exists(ppi_input(path))
  }, logical(1), USE.NAMES = FALSE)
}
ppi_dir_exists <- function(paths) dir.exists(ppi_input(paths))
ppi_normalize_path <- function(path, winslash = "\\", mustWork = NA) {
  # Keep logical project paths for searches/prefix replacement in original code.
  vapply(path, function(value) {
    relative <- ppi_relative(value)
    if (!is.na(relative)) {
      if (isTRUE(mustWork) && !ppi_file_exists(value)) stop("Path does not exist: ", value)
      return(ppi_path(relative))
    }
    normalizePath(value, winslash = winslash, mustWork = mustWork)
  }, character(1), USE.NAMES = FALSE)
}
ppi_list_files <- function(path = ".", pattern = NULL, all.files = FALSE,
                           full.names = FALSE, recursive = FALSE,
                           ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE) {
  result <- character()
  for (one in path) {
    relative <- ppi_relative(one)
    locations <- if (is.na(relative)) ppi_input(one) else c(ppi_path(relative), file.path(.ppi$data, relative))
    names <- unique(unlist(lapply(locations, function(location) {
      list.files(location, pattern = pattern, all.files = all.files, full.names = FALSE,
                 recursive = recursive, ignore.case = ignore.case,
                 include.dirs = include.dirs, no.. = no..)
    }), use.names = FALSE))
    if (!is.na(relative) && exists("ppi_registered_listing", mode = "function")) {
      names <- c(names, ppi_registered_listing(relative, recursive, pattern, ignore.case, all.files))
    }
    names <- sort(unique(names))
    if (full.names) names <- file.path(if (is.na(relative)) one else ppi_path(relative), names)
    result <- c(result, names)
  }
  result
}
ppi_source <- function(file, local = FALSE, ...) {
  environment <- if (isTRUE(local)) parent.frame() else local
  base::source(ppi_script(file, required = TRUE), local = environment, ...)
}
ppi_sys_source <- function(file, envir = baseenv(), ...) {
  base::sys.source(ppi_script(file, required = TRUE), envir = envir, ...)
}
