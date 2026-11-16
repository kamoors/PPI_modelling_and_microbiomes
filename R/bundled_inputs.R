# SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0
# Resolve manuscript input names to the processed files distributed in data/.
# This compatibility layer keeps statistical expressions independent of storage.

ppi_check_plot_versions <- function(packages) {
  # Earlier plot assemblers can silently omit legends with ggplot2 4.
  # Keep the plotting expressions and require compatible dependency versions.
  if (utils::packageVersion("ggplot2") < "4.0.0") return(invisible(NULL))
  minimum <- c(cowplot = "1.2.0", patchwork = "1.3.2")
  for (package in intersect(packages, names(minimum)))
    if (utils::packageVersion(package) < minimum[[package]])
      stop("With ggplot2 4, install ", package, " >= ", minimum[[package]],
           " before plotting. See environment/README.md.")
  invisible(NULL)
}

ppi_canonical_path <- function(path) {
  path <- gsub("\\\\", "/", path)
  aliases <- .ppi$aliases
  if (is.null(aliases) || !nrow(aliases)) return(path)
  for (i in order(nchar(aliases$legacy_prefix), decreasing = TRUE)) {
    old <- aliases$legacy_prefix[i]
    if (identical(path, old) || startsWith(path, paste0(old, "/")))
      return(paste0(aliases$package_prefix[i], substring(path, nchar(old) + 1L)))
  }
  path
}

ppi_initialize_legacy <- ppi_initialize
ppi_initialize <- function(code_root, data_root, results_root, external = list()) {
  ppi_initialize_legacy(code_root, data_root, results_root, external)
  .ppi$aliases <- read.csv(file.path(.ppi$code, "config", "path_aliases.csv"), stringsAsFactors = FALSE)
  manifest <- file.path(.ppi$data, "input_manifest.csv")
  .ppi$inputs <- if (file.exists(manifest)) read.csv(manifest, stringsAsFactors = FALSE) else
    data.frame(logical_path = character(), package_path = character())
  .ppi$inputs$logical_path <- vapply(.ppi$inputs$logical_path, ppi_canonical_path, character(1))
  if (any(grepl("^/|^[A-Za-z]:|(^|/)\\.\\.(/|$)", .ppi$inputs$package_path)))
    stop("Input manifest contains a path outside package data/.")
  .ppi$private_data <- Sys.getenv("PPI_PRIVATE_DATA_ROOT", unset = "")
  if (exists("private_data_root", inherits = TRUE) && nzchar(get("private_data_root", inherits = TRUE)))
    .ppi$private_data <- get("private_data_root", inherits = TRUE)
  .ppi$map$legacy_path <- vapply(.ppi$map$legacy_path, ppi_canonical_path, character(1))
  invisible(.ppi)
}

ppi_relative <- function(path) {
  path <- ppi_absolute(path)
  for (root in c(.ppi$work, .ppi$data)) {
    if (ppi_under(path, root)) {
      relative <- sub("^/+", "", substring(path, nchar(root) + 1L))
      return(ppi_canonical_path(gsub("/{2,}", "/", relative)))
    }
  }
  NA_character_
}
ppi_path <- function(relative = "") file.path(.ppi$work, ppi_canonical_path(relative))

ppi_private_input <- function(relative) {
  if (!nzchar(.ppi$private_data)) return(NULL)
  candidates <- file.path(.ppi$private_data, relative)
  for (i in seq_len(nrow(.ppi$aliases))) {
    prefix <- .ppi$aliases$package_prefix[i]
    if (identical(relative, prefix) || startsWith(relative, paste0(prefix, "/")))
      candidates <- c(candidates, file.path(.ppi$private_data,
        paste0(.ppi$aliases$legacy_prefix[i], substring(relative, nchar(prefix) + 1L))))
  }
  found <- candidates[file.exists(candidates)]
  if (length(found)) found[1L] else NULL
}

ppi_input <- function(path) {
  if (!is.character(path)) return(path)
  out <- vapply(path, function(value) {
    if (is.na(value) || !nzchar(value) || grepl("^[a-z]+://", value)) return(value)
    relative <- ppi_relative(value)
    if (is.na(relative)) return(ppi_absolute(value))
    generated <- if (exists("ppi_find_output", mode = "function")) ppi_find_output(relative) else NULL
    if (!is.null(generated)) return(generated)
    work <- ppi_path(relative)
    if (file.exists(work)) return(work)
    hit <- which(tolower(.ppi$inputs$logical_path) == tolower(relative))
    if (length(hit)) return(file.path(.ppi$data, .ppi$inputs$package_path[hit[1L]]))
    bundled <- file.path(.ppi$data, relative)
    if (file.exists(bundled)) return(bundled)
    private <- ppi_private_input(relative)
    if (!is.null(private)) return(private)
    bundled
  }, character(1), USE.NAMES = FALSE)
  names(out) <- names(path)
  out
}

ppi_dir_exists <- function(paths) vapply(paths, function(path) {
  if (dir.exists(ppi_input(path))) return(TRUE)
  relative <- ppi_relative(path)
  !is.na(relative) && any(startsWith(.ppi$inputs$logical_path, paste0(sub("/+$", "", relative), "/")))
}, logical(1), USE.NAMES = FALSE)

ppi_list_files <- function(path = ".", pattern = NULL, all.files = FALSE,
                           full.names = FALSE, recursive = FALSE,
                           ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE) {
  result <- character()
  for (one in path) {
    relative <- ppi_relative(one)
    locations <- if (is.na(relative)) ppi_input(one) else c(ppi_path(relative), file.path(.ppi$data, relative))
    found <- unique(unlist(lapply(locations, function(location)
      list.files(location, pattern = pattern, all.files = all.files, full.names = FALSE,
                 recursive = recursive, ignore.case = ignore.case,
                 include.dirs = include.dirs, no.. = no..)), use.names = FALSE))
    if (!is.na(relative)) {
      prefix <- if (nzchar(relative)) paste0(sub("/+$", "", relative), "/") else ""
      aliases <- .ppi$inputs$logical_path[startsWith(.ppi$inputs$logical_path, prefix)]
      suffix <- substring(aliases, nchar(prefix) + 1L)
      if (!recursive) suffix <- suffix[!grepl("/", suffix, fixed = TRUE)]
      if (!is.null(pattern)) suffix <- suffix[grepl(pattern, basename(suffix), ignore.case = ignore.case)]
      if (!all.files) suffix <- suffix[!startsWith(basename(suffix), ".")]
      found <- c(found, suffix, ppi_registered_listing(relative, recursive, pattern, ignore.case, all.files))
    }
    found <- sort(unique(found))
    if (full.names) found <- file.path(if (is.na(relative)) one else ppi_path(relative), found)
    result <- c(result, found)
  }
  unique(result)
}
