# SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0
# Reviewer-facing export names. Logical input references remain backwards compatible.
ppi_initialize_output_names <- function() {
  names_file <- file.path(.ppi$code, "config", "output_names.csv")
  .ppi$named_outputs <- file.exists(names_file)
  if (!.ppi$named_outputs) return(invisible(NULL))
  .ppi$output_names <- read.csv(names_file, stringsAsFactors = FALSE)
  .ppi$registry_file <- file.path(.ppi$results, "output_manifest.csv")
  .ppi$output_registry <- if (file.exists(.ppi$registry_file)) {
    read.csv(.ppi$registry_file, stringsAsFactors = FALSE)
  } else data.frame(logical_path = character(), output_path = character(),
                    producer = character(), section = character(), stringsAsFactors = FALSE)
  .ppi$active_group <- "00_shared_helpers"
  .ppi$current_script <- "shared_helpers"
}

ppi_code <- function(relative) file.path(.ppi$code, relative)

ppi_snake_case <- function(value) {
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

ppi_reviewer_filename <- function(filename) {
  override <- match(filename, .ppi$output_names$original_filename)
  if (!is.na(override)) return(.ppi$output_names$reviewer_filename[override])
  # FASTQ/FASTA names carry sample/model identities and pairing information.
  if (grepl("\\.(fastq|fq|fasta|fa)(\\.gz)?$", filename, ignore.case = TRUE)) return(filename)
  extension <- if (grepl("\\.(csv|tsv|txt)\\.gz$", filename, ignore.case = TRUE)) {
    sub("^.*\\.((csv|tsv|txt)\\.gz)$", "\\1", filename, ignore.case = TRUE)
  } else tools::file_ext(filename)
  stem <- if (nzchar(extension)) substr(filename, 1L, nchar(filename) - nchar(extension) - 1L) else filename
  stem <- ppi_snake_case(stem)
  if (!nzchar(stem)) stem <- "analysis_output"
  if (grepl("reaction_abundance", .ppi$current_script, fixed = TRUE))
    stem <- gsub("_clr(_|$)", "_", stem)
  stem <- sub("_$", "", stem)
  extension <- if (tolower(extension) == "rdata") "RData" else tolower(extension)
  paste0(stem, if (nzchar(extension)) paste0(".", extension) else "")
}

ppi_find_output <- function(relative) {
  if (!isTRUE(.ppi$named_outputs)) return(NULL)
  hit <- match(tolower(relative), tolower(.ppi$output_registry$logical_path))
  if (is.na(hit)) return(NULL)
  target <- file.path(.ppi$results, .ppi$output_registry$output_path[hit])
  if (file.exists(target)) target else NULL
}

ppi_named_output <- function(relative, name = NULL, directory = FALSE) {
  if (!isTRUE(.ppi$named_outputs)) return(ppi_path(relative))
  group <- .ppi$active_group
  if (isTRUE(directory)) return(file.path(.ppi$results, group))
  existing <- match(tolower(relative), tolower(.ppi$output_registry$logical_path))
  if (!is.na(existing)) return(file.path(.ppi$results, .ppi$output_registry$output_path[existing]))
  filename <- if (is.null(name)) ppi_reviewer_filename(basename(relative)) else name
  extension <- tolower(tools::file_ext(filename))
  category <- if (grepl("(^|/)source_data/", relative)) "figure_inputs" else
    if (extension %in% c("png", "pdf", "tiff", "tif", "svg", "jpeg", "jpg")) "figures" else
      if (extension %in% c("rds", "rdata")) "objects" else
        if (grepl("\\.(fastq|fq|fasta|fa)(\\.gz)?$", filename, ignore.case = TRUE)) "reads" else
        if (extension %in% c("csv", "tsv", "xlsx", "xls", "txt") ||
            grepl("\\.(csv|tsv|txt)\\.gz$", filename, ignore.case = TRUE)) "tables" else "other"
  proposed <- file.path(group, category, filename)
  if (proposed %in% .ppi$output_registry$output_path) {
    # Keep same-named exports from distinct logical folders separate.
    parent_label <- ppi_snake_case(basename(dirname(relative)))
    proposed <- file.path(group, category, parent_label, filename)
    if (proposed %in% .ppi$output_registry$output_path)
      stop("Output-name collision for ", relative, "; add an explicit distinct name.")
  }
  record <- data.frame(logical_path = relative, output_path = ppi_slashes(proposed),
                       producer = .ppi$current_script, section = group, stringsAsFactors = FALSE)
  .ppi$output_registry <- rbind(.ppi$output_registry, record)
  write.csv(.ppi$output_registry, .ppi$registry_file, row.names = FALSE)
  file.path(.ppi$results, proposed)
}

ppi_registered_listing <- function(relative, recursive, pattern, ignore.case, all.files) {
  if (!isTRUE(.ppi$named_outputs) || !nrow(.ppi$output_registry)) return(character())
  prefix <- if (nzchar(relative)) paste0(sub("/+$", "", relative), "/") else ""
  paths <- .ppi$output_registry$logical_path
  keep <- startsWith(paths, prefix)
  suffix <- substring(paths, nchar(prefix) + 1L)
  if (!recursive) keep <- keep & !grepl("/", suffix, fixed = TRUE)
  if (!is.null(pattern)) keep <- keep & grepl(pattern, basename(suffix), ignore.case = ignore.case)
  if (!all.files) keep <- keep & !startsWith(basename(suffix), ".")
  keep <- keep & file.exists(file.path(.ppi$results, .ppi$output_registry$output_path))
  suffix[keep]
}
