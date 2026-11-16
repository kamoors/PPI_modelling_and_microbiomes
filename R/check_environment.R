# SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0
dependencies <- read.csv(file.path(.ppi$code, "environment", "packages.csv"), stringsAsFactors = FALSE)
dependencies$available <- vapply(dependencies$package, function(package) {
  suppressWarnings(requireNamespace(package, quietly = TRUE))
}, logical(1))
dependencies$installed_version <- vapply(dependencies$package, function(package) {
  tryCatch(as.character(utils::packageVersion(package)), error = function(e) "MISSING")
}, character(1))
dependencies$matches_record <- dependencies$installed_version == dependencies$recorded_version
write.csv(dependencies, file.path(.ppi$results, "environment_check.csv"), row.names = FALSE)
cat("R:", R.version.string, "\n")
cat("Required package namespaces available:", sum(dependencies$available), "/", nrow(dependencies), "\n")
cat("External input root:", .ppi$data, "\nResults root:", .ppi$results, "\n")
if (any(!dependencies$available)) {
  cat("Missing/unloadable packages:", paste(dependencies$package[!dependencies$available], collapse = ", "), "\n")
}
cat("See environment_check.csv; availability alone does not validate inputs or numerical reproduction.\n")
