# Explicit opt-in installation; analyses never install packages themselves.
argument <- grep("^--file=", commandArgs(FALSE), value = TRUE)
root <- dirname(dirname(normalizePath(sub("^--file=", "", argument), winslash = "/")))
dependencies <- read.csv(file.path(root, "environment", "packages.csv"), stringsAsFactors = FALSE)
available <- vapply(dependencies$package, requireNamespace, logical(1), quietly = TRUE)
cran <- dependencies$package[!available & dependencies$repository == "CRAN"]
bioc <- dependencies$package[!available & dependencies$repository == "Bioconductor"]
if (length(cran)) install.packages(cran, repos = "https://cloud.r-project.org")
if (length(bioc)) {
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager", repos = "https://cloud.r-project.org")
  BiocManager::install(version = "3.22", ask = FALSE, update = FALSE)
  BiocManager::install(bioc, ask = FALSE, update = FALSE)
}
cat("Installed missing available packages. Compare installed versions with packages.csv before full runs.\n")
