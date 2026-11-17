# SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0
# DADA2/phyloseq preprocessing
# Requires raw reads and taxonomy references, which are not included here.
# Script 03_02 starts from the supplied processed phyloseq object.

source("scripts/helpers/setup.R")
invisible(run_analysis(
  section = "03_mouse_microbiome",
  directory = "microbiome16s",
  script = "03_01_process_mouse_16s.R",
  restricted = FALSE,
  code = {
library(dada2)
library(phyloseq)
library(Biostrings)
library(ggplot2)
path <- project_path("microbiome16s")
project_list_files(path)
fnFs <- sort(project_list_files(path, pattern = "_R1_001.fastq.gz$", full.names = TRUE))
fnRs <- sort(project_list_files(path, pattern = "_R2_001.fastq.gz$", full.names = TRUE))
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
plotQualityProfile(project_input(fnFs[c(1:2, 40:41)]))
plotQualityProfile(project_input(fnRs[c(1:2, 40:41)]))
sample.names.unique <- make.unique(sample.names)
filtFs <- file.path(path, "filtered", paste0(sample.names.unique, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names.unique, "_R_filt.fastq.gz"))
out <- filterAndTrim(project_input(fnFs), project_output(filtFs), project_input(fnRs), project_output(filtRs), truncLen = c(240,
    190), maxN = 0, maxEE = c(2, 2), truncQ = 2, rm.phix = TRUE, compress = TRUE, multithread = TRUE)
errF <- learnErrors(project_input(filtFs), multithread = TRUE)
errR <- learnErrors(project_input(filtRs), multithread = TRUE)
plotErrors(errF, nominalQ = TRUE)
derepFs <- derepFastq(project_input(filtFs))
derepRs <- derepFastq(project_input(filtRs))
names(derepFs) <- sample.names.unique
names(derepRs) <- sample.names.unique
dadaFs <- dada(derepFs, err = errF, multithread = TRUE)
dadaRs <- dada(derepRs, err = errR, multithread = TRUE)
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose = TRUE)
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
seqtab.nochim <- removeBimeraDenovo(seqtab, method = "consensus", multithread = TRUE, verbose = TRUE)
dim(seqtab.nochim)
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names.unique
track <- as.data.frame(track)
track$ratio <- track$nonchim/track$merged
taxa_og <- assignTaxonomy(seqtab.nochim, project_input(project_path("microbiome16s/silva_nr99_v138.2_toSpecies_trainset.fa.gz")),
    multithread = TRUE)
taxa_og <- addSpecies(taxa, project_input(project_path("microbiome16s/silva_v138.2_assignSpecies.fa.gz")))
taxa <- taxa_og[, -7]
unqs.mock <- seqtab.nochim["MOCK", ]
unqs.mock <- sort(unqs.mock[unqs.mock > 0], decreasing = TRUE)
cat("DADA2 inferred", length(unqs.mock), "sample sequences present in the Mock community.\n")
library(phyloseq)
library(Biostrings)
library(ggplot2)
theme_set(theme_bw())
metadata <- readxl::read_excel(project_input(project_path("Validation/16S_sample_metadata.xlsx")))
rownames(metadata) <- metadata$barcode
metadata$experiment <- 0
metadata$experiment[grep(x = metadata$sample_name, pattern = "_CP_")] <- "culture"
metadata$experiment[grep(x = metadata$sample_name, pattern = "OG_")] <- "OG_sample"
metadata$experiment[grep(x = metadata$sample_name, pattern = "_m_")] <- "mouse"
metadata$PPI <- 0
metadata$PPI[grep(x = metadata$subject_name, pattern = "PPI")] <- 1
metadata$PPI <- as.factor(metadata$PPI)
metadata$batch <- as.factor(metadata$batch)
rownames(metadata) <- as.character(metadata$barcode)
ind <- rownames(seqtab.nochim) %in% metadata$barcode
seqtab.nochim_filt <- seqtab.nochim[ind, ]
ind2 <- match(metadata$barcode, rownames(seqtab.nochim_filt))
seqtab.nochim_filt <- seqtab.nochim_filt[ind2, ]
print(rownames(seqtab.nochim_filt) == metadata$barcode)
otu <- otu_table(seqtab.nochim_filt, taxa_are_rows = FALSE)
tax <- tax_table(taxa)
samples <- sample_data(metadata)
rownames(samples) <- as.character(samples$barcode)
ps <- phyloseq(otu, tax, samples)
sample_names(ps)
dna <- DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
ps
preprocessing_out <- list(seqtab.nochim = seqtab.nochim, seqtab.nochim_filt = seqtab.nochim_filt, taxa = taxa,
    otu = otu, tax = tax, samples = samples, ps = ps)
saveRDS(preprocessing_out, project_output(project_path("microbiome16s/16S_PPI_preprocessing_withRepeats.rds"),
    name = "mouse_16s_processed_objects.rds"))
save(seqtab.nochim, seqtab.nochim_filt, taxa, otu, tax, samples, ps, metadata, file = project_output(project_path("microbiome16s/16S_PPI_preprocessing_withRepeats.RData"),
    name = "mouse_16s_processed_objects.RData"))
}))
