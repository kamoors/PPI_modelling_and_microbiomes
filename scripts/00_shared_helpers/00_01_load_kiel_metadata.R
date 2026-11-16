# SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0
# Kiel metadata preparation
# Run through run.R so input and output paths resolve within the package.
# Inputs: see data/input_manifest.csv and docs/RUN_ORDER.md.
# Writes generated files under results/; bundled inputs are read-only.

require(readr)
require(phyloseq)
require(readxl)
require(dplyr)
ppi_setwd(ppi_path("human_cohorts/kiel"))
PPI_data <- read_xlsx(ppi_input("2016-012_POPGEN.AND.FOCUS_nondietary_data_rev_date_mod_1463S_6HL.class_2018-12-20.xlsx"))
full_disease_phenotypes <- read_csv(ppi_input("phenoData.csv"))
phenotype <- read_tsv(ppi_input("Phenotype.tsv"), col_names = TRUE)
nutrition <- read_csv(ppi_input("nutrition-diabetes-full-norm.csv"), col_names = TRUE)
phenotype_all_drugs <- read_excel(ppi_input("2016-012_POPGEN.AND.FOCUS_nondietary_data_rev_date_mod_1463S_6HL.class_2018-12-20.xlsx"), 
    col_names = TRUE)
phenotype_all_drugs$old_id <- phenotype$ID[match(phenotype_all_drugs$new_id, phenotype$new_id.y)]
phenotype_all_drugs_filt <- phenotype_all_drugs[match(phenotype$ID, phenotype_all_drugs$old_id), ]
phenotype_all_drugs_filt <- phenotype_all_drugs_filt %>% relocate(old_id)
rm(phenotype_all_drugs)
phenoData <- read_excel(ppi_input("phenoData.xlsx"))
phenoData_filt <- phenoData[match(phenotype$ID, phenoData$ID), ]
rm(phenoData)
coiunts <- read_csv(ppi_input("species-counts.csv"), col_names = TRUE)
counts = as.matrix(coiunts[-1, -1])
rownames(counts) <- coiunts[-1, 1]$...1
counts_norm = as.matrix(coiunts[-1, -1])
rownames(counts_norm) <- coiunts[-1, 1]$...1
counts_norm = t(scale(as.matrix(counts_norm), center = FALSE, scale = colSums(as.matrix(counts_norm))))
ind_norm <- match(phenotype$ID, rownames(counts_norm))
counts_norm_filt <- counts_norm[ind_norm, ]
ind_4 <- match(phenotype$ID, colnames(counts))
counts_filt <- counts[, ind_4]
reac_act <- readRDS(ppi_input("reakFreqsKielCohort.RDS"))
ind_4 <- match(phenotype$ID, colnames(reac_act))
reac_act_filt <- reac_act[, ind_4]
new_id_match <- match(phenotype$new_id.y, PPI_data$new_id)
id_match_nutr <- match(phenotype$ID, nutrition$...1)
phenotype$drugs <- PPI_data$medikamente.x[new_id_match]
phenotype$alcohol <- nutrition$`EX_etoh(e)`[id_match_nutr]
phenotype$water <- nutrition$`EX_h2o(e)`[id_match_nutr]
phenotype$gramm <- nutrition$gramm[id_match_nutr]
ppi_terms <- c("omeprazol", "esomeprazol", "pantoprazol", "lansoprazol", "rabeprazol", "dexlansoprazol", 
    "pantozol", "antra", "nexium", "nexium mups", "pantopra-q", "omep", "esomep", "patoprazol", "omepazol", 
    "omeparazol")
ppi_pattern <- paste(ppi_terms, collapse = "|")
phenotype$PPI <- as.integer(grepl(ppi_pattern, phenotype$drugs, ignore.case = TRUE))
phenotype_all_drugs_filt$old_id == phenotype$ID
phenotype_all_drugs_filt$PPI <- phenotype$PPI
ind <- match(phenotype$ID, full_disease_phenotypes$ID)
full_disease_phenotypes <- full_disease_phenotypes[ind, ]
rownames(phenotype) <- phenotype$ID
statin_terms <- c("simvastatin", "atorvastatin", "pravastatin", "fluvastatin", "lovastatin", "rosuvastatin", 
    "simvabeta", "simvahexal", "simva hexa", "simva-aristo", "simvalip", "simvastin", "simvastatin al", 
    "simvastatin-ratiopharm", "simvastatin stada", "simvastatin qph", "simvastatin-1a pharma", "atorvastatin hexal", 
    "atorvastatin stada", "pravastatin abz", "prava-q", "fluvastatin actavis", "simvasstatin", "simvstatin")
aspirin_terms <- c("aspirin", "\\bASS\\b", "acetylsalicyl", "aspirin protect", "aspirin direkt", "aspirin c", 
    "aspirin komplex", "godamed")
antiplatelet_terms <- c("clopidogrel", "plavix", "clopidogel", "clopidogrll")
statin_pattern <- paste(statin_terms, collapse = "|")
aspirin_pattern <- paste(aspirin_terms, collapse = "|")
antiplatelet_pattern <- paste(antiplatelet_terms, collapse = "|")
phenotype$STATIN <- as.integer(grepl(statin_pattern, phenotype$drugs, ignore.case = TRUE))
phenotype$ASPIRIN <- as.integer(grepl(aspirin_pattern, phenotype$drugs, ignore.case = TRUE))
phenotype$ANTIPLATELET <- as.integer(grepl(antiplatelet_pattern, phenotype$drugs, ignore.case = TRUE))
phenotype$ASPIRIN_OR_ANTIPLATELET <- as.integer(phenotype$ASPIRIN == 1 | phenotype$ANTIPLATELET == 1)
phenotype_all_drugs_filt$STATIN <- phenotype$STATIN
phenotype_all_drugs_filt$ASPIRIN <- phenotype$ASPIRIN
phenotype_all_drugs_filt$ANTIPLATELET <- phenotype$ANTIPLATELET
phenotype_all_drugs_filt$ASPIRIN_OR_ANTIPLATELET <- phenotype$ASPIRIN_OR_ANTIPLATELET
data <- cbind.data.frame(phenotype$ID, phenotype$alcohol, phenotype$PPI, phenotype$STATIN, phenotype$ASPIRIN, 
    phenotype$ANTIPLATELET, phenotype$ASPIRIN_OR_ANTIPLATELET, phenotype$AGE, phenotype$GENDER, phenotype$BMI_mb, 
    phenotype$gramm)
rownames(data) <- phenotype$ID
colnames(data) <- c("ID", "alcohol", "PPI", "STATIN", "ASPIRIN", "ANTIPLATELET", "ASPIRIN_OR_ANTIPLATELET", 
    "Age", "Gender", "BMI", "Food")
data$PPI <- as.factor(data$PPI)
data$PPI <- relevel(data$PPI, ref = "0")
data$STATIN <- as.factor(data$STATIN)
data$STATIN <- relevel(data$STATIN, ref = "0")
data$ASPIRIN <- as.factor(data$ASPIRIN)
data$ASPIRIN <- relevel(data$ASPIRIN, ref = "0")
data$ANTIPLATELET <- as.factor(data$ANTIPLATELET)
data$ANTIPLATELET <- relevel(data$ANTIPLATELET, ref = "0")
data$ASPIRIN_OR_ANTIPLATELET <- as.factor(data$ASPIRIN_OR_ANTIPLATELET)
data$ASPIRIN_OR_ANTIPLATELET <- relevel(data$ASPIRIN_OR_ANTIPLATELET, ref = "0")
data$PPI <- as.factor(data$PPI)
data$PPI <- relevel(data$PPI, ref = "0")
data$Gender <- as.factor(data$Gender)
dat <- phyloseq(otu_table(counts_filt, taxa_are_rows = TRUE), sample_data(data))
focus_data <- cbind(data, phenotype_all_drugs_filt[, 4:16], phenoData_filt[, c(2, 8, 9, 18, 68)])
focus_data$glucose_mg.dl_mb <- as.numeric(focus_data$glucose_mg.dl_mb)
focus_data$diabetes_mb <- as.factor(focus_data$diabetes_mb)
focus_data$BMI.g25.30 <- as.factor(focus_data$BMI.g25.30)
focus_data$Analgesic <- as.factor(focus_data$Analgesic)
focus_data$Analgesic_NSAID <- as.factor(focus_data$Analgesic_NSAID)
focus_data$Antibiotics <- as.factor(focus_data$Antibiotics)
focus_data$Antihypertensives <- as.factor(focus_data$Antihypertensives)
