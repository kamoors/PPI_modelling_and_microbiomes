# SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0
# Oxygen-tolerance reaction classification
# Run through run.R so input and output paths resolve within the package.
# Inputs: see data/input_manifest.csv and docs/RUN_ORDER.md.
# Optional upstream reconstruction: raw reads or metabolic models are external.
# The main package workflows start from the supplied processed tables.

library(data.table)
library(ggplot2)
normalizeOxygenRequirement <- function(x) {
    x <- trimws(x)
    x <- gsub("\\s+", " ", x)
    x
}
mapOxygenClass <- function(x) {
    x <- normalizeOxygenRequirement(x)
    out <- rep(NA_character_, length(x))
    out[x %in% c("Aerobe")] <- "aerobe"
    out[x %in% c("Facultative anaerobe", "Facultative anaerobic", "Facultative", "Aerotolerant")] <- "facultative"
    out[x %in% c("Obligate anaerobe", "Obligate  anaerobe")] <- "anaerobe"
    out
}
loadModelsByMicrobeID <- function(microbeIDs, basePath = ppi_external("agora2_models")) {
    models <- list()
    for (id in microbeIDs) {
        path <- file.path(basePath, paste0(id, ".RDS"))
        if (ppi_file_exists(path)) {
            mod <- tryCatch(readRDS(ppi_input(path)), error = function(e) {
                message("Could not load model: ", id, " (", e$message, ")")
                NULL
            })
            if (!is.null(mod)) 
                models[[id]] <- mod
        }
    }
    models
}
getSpeciesReactionFrequency <- function(microbeIDs, basePath = ppi_external("agora2_models")) {
    models <- loadModelsByMicrobeID(microbeIDs = microbeIDs, basePath = basePath)
    if (length(models) == 0) 
        return(NULL)
    allReaks <- unlist(lapply(models, function(mod) mod@react_id), use.names = FALSE)
    freq <- table(allReaks)/length(models)
    freqNum <- as.numeric(freq)
    names(freqNum) <- names(freq)
    freqNum
}
getSpeciesReactionMatrix <- function(strainData, speciesNames, basePath = ppi_external("agora2_models")) {
    speciesFreqList <- list()
    for (sp in speciesNames) {
        microbeIDs <- strainData[Species == sp, MicrobeID]
        freqSp <- getSpeciesReactionFrequency(microbeIDs = microbeIDs, basePath = basePath)
        if (!is.null(freqSp)) 
            speciesFreqList[[sp]] <- freqSp
    }
    if (length(speciesFreqList) == 0) 
        return(NULL)
    allReaks <- sort(unique(unlist(lapply(speciesFreqList, names), use.names = FALSE)))
    res <- data.table(species = names(speciesFreqList))
    for (reak in allReaks) {
        res[[reak]] <- vapply(speciesFreqList, function(x) {
            if (reak %in% names(x)) 
                x[[reak]]
            else 0
        }, numeric(1))
    }
    res
}
getSpeciesOxygenScores <- function(strainData, minAnnotatedStrains = 1) {
    tmp <- copy(strainData)
    tmp[, `:=`(oxygen_norm, normalizeOxygenRequirement(get("Oxygen Requirement")))]
    tmp[, `:=`(oxygen_class, mapOxygenClass(oxygen_norm))]
    tmp <- tmp[!is.na(oxygen_class)]
    if (nrow(tmp) == 0) 
        return(NULL)
    speciesOxy <- tmp[, .(n_total = .N, n_aerobe = sum(oxygen_class == "aerobe"), n_facultative = sum(oxygen_class == 
        "facultative"), n_anaerobe = sum(oxygen_class == "anaerobe")), by = Species]
    speciesOxy[, `:=`(score_aerobe, n_aerobe/n_total)]
    speciesOxy[, `:=`(score_facultative, n_facultative/n_total)]
    speciesOxy[, `:=`(score_anaerobe, n_anaerobe/n_total)]
    speciesOxy <- speciesOxy[n_total >= minAnnotatedStrains]
    speciesOxy
}
associateReactionsWithOxygenClasses <- function(speciesReactionMat, speciesOxygen, method = "spearman", 
    minSpeciesWithReactionVariation = 3) {
    dt <- merge(speciesReactionMat, speciesOxygen[, .(Species, n_total, score_aerobe, score_facultative, 
        score_anaerobe)], by.x = "species", by.y = "Species")
    reactionCols <- setdiff(colnames(dt), c("species", "n_total", "score_aerobe", "score_facultative", 
        "score_anaerobe"))
    resList <- lapply(reactionCols, function(r) {
        x <- dt[[r]]
        xBase <- x[!is.na(x)]
        n_nonzero <- sum(xBase > 0)
        n_unique_x <- length(unique(xBase))
        if (length(xBase) < 3 || n_nonzero < minSpeciesWithReactionVariation || n_unique_x < 2) {
            return(data.table(reak = r, cor_aerobe = NA_real_, pval_aerobe = NA_real_, cor_facultative = NA_real_, 
                pval_facultative = NA_real_, cor_anaerobe = NA_real_, pval_anaerobe = NA_real_, n_species = length(xBase), 
                n_species_present = n_nonzero, meanReactionFreq = mean(xBase, na.rm = TRUE)))
        }
        getOneCor <- function(y) {
            ok <- !is.na(x) & !is.na(y)
            if (sum(ok) < 3) 
                return(c(cor = NA_real_, p = NA_real_))
            if (length(unique(y[ok])) < 2) 
                return(c(cor = NA_real_, p = NA_real_))
            ct <- suppressWarnings(cor.test(x[ok], y[ok], method = method, exact = FALSE))
            c(cor = unname(ct$estimate), p = ct$p.value)
        }
        aer <- getOneCor(dt$score_aerobe)
        fac <- getOneCor(dt$score_facultative)
        ana <- getOneCor(dt$score_anaerobe)
        data.table(reak = r, cor_aerobe = aer["cor"], pval_aerobe = aer["p"], cor_facultative = fac["cor"], 
            pval_facultative = fac["p"], cor_anaerobe = ana["cor"], pval_anaerobe = ana["p"], n_species = length(xBase), 
            n_species_present = n_nonzero, meanReactionFreq = mean(xBase, na.rm = TRUE))
    })
    res <- rbindlist(resList)
    res[, `:=`(padj_aerobe, p.adjust(pval_aerobe, method = "fdr"))]
    res[, `:=`(padj_facultative, p.adjust(pval_facultative, method = "fdr"))]
    res[, `:=`(padj_anaerobe, p.adjust(pval_anaerobe, method = "fdr"))]
    res
}
assignReactionOxygenClass <- function(assoc, corCutoff = 0.1, margin = 0.05) {
    res <- copy(assoc)
    corrNames <- c("cor_aerobe", "cor_facultative", "cor_anaerobe")
    classNames <- c("aerobe", "facultative", "anaerobe")
    res[, `:=`(assigned_class, "mixed")]
    res[, `:=`(best_cor, NA_real_)]
    res[, `:=`(second_cor, NA_real_)]
    res[, `:=`(best_minus_second, NA_real_)]
    for (i in seq_len(nrow(res))) {
        vals <- as.numeric(res[i, ..corrNames])
        if (all(is.na(vals))) 
            next
        vals2 <- vals
        vals2[is.na(vals2)] <- -Inf
        ord <- order(vals2, decreasing = TRUE)
        bestIdx <- ord[1]
        secondIdx <- ord[2]
        bestVal <- vals2[bestIdx]
        secondVal <- vals2[secondIdx]
        if (!is.finite(secondVal)) 
            secondVal <- NA_real_
        res[i, `:=`(best_cor, bestVal)]
        res[i, `:=`(second_cor, secondVal)]
        res[i, `:=`(best_minus_second, ifelse(is.na(secondVal), NA_real_, bestVal - secondVal))]
        if (!is.finite(bestVal)) 
            next
        if (bestVal >= corCutoff) {
            if (is.na(secondVal) || (bestVal - secondVal) >= margin) 
                res[i, `:=`(assigned_class, classNames[bestIdx])]
        }
    }
    res
}
compareReactionClasses <- function(assocPPI, class1, class2, fcCol = "beta_meta") {
    x1 <- assocPPI[assigned_class == class1, get(fcCol)]
    x2 <- assocPPI[assigned_class == class2, get(fcCol)]
    x1 <- x1[!is.na(x1)]
    x2 <- x2[!is.na(x2)]
    if (length(x1) == 0 || length(x2) == 0) 
        return(NULL)
    wt <- wilcox.test(x1, x2)
    list(class1 = class1, class2 = class2, mean1 = mean(x1), mean2 = mean(x2), median1 = median(x1), 
        median2 = median(x2), n1 = length(x1), n2 = length(x2), p.value = wt$p.value, test = wt)
}
computePairwiseTests <- function(assocPPI, classes = c("aerobe", "facultative", "anaerobe", "mixed"), 
    fcCol = "beta_meta") {
    combs <- t(combn(classes, 2))
    res <- lapply(1:nrow(combs), function(i) {
        cl1 <- combs[i, 1]
        cl2 <- combs[i, 2]
        x1 <- assocPPI[assigned_class == cl1, get(fcCol)]
        x2 <- assocPPI[assigned_class == cl2, get(fcCol)]
        x1 <- x1[!is.na(x1)]
        x2 <- x2[!is.na(x2)]
        if (length(x1) == 0 || length(x2) == 0) 
            return(NULL)
        wt <- wilcox.test(x1, x2)
        data.table(class1 = cl1, class2 = cl2, p = wt$p.value, n1 = length(x1), n2 = length(x2))
    })
    res <- rbindlist(res)
    res[, `:=`(padj, p.adjust(p, method = "fdr"))]
    res
}
plotReactionClasses_withStats <- function(assocPPI, includeClasses = c("aerobe", "facultative", "anaerobe", 
    "mixed"), fcCol = "beta_meta", title = "Reaction classes by oxygen association") {
    dt <- copy(assocPPI)
    dt <- dt[assigned_class %in% includeClasses]
    dt <- dt[!is.na(get(fcCol))]
    if (nrow(dt) == 0) {
        message("No reactions available for plotting.")
        return(NULL)
    }
    counts <- sapply(includeClasses, function(cl) sum(dt$assigned_class == cl))
    dt[, `:=`(class_label, factor(assigned_class, levels = includeClasses, labels = paste0(includeClasses, 
        "\n(n=", counts, ")")))]
    stats <- computePairwiseTests(dt, classes = includeClasses, fcCol = fcCol)
    y_max <- max(dt[[fcCol]], na.rm = TRUE)
    y_min <- min(dt[[fcCol]], na.rm = TRUE)
    y_range <- y_max - y_min
    stats[, `:=`(x1, match(class1, includeClasses))]
    stats[, `:=`(x2, match(class2, includeClasses))]
    stats <- stats[order(padj)]
    stats[, `:=`(y, y_max + (1:.N) * 0.05 * y_range)]
    stats[, `:=`(label, paste0("FDR=", signif(padj, 2)))]
    p <- ggplot(dt, aes(x = class_label, y = get(fcCol))) + geom_boxplot(outlier.shape = NA, width = 0.6) + 
        geom_jitter(width = 0.15, alpha = 0.35) + labs(x = NULL, y = "Reaction PPI association") + theme_bw()
    for (i in 1:nrow(stats)) {
        p <- p + annotate("segment", x = stats$x1[i], xend = stats$x2[i], y = stats$y[i], yend = stats$y[i]) + 
            annotate("text", x = (stats$x1[i] + stats$x2[i])/2, y = stats$y[i], label = stats$label[i], 
                vjust = -0.3, size = 3)
    }
    invisible(list(plot = p, stats = stats))
}
determineOxygenClasses <- function(basePath = ppi_external("agora2_models"), mappingFile = "agora2_mapping.csv", 
    minAnnotatedStrains = 1, minSpeciesWithReactionVariation = 3, correlationMethod = "spearman", corCutoff = 0.1, 
    margin = 0.05, reactionAssignmentFile = "reaction_assignment.csv") {
    strainData <- fread(ppi_input(mappingFile))
    speciesOxygen <- getSpeciesOxygenScores(strainData = strainData, minAnnotatedStrains = minAnnotatedStrains)
    if (is.null(speciesOxygen) || nrow(speciesOxygen) == 0) 
        stop("No species with sufficient oxygen annotation found.")
    speciesReactionMat <- getSpeciesReactionMatrix(strainData = strainData, speciesNames = speciesOxygen$Species, 
        basePath = basePath)
    if (is.null(speciesReactionMat)) 
        stop("No species models could be loaded.")
    speciesOxygen <- speciesOxygen[Species %in% speciesReactionMat$species]
    assoc <- associateReactionsWithOxygenClasses(speciesReactionMat = speciesReactionMat, speciesOxygen = speciesOxygen, 
        method = correlationMethod, minSpeciesWithReactionVariation = minSpeciesWithReactionVariation)
    assocAssigned <- assignReactionOxygenClass(assoc = assoc, corCutoff = corCutoff, margin = margin)
    write.csv(assocAssigned, ppi_output(reactionAssignmentFile))
}
compareReactionFrequencies_threeOxygenClasses <- function(fcFile = ppi_path("human_cohorts/results/Focus_metacardis_react_ab_meta_analysis_focusNewPPI.csv"), 
    filterSignificantFC = FALSE, pCutoff = 0.1, fcCol = "beta_meta", fdrCol = "fdr_meta", reactionCol = "rxn_name", 
    reactionAssignmentFile = "reaction_assignment.csv") {
    assocAssigned <- fread(ppi_input(reactionAssignmentFile))
    fcs <- fread(ppi_input(fcFile))
    assocPPI <- merge(assocAssigned, fcs, by.x = "reak", by.y = reactionCol)
    if (filterSignificantFC) 
        assocPPI <- assocPPI[get(fdrCol) < pCutoff]
    cat("\n")
    cat("Number of reactions tested:", nrow(assocAssigned), "\n")
    cat("Number of reactions with fold-changes:", nrow(assocPPI), "\n")
    cat("\nReaction-class assignment counts:\n")
    print(table(assocAssigned$assigned_class, useNA = "ifany"))
    write.csv(assocAssigned, ppi_output("reaction_assignment.csv", name = "oxygen_reaction_classes.csv"))
    classSummary <- assocPPI[, .(n = .N, mean_fc = mean(get(fcCol), na.rm = TRUE), median_fc = median(get(fcCol), 
        na.rm = TRUE)), by = assigned_class][order(match(assigned_class, c("aerobe", "facultative", "anaerobe", 
        "mixed")))]
    cat("\nFold-change summary by class:\n")
    print(classSummary)
    pairTests <- list(aerobe_vs_facultative = compareReactionClasses(assocPPI, "aerobe", "facultative", 
        fcCol = fcCol), aerobe_vs_anaerobe = compareReactionClasses(assocPPI, "aerobe", "anaerobe", fcCol = fcCol), 
        facultative_vs_anaerobe = compareReactionClasses(assocPPI, "facultative", "anaerobe", fcCol = fcCol), 
        anaerobe_vs_mixed = compareReactionClasses(assocPPI, "anaerobe", "mixed", fcCol = fcCol), facultative_vs_mixed = compareReactionClasses(assocPPI, 
            "facultative", "mixed", fcCol = fcCol), aerobe_vs_mixed = compareReactionClasses(assocPPI, 
            "aerobe", "mixed", fcCol = fcCol))
    cat("\nPairwise Wilcoxon tests:\n")
    for (nm in names(pairTests)) {
        pt <- pairTests[[nm]]
        if (!is.null(pt)) {
            cat("\n", nm, "\n", sep = "")
            cat("means:", pt$mean1, "vs", pt$mean2, "\n")
            cat("medians:", pt$median1, "vs", pt$median2, "\n")
            cat("n:", pt$n1, "vs", pt$n2, "\n")
            cat("p:", pt$p.value, "\n")
        }
    }
    plotReactionClasses_withStats(assocPPI = assocPPI, includeClasses = c("aerobe", "facultative", "anaerobe", 
        "mixed"), fcCol = fcCol, title = ifelse(filterSignificantFC, "Reaction classes by oxygen association (significant PPI reactions)", 
        "Reaction classes by oxygen association"))
}
reactionAssignmentFile <- "reaction_assignment.csv"
determineOxygenClasses(basePath = ppi_external("agora2_models"), mappingFile = "data/agora2_mapping.csv", 
    minAnnotatedStrains = 2, minSpeciesWithReactionVariation = 3, correlationMethod = "spearman", corCutoff = 0.2, 
    margin = 0.2, reactionAssignmentFile)
p1 <- compareReactionFrequencies_threeOxygenClasses(fcFile = "data/Focus_metacardis_react_ab_meta_analysis_focusNewPPI.csv", 
    filterSignificantFC = FALSE, pCutoff = 0.1, fcCol = "beta_meta", fdrCol = "fdr_meta", reactionCol = "rxn_name", 
    reactionAssignmentFile)
saveRDS(p1$plot, ppi_output("anaerobe_plot.RDS", name = "oxygen_reaction_effect_plot.rds"))
