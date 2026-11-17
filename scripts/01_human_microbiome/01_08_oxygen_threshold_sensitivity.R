# SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0
# Oxygen-classification threshold sensitivity

source("scripts/helpers/setup.R")
invisible(run_analysis(
  section = "01_human_microbiome",
  directory = "human_cohorts/results/figs_PPI_Paper/figs/robustness",
  script = "01_08_oxygen_threshold_sensitivity.R",
  restricted = FALSE,
  code = {
project_setwd(project_path("human_cohorts/results/figs_PPI_Paper/figs/robustness"))
library(data.table)
library(ggplot2)
runThresholdSensitivity <- function(reactionAssignmentFile = project_path("human_cohorts/results/figs_PPI_Paper/figs/reaction_assignment.csv"),
    fcFile = project_path("human_cohorts/results/Focus_metacardis_react_ab_meta_analysis_focusNewPPI.csv"),
    corCutoffs = c(0.1, 0.15, 0.2, 0.25, 0.3), margins = c(0.1, 0.15, 0.2, 0.25, 0.3), referenceCorCutoff = 0.2,
    referenceMargin = 0.2) {
    assoc <- fread(project_input(reactionAssignmentFile), encoding = "Latin-1", check.names = TRUE)[, .(reak,
        cor_aerobe, cor_facultative, cor_anaerobe)]
    effects <- fread(project_input(fcFile), encoding = "Latin-1")[, .(rxn_name, beta_meta, fdr_meta)]
    reference <- assignReactionOxygenClass(assoc, corCutoff = referenceCorCutoff, margin = referenceMargin)[,
        .(reak, reference_class = assigned_class)]
    thresholdGrid <- CJ(cor_cutoff = corCutoffs, margin = margins)
    classResults <- vector("list", nrow(thresholdGrid))
    testResults <- vector("list", nrow(thresholdGrid))
    stabilityResults <- vector("list", nrow(thresholdGrid))
    for (i in seq_len(nrow(thresholdGrid))) {
        currentCor <- thresholdGrid$cor_cutoff[i]
        currentMargin <- thresholdGrid$margin[i]
        assigned <- assignReactionOxygenClass(assoc, corCutoff = currentCor, margin = currentMargin)
        dat <- merge(assigned, effects, by.x = "reak", by.y = "rxn_name", all = FALSE)
        dat <- merge(dat, reference, by = "reak", all.x = TRUE)
        classResults[[i]] <- dat[!is.na(beta_meta), .(n = .N, mean_beta = mean(beta_meta), median_beta = median(beta_meta),
            q25_beta = quantile(beta_meta, 0.25), q75_beta = quantile(beta_meta, 0.75), significant_positive = sum(fdr_meta <
                0.1 & beta_meta > 0, na.rm = TRUE), significant_negative = sum(fdr_meta < 0.1 & beta_meta <
                0, na.rm = TRUE)), by = assigned_class][, `:=`(cor_cutoff = currentCor, margin = currentMargin)]
        currentTests <- computePairwiseTests(dat, classes = c("aerobe", "facultative", "anaerobe", "mixed"),
            fcCol = "beta_meta")
        setnames(currentTests, "padj", "p_bh")
        currentTests[, `:=`(cor_cutoff = currentCor, margin = currentMargin, comparison = paste(class1,
            class2, sep = "_vs_"))]
        testResults[[i]] <- currentTests
        referenceNonMixed <- dat$reference_class != "mixed"
        stabilityResults[[i]] <- data.table(cor_cutoff = currentCor, margin = currentMargin, overall_agreement = mean(dat$assigned_class ==
            dat$reference_class), reference_nonmixed_retained = mean(dat$assigned_class[referenceNonMixed] ==
            dat$reference_class[referenceNonMixed]), n_reference_nonmixed = sum(referenceNonMixed))
    }
    list(class_summary = rbindlist(classResults), pairwise_tests = rbindlist(testResults), stability = rbindlist(stabilityResults))
}
library(data.table)
library(ggplot2)
assignReactionOxygenClass <- function(assoc, corCutoff = 0.1, margin = 0.05) {
    res <- copy(as.data.table(assoc))
    correlationColumns <- c("cor_aerobe", "cor_facultative", "cor_anaerobe")
    classNames <- c("aerobe", "facultative", "anaerobe")
    res[, `:=`(assigned_class = "mixed", best_cor = NA_real_, second_cor = NA_real_, best_minus_second = NA_real_)]
    for (i in seq_len(nrow(res))) {
        values <- as.numeric(res[i, ..correlationColumns])
        if (all(is.na(values))) {
            next
        }
        values[is.na(values)] <- -Inf
        orderValues <- order(values, decreasing = TRUE)
        bestIndex <- orderValues[1]
        secondIndex <- orderValues[2]
        bestValue <- values[bestIndex]
        secondValue <- values[secondIndex]
        if (!is.finite(secondValue)) {
            secondValue <- NA_real_
        }
        difference <- ifelse(is.na(secondValue), NA_real_, bestValue - secondValue)
        res[i, `:=`(best_cor = bestValue, second_cor = secondValue, best_minus_second = difference)]
        passesThreshold <- is.finite(bestValue) && bestValue >= corCutoff && (is.na(secondValue) || difference >=
            margin)
        if (passesThreshold) {
            res[i, `:=`(assigned_class, classNames[bestIndex])]
        }
    }
    res
}
computePairwiseTests <- function(assocPPI, classes = c("aerobe", "facultative", "anaerobe", "mixed"),
    fcCol = "beta_meta") {
    comparisons <- t(combn(classes, 2))
    results <- lapply(seq_len(nrow(comparisons)), function(i) {
        class1 <- comparisons[i, 1]
        class2 <- comparisons[i, 2]
        values1 <- assocPPI[assigned_class == class1, get(fcCol)]
        values2 <- assocPPI[assigned_class == class2, get(fcCol)]
        values1 <- values1[!is.na(values1)]
        values2 <- values2[!is.na(values2)]
        if (length(values1) == 0 || length(values2) == 0) {
            return(NULL)
        }
        test <- wilcox.test(values1, values2, exact = FALSE)
        data.table(class1 = class1, class2 = class2, p = test$p.value, n1 = length(values1), n2 = length(values2))
    })
    results <- rbindlist(results)
    results[, `:=`(padj, p.adjust(p, method = "BH"))]
    results
}
thresholdSensitivity <- runThresholdSensitivity()
fwrite(thresholdSensitivity$class_summary, project_output("oxygen_threshold_sensitivity_class_summary.csv",
    name = "oxygen_threshold_class_summary.csv"))
fwrite(thresholdSensitivity$pairwise_tests, project_output("oxygen_threshold_sensitivity_pairwise_tests.csv",
    name = "oxygen_threshold_pairwise_tests.csv"))
fwrite(thresholdSensitivity$stability, project_output("oxygen_threshold_sensitivity_stability.csv", name = "oxygen_threshold_classification_stability.csv"))
medianComparison <- dcast(thresholdSensitivity$class_summary, cor_cutoff + margin ~ assigned_class, value.var = "median_beta")
medianComparison[, `:=`(ordering_preserved, aerobe > anaerobe & facultative > anaerobe)]
print(medianComparison)
cat("\nOrdering preserved in all threshold combinations:", all(medianComparison$ordering_preserved),
    "\n")
robustnessSummary <- thresholdSensitivity$pairwise_tests[, .(largest_BH_adjusted_P = max(p_bh, na.rm = TRUE)),
    by = comparison][order(largest_BH_adjusted_P)]
print(robustnessSummary)
thresholdPlot <- ggplot(thresholdSensitivity$class_summary, aes(x = margin, y = median_beta, colour = assigned_class,
    group = assigned_class)) + geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") + geom_line(linewidth = 0.8) +
    geom_point(size = 2) + geom_vline(xintercept = 0.2, linetype = "dotted", colour = "black") + facet_wrap(~cor_cutoff,
    labeller = label_both) + scale_colour_manual(values = c(aerobe = "#0072B2", facultative = "#E69F00",
    anaerobe = "#009E73", mixed = "#7F7F7F")) + labs(x = "Minimum correlation difference", y = "Median reaction meta-analysis effect",
    colour = "Reaction class") + theme_classic() + theme(legend.position = "bottom")
ggsave(project_output("oxygen_threshold_sensitivity.png", name = "oxygen_threshold_sensitivity.png"), thresholdPlot,
    width = 10, height = 6, dpi = 300)
}))
