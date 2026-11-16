# SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0
# Fixed-effect predicted-secretion meta-analysis from cohort-level HC3 tables.
# Use the current Kiel result and the filtered MetaCardis/Results table.
# The identically named standalone table in Results/Kiel_Metacardis is a
# different version and is deliberately not an input here.
kiel <- read.csv(ppi_input(ppi_path('human_cohorts/results/Focus_mgpipe_net_secretion_HC3_FULL_NEWPPI_clean_results.csv')),check.names=FALSE)
metacardis <- read.csv(ppi_input(ppi_path('human_cohorts/metacardis_results/Metacardis_mgsec_fullFormula_totalMedicationAdjusted_noPPIExcludedByMetabolomics_preBoot.csv')),check.names=FALSE)
# Preserve the legacy row index for exact compatibility with the final table.
if(''%in%names(metacardis))names(metacardis)[names(metacardis)==''] <- 'X'
stopifnot(!anyDuplicated(kiel$metabolite),!anyDuplicated(metacardis$rxn))
kiel <- kiel[,c('metabolite','fullName','hc3_coef','hc3_se','hc3_p','hc3_fdr')]
result <- merge(kiel,metacardis,by.x='metabolite',by.y='rxn',suffixes=c('.focus','.metac'),sort=TRUE)
result$w1 <- 1/result$hc3_se.focus^2
result$w2 <- 1/result$hc3_se.metac^2
result$beta_meta <- (result$w1*result$hc3_coef.focus+result$w2*result$hc3_coef.metac)/(result$w1+result$w2)
result$se_meta <- sqrt(1/(result$w1+result$w2))
result$z_meta <- result$beta_meta/result$se_meta
result$p_meta <- 2*pnorm(-abs(result$z_meta))
result$fdr_meta <- p.adjust(result$p_meta,method='fdr')
result$same_direction <- sign(result$hc3_coef.focus)==sign(result$hc3_coef.metac)
result$sig_class <- ifelse(result$hc3_fdr.focus<0.05 & result$hc3_fdr.metac<0.05,'Significant in both',
 ifelse(result$hc3_fdr.focus<0.05,'Kiel only',ifelse(result$hc3_fdr.metac<0.05,'Metacardis only','Not significant')))
write.csv(result,ppi_output(ppi_path('human_cohorts/results/Kiel_Metacardis_mg_sec_comparison_with_meta_focus_newPPI.csv'),name='predicted_secretion_meta_analysis.csv'),row.names=FALSE)
write.csv(result,ppi_output(ppi_path('figures/source_data/Cohort_Microbiome_Function/Kiel_MetaCardis_predicted_secretion_meta.csv'),name='predicted_secretion_meta_analysis_figure_input.csv'),row.names=FALSE)
