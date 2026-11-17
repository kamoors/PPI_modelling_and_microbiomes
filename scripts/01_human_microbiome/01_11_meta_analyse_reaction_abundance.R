# SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0
# Fixed-effect reaction-abundance meta-analysis from cohort-level HC3 tables.
# Uses HC3 standard errors from each cohort.
# Reported non-users with metabolomics-detected PPI were excluded in the
# MetaCardis source fits; no participant data are needed for this step.

source("scripts/helpers/setup.R")
invisible(run_analysis(
  section = "01_human_microbiome",
  directory = "human_cohorts/results",
  script = "01_11_meta_analyse_reaction_abundance.R",
  restricted = FALSE,
  code = {
kiel <- read.csv(project_input(project_path('human_cohorts/results/Focus_mgpipe_reaction_abundance_LOG_HC3_FULL_withNewPPI_results.csv')),check.names=FALSE)
metacardis <- read.csv(project_input(project_path('human_cohorts/results/MetaCardis_mgpipe_reaction_abundance_cleanPPI_LOG_HC3_FULL_clean_results.csv')),check.names=FALSE)
names(kiel)[names(kiel)=='metabolite'] <- 'rxn_name'
names(metacardis)[names(metacardis)=='rxn'] <- 'rxn_name'
stopifnot(!anyDuplicated(kiel$rxn_name),!anyDuplicated(metacardis$rxn_name))
result <- merge(kiel,metacardis,by='rxn_name',suffixes=c('.f','.m'),sort=TRUE)
annotations <- read.csv(project_input(project_path('annotations/reaction_meta_annotations.csv')),check.names=FALSE)
result <- merge(result,annotations,by='rxn_name',all.x=TRUE,sort=TRUE)
result$w_focus <- 1/result$hc3_se.f^2
result$w_meta <- 1/result$hc3_se.m^2
result$beta_meta <- (result$w_focus*result$hc3_coef.f+result$w_meta*result$hc3_coef.m)/(result$w_focus+result$w_meta)
result$se_meta <- sqrt(1/(result$w_focus+result$w_meta))
result$z_meta <- result$beta_meta/result$se_meta
result$p_meta <- 2*pnorm(-abs(result$z_meta))
result$fdr_meta <- p.adjust(result$p_meta,method='fdr')
result$same_direction <- sign(result$hc3_coef.f)==sign(result$hc3_coef.m)
result$sig_class <- ifelse(result$hc3_fdr.f<0.05 & result$hc3_fdr.m<0.05,'Significant in both',
 ifelse(result$hc3_fdr.f<0.05,'Kiel only',ifelse(result$hc3_fdr.m<0.05,'MetaCardis only','Not significant')))
write.csv(result,project_output(project_path('human_cohorts/results/Focus_metacardis_react_ab_meta_analysis_focusNewPPI.csv'),name='reaction_abundance_meta_analysis.csv'),row.names=FALSE)
# Register the source-data alias so subsequent Figure 1 runs consume this result.
write.csv(result,project_output(project_path('figures/source_data/Cohort_Microbiome_Function/reaction_abundance_meta_analysis.csv'),name='reaction_abundance_meta_analysis_figure_input.csv'),row.names=FALSE)
}))
