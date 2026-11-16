# SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0
# Regression check: the two Figure 3 distance panels must use different contexts.
file_argument<-grep("^--file=",commandArgs(FALSE),value=TRUE)
root<-dirname(dirname(normalizePath(sub("^--file=","",file_argument),winslash="/")))
required<-c("ggplot2","dplyr")
stopifnot(all(vapply(required,requireNamespace,logical(1),quietly=TRUE)))
expressions<-parse(file.path(root,"scripts/03_mouse_microbiome/03_02_plot_mouse_microbiome.R"),encoding="UTF-8")
hit<-Filter(function(x)is.call(x)&&identical(x[[1]],as.name("<-"))&&
    identical(x[[2]],as.name("make_distance_panel")),as.list(expressions))
stopifnot(length(hit)==1L)
e<-new.env(parent=.GlobalEnv)
e$results_df<-read.csv(file.path(root,"data/figure_inputs/microbiome16s/Figure3_BrayCurtis_plot_data.csv"))
e$base_theme<-ggplot2::theme_bw()
e$group_cols<-c(nonrecipient="white",recipient="grey70")
eval(hit[[1]],e)
culture<-e$make_distance_panel("Ex vivo")$data
mouse<-e$make_distance_panel("Mouse")$data
stopifnot(identical(unique(culture$sample_type),"Ex vivo"),identical(unique(mouse$sample_type),"Mouse"),
    nrow(culture)==144L,nrow(mouse)==96L,!identical(culture,mouse))
p<-function(df,donor) {
    df<-df[df$donor_group==donor,]
    wilcox.test(df$distance[df$group=="recipient"],df$distance[df$group=="nonrecipient"])$p.value
}
stopifnot(abs(p(culture,"Ctrl_1")-2/choose(24,4))<1e-12,
          abs(p(mouse,"Ctrl_1")-2/choose(16,2))<1e-12)
cat("Figure 3 contexts stay separate; recipient counts and checked donor-comparison P values match.\n")
