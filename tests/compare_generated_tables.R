# SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0
# Compare generated CSV result tables with matching saved tables in data/.
# Run analyses first, then: Rscript tests/compare_generated_tables.R results
file_argument <- grep("^--file=",commandArgs(FALSE),value=TRUE)
code_root <- dirname(dirname(normalizePath(sub("^--file=","",file_argument),winslash="/")))
if (.Platform$OS.type=="windows") invisible(Sys.setlocale("LC_CTYPE","English_United States.utf8"))
args <- commandArgs(TRUE)
results <- if(length(args)) normalizePath(args[1],winslash="/",mustWork=TRUE) else file.path(code_root,"results")
registry_file <- file.path(results,"output_manifest.csv")
if(!file.exists(registry_file))stop("Run an analysis before comparing results.")
registry <- read.csv(registry_file,stringsAsFactors=FALSE)
manifest <- read.csv(file.path(code_root,"data/input_manifest.csv"),stringsAsFactors=FALSE)
manifest <- manifest[manifest$kind %in% c("analysis_selection_result","aggregate_human_result","aggregate_model_statistics",
    "aggregate_figure_input","aggregate_experimental_result"),]
read_table <- function(path) read.csv(path,check.names=FALSE,stringsAsFactors=FALSE)
strip_index <- function(x) {
 if(ncol(x)&&names(x)[1] %in% c("","X","...1")&&identical(as.character(x[[1]]),as.character(seq_len(nrow(x)))))x<-x[,-1,drop=FALSE]
 x
}
report <- list()
for(i in seq_len(nrow(registry))) {
 logical <- registry$logical_path[i]
 if(!grepl("\\.csv(\\.gz)?$",logical,ignore.case=TRUE))next
 if(grepl("source_manifest|resolved_input_paths",basename(logical)))next
 key <- sub("/NMR_annotated_limma_dupCor/","/NMR_exvivo_annotated_limma_dupCor/",logical,fixed=TRUE)
 hits <- which(tolower(manifest$logical_path)==tolower(key))
 if(!length(hits))hits<-which(tolower(basename(manifest$logical_path))==tolower(basename(key)))
 paths <- unique(manifest$package_path[hits])
 if(length(paths)!=1L)next
 reference <- file.path(code_root,"data",paths)
 generated <- file.path(results,registry$output_path[i])
 if(!file.exists(generated))next
 a<-strip_index(read_table(reference));b<-strip_index(read_table(generated))
 original_generated_rows<-nrow(b)
 keys <- intersect(c("metabolite","gene","reaction","rxn_name","rxn","analysis_group","feature"),names(a))
 keys <- keys[keys %in% names(b)]
 aligned<-FALSE
 same_keys<-TRUE
 for(key in keys)if(!anyNA(a[[key]])&&!anyNA(b[[key]])&&!anyDuplicated(a[[key]])&&!anyDuplicated(b[[key]])) {
   same_keys<-setequal(a[[key]],b[[key]])
   b<-b[match(a[[key]],b[[key]]),,drop=FALSE];aligned<-TRUE;break
 }
 shape<-nrow(a)==original_generated_rows&&same_keys&&identical(dim(a),dim(b))&&setequal(names(a),names(b))
 numeric <- names(a)[vapply(a,is.numeric,logical(1))]
 failures<-0L
 for(field in numeric) {
  x<-a[[field]];y<-b[[field]]
  if(is.null(y)||length(x)!=length(y)||!is.numeric(y)){failures<-failures+max(1L,length(x));next}
  good<-(is.na(x)&is.na(y))|(is.finite(x)&is.finite(y)&abs(x-y)<=1e-12+1e-8*abs(x))|
       (!is.na(x)&!is.na(y)&!is.finite(x)&!is.finite(y)&x==y)
  good[is.na(good)]<-FALSE;failures<-failures+sum(!good)
 }
 report[[length(report)+1L]]<-data.frame(table=basename(logical),reference=paths,
    reference_rows=nrow(a),generated_rows=original_generated_rows,numeric_columns=length(numeric),
    numeric_mismatches=failures,status=if(shape&&failures==0L)"NUMERIC_MATCH" else "DIFFERS")
}
if(!length(report))stop("No generated CSV tables have an unambiguous bundled reference; this does not validate the run.")
report<-do.call(rbind,report)
write.csv(report,file.path(results,"reference_comparison.csv"),row.names=FALSE)
print(report[,c("table","numeric_columns","numeric_mismatches","status")],row.names=FALSE)
cat("Comparison covers matching saved CSV result tables only, not every output or database-dependent enrichment.\n")
if(any(report$status=="DIFFERS"))quit(status=1L)
