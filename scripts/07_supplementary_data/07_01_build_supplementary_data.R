# SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0
# Assemble the published Supplementary Data tables. This is a reporting step.
# Statistical refitting is performed by scripts/01 through scripts/06.
if (!requireNamespace("openxlsx", quietly=TRUE)) stop("Install openxlsx first.")
names <- c("Microbiome_and_Cohort_Analyses", "Experimental_Metabolomics",
           "Mouse_Transcriptomics", "CrossLayer_Integration")
out_dir <- ppi_path("supplementary_data")
for (i in seq_along(names)) {
    tables <- readRDS(ppi_input(ppi_path(paste0("reference/supplementary_data_",i,"_tables.rds"))))
    wb <- openxlsx::createWorkbook(creator="Karlis A. Moors")
    header <- openxlsx::createStyle(fontColour="#FFFFFF",fgFill="#17365D",textDecoration="bold",wrapText=TRUE,valign="top")
    for (sheet in names(tables)) {
        openxlsx::addWorksheet(wb,sheet,gridLines=FALSE)
        tab <- tables[[sheet]]
        # The Excel writer inserts CR before LF. Normalize existing CRLF first
        # so multiline labels do not acquire duplicate carriage returns.
        tab[] <- lapply(tab, function(x) if (is.character(x)) gsub("\r\n", "\n", x, fixed=TRUE) else x)
        openxlsx::writeData(wb,sheet,tab,headerStyle=header,withFilter=nrow(tab)>0)
        openxlsx::freezePane(wb,sheet,firstRow=TRUE)
        openxlsx::setColWidths(wb,sheet,cols=seq_len(ncol(tab)),widths=20)
        openxlsx::setRowHeights(wb,sheet,rows=1,heights=30)
    }
    filename <- paste0("supplementary_data_",i,"_",tolower(names[i]),".xlsx")
    openxlsx::saveWorkbook(wb,ppi_output(file.path(out_dir,filename),name=filename),overwrite=TRUE)
}
# Copy the mixed-type NMR workbook without rewriting cell values or metadata.
source <- ppi_input(ppi_path("reference/supplementary_data_5_nmr.xlsx"))
target <- ppi_output(file.path(out_dir,"supplementary_data_5_nmr_source_data.xlsx"),
                     name="supplementary_data_5_nmr_source_data.xlsx")
stopifnot(file.copy(source,target,overwrite=TRUE))
message("Five supplementary workbooks assembled from published table snapshots.")
