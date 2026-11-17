# SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0
# Serum import function

process_metabolite_sheet <- function(data, metid, meta_cols = 9) {
    require(dplyr)
    metadata <- data %>% dplyr::select(1:meta_cols)
    feat_data <- data %>% dplyr::select((meta_cols + 1):ncol(data))
    all_features <- feat_data
    subset_and_rename <- function(metid_sub, rank_primary) {
        keep_ids <- intersect(colnames(feat_data), metid_sub$id)
        df <- feat_data %>% dplyr::select(all_of(keep_ids))
        new_names <- vapply(keep_ids, function(orig) {
            i <- match(orig, metid_sub$id)
            nm <- metid_sub[[rank_primary]][i]
            if (is.na(nm) || nm == "") {
                nm <- metid_sub$Metabolite_Rank2[i]
            }
            if ((is.na(nm) || nm == "") && "Metabolite_Rank2_1" %in% names(metid_sub)) {
                nm <- metid_sub$Metabolite_Rank2_1[i]
            }
            if (is.na(nm) || nm == "") {
                nm <- orig
            }
            nm
        }, FUN.VALUE = character(1), USE.NAMES = FALSE)
        colnames(df) <- new_names
        df
    }
    high_met <- metid %>% filter(Rank != "Rank3")
    high_conf_df <- subset_and_rename(high_met, "Metabolite_Rank1")
    low_met <- metid %>% filter(Rank == "Rank3")
    low_conf_df <- subset_and_rename(low_met, "Metabolite_Rank3")
    list(metadata = metadata, all_features = all_features, high_conf = high_conf_df, low_conf = low_conf_df)
}
