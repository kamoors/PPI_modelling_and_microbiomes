# SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0
# Caecum/faeces LC-MS import function

process_metabolite_sheet <- function(data, metid, meta_cols = 9) {
    require(dplyr)
    metadata <- data %>% dplyr::select(1:meta_cols)
    feat_data <- data %>% dplyr::select((meta_cols + 1):ncol(data))
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
        bad <- is.na(new_names) | trimws(new_names) == ""
        if (any(bad)) {
            print("BAD NEW NAMES")
            print(tibble::tibble(keep_id = keep_ids[bad], new_name = new_names[bad], id_match = match(keep_ids[bad],
                metid_sub$id), rank_primary = rank_primary, primary_value = metid_sub[[rank_primary]][match(keep_ids[bad],
                metid_sub$id)], rank2_value = metid_sub$Metabolite_Rank2[match(keep_ids[bad], metid_sub$id)],
                rank2_1_value = if ("Metabolite_Rank2_1" %in% names(metid_sub)) {
                  metid_sub$Metabolite_Rank2_1[match(keep_ids[bad], metid_sub$id)]
                }
                else {
                  NA
                }, rank3_value = if ("Metabolite_Rank3" %in% names(metid_sub)) {
                  metid_sub$Metabolite_Rank3[match(keep_ids[bad], metid_sub$id)]
                }
                else {
                  NA
                }, rank = metid_sub$Rank[match(keep_ids[bad], metid_sub$id)]))
        }
        colnames(df) <- new_names
        df
    }
    high_met <- metid %>% filter(Rank != "Rank3")
    high_conf_df <- subset_and_rename(high_met, "Metabolite_Rank1")
    low_met <- metid %>% filter(Rank == "Rank3")
    low_conf_df <- subset_and_rename(low_met, "Metabolite_Rank3")
    list(metadata = metadata, high_conf = high_conf_df, low_conf = low_conf_df)
}
