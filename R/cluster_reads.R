cluster_reads <- function(x, chr, start, end) {
    methy_data <- query_methy(x, chr, start, end ) %>%
        filter(pos >= start & pos < end)

    read_stats <- methy_data %>%
        group_by(read_name) %>%
        summarise(
            start = min(pos),
            end = max(pos),
            span = end - start,
            strand = unique(strand)
        ) %>%
        arrange(strand)

    max_span <- max(read_stats$span)

    keep_reads <- read_stats$read_name[read_stats$span > 0.9 * max_span]

    methy_data <- methy_data %>%
        dplyr::filter(read_name %in% keep_reads)

    mod_mat <- methy_data %>%
        dplyr::select(read_name, pos, mod_prob) %>%
        dplyr::arrange(pos) %>%
        tidyr::pivot_wider(names_from = pos, values_from = mod_prob) %>%
        df_to_matrix()

    mod_mat_filled <- mod_mat[order(rownames(mod_mat)), ]

    col_missingness <- mat_col_map(mod_mat_filled, missingness)
    mod_mat_filled <- mod_mat_filled[, col_missingness < 0.8]

    row_missingness <- mat_row_map(mod_mat_filled, missingness)
    mod_mat_filled <- mod_mat_filled[row_missingness < 0.2, ]
    for (i in 1:nrow(mod_mat_filled)) {
        mod_mat_filled[i, is.na(mod_mat_filled[i, ])] <- mean(mod_mat_filled[i, ], na.rm = TRUE)
    }

    dbsc <- dbscan::hdbscan(mod_mat_filled, minPts = 15)
    clust_df <- data.frame(read_name = rownames(mod_mat_filled), cluster_id = dbsc$cluster)

    out_df <- as.data.frame(mod_mat) %>%
        tibble::rownames_to_column("read_name") %>%
        tidyr::pivot_longer(
            cols = !dplyr::contains("read_name"),
            names_to = "pos",
            values_to = "methy_prob"
        )

    out_df %>%
        dplyr::group_by(read_name) %>%
        dplyr::summarise(mean = mean(methy_prob, na.rm = TRUE)) %>%
        dplyr::left_join(clust_df, by = "read_name") %>%
        dplyr::left_join(read_stats, by = "read_name") %>%
        dplyr::arrange(cluster_id)
}
