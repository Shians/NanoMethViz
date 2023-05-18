#' Cluster reads based on methylation
#'
#' @param x a ModBamResult object.
#' @param chr the chromosome name where to find the region.
#' @param start the start position of the region.
#' @param end the end position of the region.
#' @param min_pts the minimum number of points needed to form a cluster (default = 10).
#'
#' @return A tibble with information about each read's cluster assignment and read statistics.
#'
#' @import tidyr
#' @import dplyr
#' @import dbscan
#' @importFrom tibble rownames_to_column
cluster_reads <- function(x, chr, start, end, min_pts = 10) {
    # query data
    methy_data <- query_methy(x, chr, start, end)

    if (nrow(methy_data) == 0) {
        stop(glue::glue("no reads containing methylation data found in specified region"))
    }

    methy_data <- methy_data %>%
        dplyr::filter(pos >= start & pos < end)

    read_stats <- get_read_stats(methy_data)

    # identify the read names whose span is at least 90% the length of maximum span
    # filter methylation data for only those reads that meet the above condition of span
    max_span <- max(read_stats$span)
    keep_reads <- read_stats$read_name[read_stats$span > 0.9 * max_span]
    methy_data <- methy_data %>%
        dplyr::filter(read_name %in% keep_reads)

    # convert methylation data into a matrix with one row for each read name
    mod_mat <- methy_data %>%
        dplyr::select(read_name, pos, mod_prob) %>%
        dplyr::arrange(pos) %>%
        tidyr::pivot_wider(names_from = pos, values_from = mod_prob) %>%
        df_to_matrix()

    # pre-check before filtering
    if (nrow(mod_mat) < min_pts) {
        stop(glue::glue("fewer reads available ({nrow(mod_mat)} reads) than minimum cluster size 'min_pts' ({min_pts})"))
    }

    # remove positions with high missingness (>60%) then reads with high missingness (>30%)
    mod_mat_filled <- mod_mat[order(rownames(mod_mat)), ]
    col_missingness <- mat_col_map(mod_mat_filled, missingness)
    mod_mat_filled <- mod_mat_filled[, col_missingness < 0.6]
    row_missingness <- mat_row_map(mod_mat_filled, missingness)
    mod_mat_filled <- mod_mat_filled[row_missingness < 0.3, ]

    # fill in missing values with mean methylation probability across that read
    for (i in 1:nrow(mod_mat_filled)) {
        mod_mat_filled[i, is.na(mod_mat_filled[i, ])] <- mean(mod_mat_filled[i, ], na.rm = TRUE)
    }

    # post-check before filtering
    if (nrow(mod_mat_filled) < min_pts) {
        stop(glue::glue("fewer reads available ({nrow(mod_mat_filled)} reads) than minimum cluster size 'min_pts' ({min_pts})"))
    }

    # cluster reads using HDBSCAN algorithm with specified minimum number of points
    dbsc <- dbscan::hdbscan(mod_mat_filled, minPts = min_pts)
    clust_df <- data.frame(read_name = rownames(mod_mat_filled), cluster_id = dbsc$cluster)

    # merge and process results of cluster analysis and read statistics
    clust_df %>%
        dplyr::left_join(read_stats, by = "read_name") %>%
        dplyr::arrange(cluster_id) %>%
        dplyr::mutate(
            cluster_id = as.factor(cluster_id),
            start = as.integer(start),
            end = as.integer(end),
            span = as.integer(span)
        )
}

# summarize read statistics (start, end, strand) based on same read name
get_read_stats <- function(methy_data) {
    methy_data %>%
        group_by(read_name) %>%
        summarise(
            start = min(pos),
            end = max(pos),
            mean = mean(mod_prob, na.rm = TRUE),
            span = end - start,
            strand = unique(strand)
        ) %>%
        arrange(strand)
}

