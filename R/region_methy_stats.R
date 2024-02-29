#' Calculate region methylation statistics
#'
#' Calculate the average methylation probability and prevalence based on
#' specified probability threshold.
#'
#' @param nmr the NanoMethResult object.
#' @param regions the table of regions to query statistics for.
#' @param threshold the threshold to use for determining methylation calls
#'   for the calculation of prevalence.
#'
#' @return table of regions with additional columns of methylation summary
#'   statistics.
#' @export
#'
#' @examples
#' nmr <- load_example_nanomethresult()
#' gene_anno <- exons_to_genes(NanoMethViz::exons(nmr))
#' region_methy_stats(nmr, gene_anno)
region_methy_stats <- function(nmr, regions, threshold = 0.5) {

    # calculate methylation statistics for each region using
    methy_stats <- purrr::map_df(
        seq_len(nrow(regions)),
        get_region_methy_stats,
        nmr = nmr,
        regions = regions,
        threshold = threshold
    )

    dplyr::bind_cols(regions, methy_stats)
}

get_region_methy_stats <- function(index, nmr, regions, threshold) {

    # query methylation data for the region
    chr <- regions$chr[index]
    start <- regions$start[index]
    end <- regions$end[index]
    methy <- query_methy(nmr, chr, start, end, force = TRUE)

    # calculate methylation statistics if methylation data is available
    if (nrow(methy) != 0) {
        methy %>%
            summarise(
                mean_methy_prob = mean(sigmoid(.data$statistic)),
                prevalence = mean(sigmoid(.data$statistic) > threshold)
            )
    } else {
        tibble(mean_methy_prob = NA_real_, prevalence = NA_real_)
    }
}
