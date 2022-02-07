#' Calculate region methylation statistics
#'
#' Calculate the average methylation probability and prevalence based on
#' specified probability threshold.
#'
#' @param nmr the NanoMethResult object
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
    methy_stats <- purrr::map_df(
        1:nrow(regions),
        get_region_methy_stats,
        nmr = nmr,
        regions = regions,
        threshold = threshold
    )

    dplyr::bind_cols(regions, methy_stats)
}

get_region_methy_stats <- function(index, nmr, regions, threshold) {
    methy <- query_methy(
        nmr,
        regions$chr[index],
        regions$start[index],
        regions$end[index],
        force = TRUE
    )

    if (nrow(methy) != 0) {
        methy %>%
            summarise(
                mean_methy_prob = mean(e1071::sigmoid(.data$statistic)),
                prevalence = mean(e1071::sigmoid(.data$statistic) > threshold)
            )
    } else {
        tibble(mean_methy_prob = NA_real_, prevalence = NA_real_)
    }
}
