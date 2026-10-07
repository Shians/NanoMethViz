#' Plot violin for regions
#'
#' This function plots a violin plot of the methylation proportion for each
#' region in the regions table. The methylation proportion is calculated as the
#' mean of the modification probability within each region, and the violin shows
#' the distribution across groups. Regions are grouped and coloured by the
#' `group_col` column in the `regions` table or `samples(x)`.
#'
#' @param x the NanoMethResult object.
#' @param regions a table of regions containing at least columns chr, strand,
#'   start and end. Any additional columns can be used for grouping.
#' @param binary_threshold the modification probability such that calls with
#'   modification probability above the threshold are considered methylated, and
#'   those with probability equal or below are considered unmethylated.
#' @param group_col the column to group and colour violins by. This column can
#'   be from the regions table or samples(x). If NULL, a single violin is drawn
#'   for all regions and samples.
#' @param palette the ggplot colour palette used for groups.
#'
#' @return a ggplot object containing the methylation violin plot.
#'
#' @examples
#' nmr <- load_example_nanomethresult()
#' gene_anno <- exons_to_genes(NanoMethViz::exons(nmr))
#' plot_violin(nmr, gene_anno)
#' plot_violin(nmr, gene_anno, group_col = "sample")
#'
#' @export
plot_violin <- function(
    x,
    regions,
    binary_threshold = 0.5,
    group_col = "group",
    palette = ggplot2::scale_colour_brewer(palette = "Set1")
) {
    if (!is.null(group_col)) {
        avail_columns <- c(colnames(samples(x)), colnames(regions))
        assertthat::assert_that(
            group_col %in% avail_columns,
            msg = glue::glue("'{group_col}' could not be found in columns of 'regions' or samples(x)")
        )
    }

    # grouped regions crashes downstream operations
    regions <- ungroup(regions)

    regions$methy_data <- purrr::map(
        seq_len(nrow(regions)),
        ~query_methy(
            x,
            regions$chr[.x],
            regions$start[.x],
            regions$end[.x],
            force = TRUE
        )
    )

    # remove regions with no data
    regions <- regions %>%
        dplyr::filter(purrr::map_lgl(.data$methy_data, function(x) nrow(x) != 0))

    # summarise each region per sample, keeping any annotation columns of the
    # regions table so they remain available for grouping
    region_cols <- setdiff(
        colnames(regions),
        c("chr", "strand", "start", "end", "methy_data")
    )
    region_data <- regions %>%
        dplyr::mutate(.region_id = dplyr::row_number()) %>%
        dplyr::select(!dplyr::any_of(c("chr", "strand", "start", "end"))) %>%
        tidyr::unnest("methy_data") %>%
        dplyr::summarise(
            methy_prop = mean(.data$mod_prob > binary_threshold),
            .by = dplyr::all_of(c(".region_id", region_cols, "sample"))) %>%
        dplyr::inner_join(samples(x), by = "sample", multiple = "all")

    if (!is.null(group_col)) {
        aes_spec <- ggplot2::aes(
                x = .data[[group_col]],
                y = .data$methy_prop,
                col = .data[[group_col]])
    } else {
        # a single violin over all regions and samples
        aes_spec <- ggplot2::aes(
                x = "",
                y = .data$methy_prop)
    }

    # draw_quantiles was deprecated in ggplot2 4.0.0
    if (utils::packageVersion("ggplot2") >= "4.0.0") {
        violin <- ggplot2::geom_violin(quantiles = 0.5, quantile.linetype = 1)
    } else {
        violin <- ggplot2::geom_violin(draw_quantiles = 0.5)
    }

    ggplot2::ggplot(region_data, aes_spec) +
        violin +
        ggplot2::scale_y_continuous(limits = c(0, 1)) +
        ggplot2::theme_bw() +
        palette
}
