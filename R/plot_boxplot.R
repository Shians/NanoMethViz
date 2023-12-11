#' Plot boxplots for regions
#'
#' @param x the NanoMethResult object.
#' @param regions a table of regions containing at least columns chr, strand,
#'   start and end. Any additional columns can be used for grouping.
#' @param binary_threshold the modification probability such that calls with
#'   modification probability above the threshold are considered methylated, and
#'   those with probability equal or below are considered unmethylated.
#' @param group_col the column to group aggregated trends by. This column can
#'   be in from the regions table or samples(x).
#' @param palette the ggplot colour palette used for groups.
#'
#' @return a ggplot object containing the methylation boxplots.
#'
#' @examples
#' nmr <- load_example_nanomethresult()
#' gene_anno <- exons_to_genes(NanoMethViz::exons(nmr))
#' plot_boxplot(nmr, gene_anno)
#' plot_boxplot(nmr, gene_anno, group_col = "sample")
#' plot_boxplot(nmr, gene_anno, group_col = "group")
#'
#' @export
plot_violin <- function(
    nmr,
    gene_anno,
    binary_threshold = 0.5,
    group_col = NULL,
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

    query_row_methy <- function(i, methy, regions, flank) {
        # query a larger region such that smoothing doesn't
        # misbehave around ends
        flank <- flank * 1.1
        query_methy(
            methy,
            regions$chr[i],
            regions$start[i] - flank,
            regions$end[i] + flank,
            force = TRUE)
    }

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
        dplyr::filter(purrr::map_lgl(methy_data, function(x) nrow(x) != 0))

    region_data <- regions %>%
        dplyr::select(!dplyr::any_of(c("chr", "strand", "start", "end"))) %>%
        tidyr::unnest("methy_data") %>%
        dplyr::summarise(methy_prop = mean(.data$mod_prob > binary_threshold), .by = c("gene_id", "symbol", "sample", "chr", "strand")) %>%
        dplyr::inner_join(samples(x), by = "sample", multiple = "all")

    if (!is.null(group_col)) {
        aes_spec <- ggplot2::aes(
                x = .data[[group_col]],
                y = .data$methy_prop,
                col = .data$group)
    } else {
        aes_spec <- ggplot2::aes(
                x = .data$binned_pos,
                y = .data$methy_prop)
    }

    ggplot2::ggplot(region_data, aes_spec) +
        ggplot2::geom_violin() +
        ggplot2::scale_y_continuous(limits = c(0, 1)) +
        ggplot2::theme_bw() +
        palette
}
