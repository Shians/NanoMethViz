#' Plot gene methylation heatmap
#'
#' @param x the NanoMethResult object.
#' @param gene the gene symbol for the gene to plot.
#' @param ... additional arguments
#'
#' @return a ggplot object of the heatmap
#'
#' @examples
#' nmr <- load_example_nanomethresult()
#' plot_gene_heatmap(nmr, "Peg3")
#'
#' @export
setGeneric("plot_gene_heatmap", function(x, gene, ...) {
    standardGeneric("plot_gene_heatmap")
})

#' @rdname plot_gene_heatmap
#'
#' @param window_prop the size of flanking region to plot. Can be a vector of two
#'   values for left and right window size. Values indicate proportion of gene
#'   length.
#' @param pos_style the style for plotting the base positions along the x-axis.
#'   Defaults to "to_scale", plotting (potentially) overlapping squares
#'   along the genomic position to scale. The "compact" options plots only the
#'   positions with measured modification.
#'
#' @return a ggplot plot containing the heatmap.
#'
#' @examples
#' nmr <- load_example_nanomethresult()
#' plot_gene_heatmap(nmr, "Peg3")
#'
#' @export
setMethod(
    "plot_gene_heatmap",
    signature(x = "NanoMethResult", gene = "character"),
    function(
        x,
        gene,
        window_prop = 0.3,
        pos_style = c("to_scale", "compact")
    ) {
        pos_style <- match.arg(pos_style)

        .plot_gene_heatmap(
            x = x,
            gene = gene,
            window_prop = window_prop,
            pos_style = pos_style
        )
    }
)

.plot_gene_heatmap <- function(
    x,
    gene,
    window_prop,
    pos_style
) {
    assertthat::assert_that(
        nrow(exons(x)) > 0,
        msg = "exons(x) is empty, gene cannot be queried"
    )

    if (length(window_prop) == 1) {
        # convert to two sided window
        window_prop <- c(window_prop, window_prop)
    }

    methy_data <- query_methy_gene(
        nmeth_results, gene, window_prop = window_prop)

    # add sample information
    methy_data <- dplyr::left_join(
        methy_data,
        NanoMethViz::samples(nmeth_results),
        by = "sample"
    )

    read_data <- methy_data %>%
        dplyr::group_by(read_name) %>%
        dplyr::summarise(start = min(pos), end = max(pos))

    group_data <- methy_data %>%
            dplyr::group_by(read_name) %>%
            dplyr::summarise(group = unique(group))

    # get read starts and ends
    read_data <- dplyr::left_join(read_data, group_data, by = "read_name")

    # get grouping indices to pack reads
    append_read_group <- function(x, k) {
        x$read_group <- paste0(k, stacked_interval_inds(x))
        x
    }
    grouping_data <- read_data %>%
        group_by(group) %>%
        group_modify(append_read_group) %>%
        ungroup()

    methy_data$mod_prob <- e1071::sigmoid(methy_data$statistic)

    methy_data <- left_join(
        methy_data,
        dplyr::select(grouping_data, read_name, read_group),
        by = "read_name"
    )

    # heatmap theme
    theme_methy_heatmap <- function() {
        ggthemes::theme_tufte() +
            ggplot2::theme(
                axis.ticks.y = element_blank(),
                axis.title.y = element_blank(),
                axis.text.y = element_blank()
            )
    }

    if (pos_style == "compact") {
        # only plots sites with measured modification, evenly spaced
        ggplot(methy_data, aes(x = factor(pos), y = read_group, fill = mod_prob)) +
            scale_fill_gradient(low = "blue", high = "red") +
            geom_raster() +
            facet_wrap(~group, scales = "free_y", nrow = 2) +
            theme_methy_heatmap() +
            ggplot2::theme(
                axis.ticks.x = element_blank(),
                axis.text.x = element_blank(),
            ) +
            xlab("Site")
    } else if (pos_style == "to_scale") {
        # plots all sites in range, evenly spaced with square geoms
        # data will overlap
        ggplot(methy_data, aes(y = read_group)) +
            scale_colour_gradient(low = "blue", high = "red") +
            ggplot2::geom_errorbarh(
                ggplot2::aes(
                    xmin = .data$start,
                    xmax = .data$end
                ),
                data = dplyr::left_join(read_data, grouping_data)
            ) +
            geom_point(aes(x = pos, col = mod_prob), alpha = 0.33, shape = 15) +
            facet_wrap(~group, scales = "free_y", nrow = 2) +
            theme_methy_heatmap() +
            xlab("Position")
    }
}
