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
#' @importFrom scico scale_colour_scico
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
    pos_style,
    xlim = NA
) {
    assertthat::assert_that(
        nrow(exons(x)) > 0,
        msg = "exons(x) is empty, gene cannot be queried"
    )

    if (!anyNA(xlim)) {
        assertthat::assert_that(
            is.numeric(xlim),
            length(xlim) == 2,
            xlim[1] < xlim[2]
        )
    }

    if (length(window_prop) == 1) {
        # convert to two sided window
        window_prop <- c(window_prop, window_prop)
    }

    methy_data <- query_methy_gene(x, gene, window_prop = window_prop)

    # add sample information
    methy_data <- dplyr::left_join(
        methy_data,
        NanoMethViz::samples(x),
        by = "sample"
    )

    read_data <- methy_data %>%
        dplyr::group_by(.data$read_name) %>%
        dplyr::summarise(start = min(.data$pos), end = max(.data$pos))

    group_data <- methy_data %>%
            dplyr::group_by(.data$read_name) %>%
            dplyr::summarise(group = unique(.data$group))

    # get read starts and ends
    read_data <- dplyr::left_join(read_data, group_data, by = "read_name")

    # get grouping indices to pack reads
    append_read_group <- function(x, k) {
        x$read_group <- paste0(k, stacked_interval_inds(x))
        x
    }
    grouping_data <- read_data %>%
        dplyr::group_by(.data$group) %>%
        dplyr::group_modify(append_read_group) %>%
        dplyr::ungroup()

    methy_data$mod_prob <- e1071::sigmoid(methy_data$statistic)

    methy_data <- dplyr::left_join(
        methy_data,
        dplyr::select(grouping_data, "read_name", "read_group"),
        by = "read_name"
    )

    # heatmap theme
    theme_methy_heatmap <- function() {
        theme_minimal() +
            ggplot2::theme(
                axis.ticks.y = ggplot2::element_blank(),
                axis.title.y = ggplot2::element_blank(),
                axis.text.y = ggplot2::element_blank()
            )
    }

    if (pos_style == "compact") {
        # only plots sites with measured modification, evenly spaced
        p <- ggplot(methy_data,
            aes(x = factor(.data$pos),
                y = .data$read_group,
                fill = .data$mod_prob)) +
            scico::scale_colour_scico(palette = 'imola', direction = -1) +
            ggplot2::geom_raster() +
            ggplot2::facet_wrap(~group, scales = "free_y", nrow = 2) +
            theme_methy_heatmap() +
            ggplot2::theme(
                axis.ticks.x = ggplot2::element_blank(),
                axis.text.x = ggplot2::element_blank(),
            ) +
            ggplot2::xlab("Site")
    } else if (pos_style == "to_scale") {
        # plots all sites in range, evenly spaced with square geoms
        # data will overlap
        p <- ggplot(methy_data, aes(y = .data$read_group)) +
            ggplot2::geom_errorbarh(
                ggplot2::aes(
                    xmin = .data$start,
                    xmax = .data$end
                ),
                data = dplyr::left_join(
                    read_data, grouping_data,
                    by = c("read_name", "start", "end", "group"))
            ) +
            ggplot2::geom_point(
                aes(x = .data$pos, col = .data$mod_prob), alpha = 0.33, shape = 15) +
            scico::scale_colour_scico(palette = 'imola', direction = -1) +
            ggplot2::facet_wrap(~group, scales = "free_y", nrow = 2, strip.position = "left") +
            theme_methy_heatmap() +
            ggplot2::xlab("Position")
    }

    if (!anyNA(xlim)) {
        p <- p + ggplot2::xlim(xlim[1], xlim[2])
    }

    p
}
