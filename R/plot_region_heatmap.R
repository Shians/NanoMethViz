#' Plot region methylation heatmap
#'
#' @param x the NanoMethResult object.
#' @param chr the chromosome to plot.
#' @param start the start of the plotting region.
#' @param end the end of the plotting region.
#' @param ... additional arguments.
#'
#' @return a ggplot object of the heatmap.
#'
#' @examples
#' nmr <- load_example_nanomethresult()
#' plot_region_heatmap(nmr, "chr7", 6703892, 6730431)
#'
#' @export
setGeneric("plot_region_heatmap", function(x, chr, start, end, ...) {
    standardGeneric("plot_region_heatmap")
})

#' @rdname plot_region_heatmap
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
#' plot_region_heatmap(nmr, "chr7", 6703892, 6730431)
#'
#' @export
setMethod(
    "plot_region_heatmap",
    signature(
        x = "NanoMethResult",
        chr = "character",
        start = "numeric",
        end = "numeric"),

    function(
        x,
        chr,
        start,
        end,
        pos_style = c("to_scale", "compact"),
        window_prop = 0.3
    ) {
        pos_style <- match.arg(pos_style)

        .plot_region_heatmap(
            x = x,
            chr = chr,
            start = start,
            end = end,
            pos_style = pos_style,
            window_prop = window_prop
        )
    }
)

#' @rdname plot_region_heatmap
#'
#' @export
setMethod("plot_region_heatmap",
    signature(
        x = "NanoMethResult",
        chr = "factor",
        start = "numeric",
        end = "numeric"),

    function(
        x,
        chr,
        start,
        end,
        pos_style = c("to_scale", "compact"),
        window_prop = 0.3
    ) {
        chr <- as.character(chr)
        plot_region(
            x = x,
            chr = chr,
            start = start,
            end = end,
            pos_style = pos_style,
            window_prop = window_prop
        )
    }
)

#' @importFrom scico scale_colour_scico
.plot_region_heatmap <- function(
    x,
    chr,
    start,
    end,
    window_prop,
    pos_style,
    xlim = NA
) {
    if (!anyNA(xlim)) {
        assertthat::assert_that(
            is.numeric(xlim),
            length(xlim) == 2,
            xlim[1] < xlim[2]
        )
    }

    if (length(window_prop) == 1) {
        window_prop <- c(window_prop, window_prop)
    }

    window_left <- (end - start) * window_prop[1]
    window_right <- (end - start) * window_prop[2]

    methy_data <- query_methy(x, chr, start - window_left, end + window_right)

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
        theme_bw() +
            ggplot2::theme(
                axis.ticks.y = ggplot2::element_blank(),
                axis.title.y = ggplot2::element_blank(),
                axis.text.y = ggplot2::element_blank()
            )
    }

    if (pos_style == "compact") {
        # only plots sites with measured modification, evenly spaced
        p <- ggplot2::ggplot(methy_data, aes(x = factor(.data$pos), y = .data$read_group, fill = .data$mod_prob)) +
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
        p <- ggplot2::ggplot(methy_data, aes(y = .data$read_group)) +
            ggplot2::geom_errorbarh(
                ggplot2::aes(
                    xmin = .data$start,
                    xmax = .data$end
                ),
                data = dplyr::left_join(
                    read_data,
                    grouping_data,
                    by = c("read_name", "start", "end", "group"))
            ) +
            ggplot2::geom_point(aes(x = .data$pos, col = .data$mod_prob), alpha = 0.33, shape = 15) +
            scico::scale_colour_scico(palette = 'imola', direction = -1) +
            ggplot2::facet_wrap(~group, scales = "free_y", nrow = 2) +
            theme_methy_heatmap() +
            ggplot2::xlab("Position")
    }

    if (!anyNA(xlim)) {
        p <- p + ggplot2::xlim(xlim[1], xlim[2])
    }

    p
}
