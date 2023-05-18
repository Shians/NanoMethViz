#' @rdname plot_region_heatmap
#'
#' @param window_prop the size of flanking region to plot. Can be a vector of two
#'   values for left and right window size. Values indicate proportion of gene
#'   length.
#' @param pos_style the style for plotting the base positions along the x-axis.
#'   Defaults to "to_scale", plotting (potentially) overlapping squares
#'   along the genomic position to scale. The "compact" options plots only the
#'   positions with measured modification.
#' @param subsample the number of read of packed read rows to subsample to.
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
        window_prop = 0,
        subsample = 50

    ) {
        pos_style <- match.arg(pos_style)

        .plot_region_heatmap(
            x = x,
            chr = chr,
            start = start,
            end = end,
            pos_style = pos_style,
            window_prop = window_prop,
            subsample = subsample
        )
    }
)

#' @rdname plot_region_heatmap
#'
#' @export
setMethod(
    "plot_region_heatmap",
    signature(
        x = "ModBamResult",
        chr = "character",
        start = "numeric",
        end = "numeric"),

    function(
        x,
        chr,
        start,
        end,
        pos_style = c("to_scale", "compact"),
        window_prop = 0,
        subsample = 50

    ) {
        pos_style <- match.arg(pos_style)

        .plot_region_heatmap(
            x = x,
            chr = chr,
            start = start,
            end = end,
            pos_style = pos_style,
            window_prop = window_prop,
            subsample = subsample
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
        window_prop = 0,
        subsample = 50
    ) {
        chr <- as.character(chr)
        .plot_region_heatmap(
            x = x,
            chr = chr,
            start = start,
            end = end,
            pos_style = pos_style,
            window_prop = window_prop,
            subsample = subsample
        )
    }
)

#' @rdname plot_region_heatmap
#'
#' @export
setMethod("plot_region_heatmap",
    signature(
        x = "ModBamResult",
        chr = "factor",
        start = "numeric",
        end = "numeric"),

    function(
        x,
        chr,
        start,
        end,
        pos_style = c("to_scale", "compact"),
        window_prop = 0,
        subsample = 50
    ) {
        chr <- as.character(chr)
        .plot_region_heatmap(
            x = x,
            chr = chr,
            start = start,
            end = end,
            pos_style = pos_style,
            window_prop = window_prop,
            subsample = subsample
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
    subsample
) {
    if (length(window_prop) == 1) {
        window_prop <- c(window_prop, window_prop)
    }

    window_left <- (end - start) * window_prop[1]
    window_right <- (end - start) * window_prop[2]

    plot_left <- round(start - window_left)
    plot_right <- round(end + window_right)

    methy_data <- query_methy(x, chr, plot_left, plot_right)

    # add sample information
    methy_data <- dplyr::left_join(
        NanoMethViz::samples(x),
        methy_data,
        by = "sample",
        multiple = "all"
    )

    read_data <- methy_data %>%
        dplyr::group_by(.data$read_name) %>%
        dplyr::summarise(start = min(.data$pos), end = max(.data$pos))

    group_data <- methy_data %>%
            dplyr::group_by(.data$read_name) %>%
            dplyr::summarise(group = unique(.data$group))

    # get read starts and ends
    read_data <- dplyr::left_join(read_data, group_data, by = "read_name", multiple = "all")

    # get grouping indices to pack reads
    append_read_group <- function(x, k) {
        x$read_group <- paste0(k, stacked_interval_inds(x))
        x
    }
    grouping_data <- read_data %>%
        dplyr::group_by(.data$group) %>%
        dplyr::group_modify(append_read_group) %>%
        dplyr::ungroup()

    methy_data$mod_prob <- sigmoid(methy_data$statistic)

    methy_data <- dplyr::left_join(
        methy_data,
        dplyr::select(grouping_data, "read_name", "read_group"),
        by = "read_name",
        multiple = "all"
    )

    plot_heatmap_internal(
        methy_data = methy_data,
        read_data = read_data,
        grouping_data = grouping_data,
        pos_style = pos_style,
        subsample = subsample
    )
}
