#' @rdname plot_gene_heatmap
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
        pos_style = c("to_scale", "compact"),
        subsample = 50
    ) {
        pos_style <- match.arg(pos_style)

        .plot_gene_heatmap(
            x = x,
            gene = gene,
            window_prop = window_prop,
            pos_style = pos_style,
            subsample = subsample
        )
    }
)

#' @rdname plot_gene_heatmap
#'
#' @export
setMethod(
    "plot_gene_heatmap",
    signature(x = "ModBamResult", gene = "character"),
    function(
        x,
        gene,
        window_prop = 0.3,
        pos_style = c("to_scale", "compact"),
        subsample = 50
    ) {
        pos_style <- match.arg(pos_style)

        .plot_gene_heatmap(
            x = x,
            gene = gene,
            window_prop = window_prop,
            pos_style = pos_style,
            subsample = subsample
        )
    }
)

.plot_gene_heatmap <- function(
    x,
    gene,
    window_prop,
    pos_style,
    subsample
) {
    assertthat::assert_that(
        nrow(exons(x)) > 0,
        msg = "exons(x) is empty, gene cannot be queried"
    )

    if (length(window_prop) == 1) {
        # convert to two sided window
        window_prop <- c(window_prop, window_prop)
    }

    # query_gene_methy
    if (!gene %in% exons(x)$symbol) {
        stop(glue::glue("gene {gene} not found in exon annotation"))
    }

    pos_range <- gene_pos_range(x, gene)

    gene_width <- pos_range[2] - pos_range[1]
    window_left <- gene_width * window_prop[1]
    window_right <- gene_width * window_prop[2]

    plot_left <- pos_range[1] - window_left
    plot_right <- pos_range[2] + window_right

    chr <- exons(x) %>%
        dplyr::filter(.data$symbol == gene) %>%
        dplyr::slice(1) %>%
        dplyr::pull(chr)

    methy_data <- query_methy(
        x,
        chr = chr,
        start = plot_left,
        end = plot_right
    )

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


    methy_data$mod_prob <- e1071::sigmoid(methy_data$statistic)

    methy_data <- dplyr::left_join(
        methy_data,
        dplyr::select(grouping_data, "read_name", "read_group"),
        by = "read_name",
        multiple = "all"
    )

    # subsample reads down to a certain number of read groups
    subsample_groups <- function(x, key, subsample) {
        if (nrow(x) > subsample) {
            x <- dplyr::sample_n(x, subsample)
        }
        x
    }
    methy_data <- methy_data %>%
        dplyr::nest_by(.data$group, .data$read_group) %>%
        dplyr::group_by(.data$group) %>%
        dplyr::group_modify(subsample_groups, subsample = subsample)
    methy_data <- tidyr::unnest(methy_data, "data")
    read_data <- read_data %>%
        dplyr::filter(.data$read_name %in% methy_data$read_name)

    if (pos_style == "compact") {
        # only plots sites with measured modification, evenly spaced
        p <- ggplot(methy_data,
            aes(x = factor(.data$pos),
                y = .data$read_group,
                fill = .data$mod_prob)) +
            ggplot2::geom_raster() +
            heatmap_fill_scale +
            ggplot2::facet_wrap(~group, scales = "free_y", ncol = 1, strip.position = "right") +
            theme_methy_heatmap +
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
                    by = c("read_name", "start", "end", "group"),
                    multiple = "all"
                ),
                alpha = 1,
                color = "darkgray",
                linewidth = 1.2,
                height = 0
            ) +
            ggplot2::geom_point(
                aes(x = .data$pos, col = .data$mod_prob), alpha = 1, shape = 15) +
            heatmap_col_scale +
            ggplot2::facet_wrap(~group, scales = "free_y", ncol = 1, strip.position = "right") +
            theme_methy_heatmap +
            ggplot2::xlab("Position")
    }

    p
}
