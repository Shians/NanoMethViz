#' @rdname plot_region
#'
#' @param anno_regions the data.frame of regions to be annotated.
#' @param binary_threshold the modification probability such that calls with
#'   modification probability above the threshold are set to 1 and probabilities
#'   equal to or below the threshold are set to 0.
#' @param avg_method the average method for pre-smoothing at each genomic position.
#'   Data is pre-smoothed at each genomic position before the smoothed aggregate line
#'   is generated for performance reasons. The default is "mean" which corresponds
#'   to the average methylation fraction. The alternative "median" option is
#'   closer to an average within the more common methylation state.
#' @param spaghetti whether or not individual reads should be shown.
#' @param heatmap whether or not read-methylation heatmap should be shown.
#' @param heatmap_subsample how many packed rows of reads to subsample to.
#' @param smoothing_window the window size for smoothing the trend line.
#' @param window_prop the size of flanking region to plot. Can be a vector of two
#'   values for left and right window size. Values indicate proportion of gene
#'   length.
#' @param palette the ggplot colour palette used for groups.
#' @param line_size the size of the lines.
#' @param mod_scale the scale range for modification probabilities. Default c(0, 1), set to "auto" for automatic
#'   limits.
#' @param span DEPRECATED, use smoothing_window instead. Will be removed in next version.
#'
#' @details
#' This function plots the methylation data for a given region. The region is specified by
#' chromosome, start and end positions. The basic plot contains a smoothed line plot of
#' the methylation of each experimental group. Since V3.0.0 NanoMethViz has changed the
#' smoothing strategy from a loess smoothing to a weighted moving average. This is because
#' the loess smoothing was too computationally expensive for large datasets and had a
#' span parameter that was difficult to tune. The new smoothing strategy is controlled
#' by the smoothing_window argument.
#'
#' @examples
#' nmr <- load_example_nanomethresult()
#' plot_region(nmr, "chr7", 6703892, 6730431)
#'
#' @export
setMethod("plot_region",
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
        anno_regions = NULL,
        binary_threshold = NULL,
        avg_method = c("mean", "median"),
        spaghetti = FALSE,
        heatmap = FALSE,
        heatmap_subsample = 50,
        smoothing_window = 2000,
        window_prop = 0,
        palette = ggplot2::scale_colour_brewer(palette = "Set1"),
        line_size = 1,
        mod_scale = c(0, 1),
        span = NULL
    ) {
        if (!missing("span")) {
            warning("the 'span' argument has been deprecated, please use 'smoothing_window' instead")
        }
        avg_method <- match.arg(avg_method)

        plot_region_impl(
            x = x,
            chr = chr,
            start = start,
            end = end,
            anno_regions = anno_regions,
            binary_threshold = binary_threshold,
            spaghetti = spaghetti,
            avg_method = avg_method,
            heatmap = heatmap,
            heatmap_subsample = heatmap_subsample,
            smoothing_window = smoothing_window,
            window_prop = window_prop,
            palette = palette,
            line_size = line_size,
            mod_scale = mod_scale
        )
    }
)

#' @rdname plot_region
#' @export
setMethod("plot_region",
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
        anno_regions = NULL,
        binary_threshold = NULL,
        avg_method = c("mean", "median"),
        spaghetti = FALSE,
        heatmap = FALSE,
        heatmap_subsample = 50,
        smoothing_window = 2000,
        window_prop = 0,
        palette = ggplot2::scale_colour_brewer(palette = "Set1"),
        line_size = 1,
        mod_scale = c(0, 1),
        span = NULL
    ) {
        if (!missing("span")) {
            warning("the 'span' argument has been deprecated, please use 'smoothing_window' instead")
        }
        avg_method <- match.arg(avg_method)

        plot_region_impl(
            x = x,
            chr = chr,
            start = start,
            end = end,
            anno_regions = anno_regions,
            binary_threshold = binary_threshold,
            spaghetti = spaghetti,
            avg_method = avg_method,
            heatmap = heatmap,
            heatmap_subsample = heatmap_subsample,
            smoothing_window = smoothing_window,
            window_prop = window_prop,
            palette = palette,
            line_size = line_size,
            mod_scale = mod_scale
        )
    }
)

#' @rdname plot_region
#' @export
setMethod("plot_region",
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
        anno_regions = NULL,
        binary_threshold = NULL,
        avg_method = c("mean", "median"),
        spaghetti = FALSE,
        heatmap = FALSE,
        heatmap_subsample = 50,
        smoothing_window = 2000,
        window_prop = 0,
        palette = ggplot2::scale_colour_brewer(palette = "Set1"),
        line_size = 1,
        mod_scale = c(0, 1),
        span = NULL
    ) {
        if (!missing("span")) {
            warning("the 'span' argument has been deprecated, please use 'smoothing_window' instead")
        }
        chr <- as.character(chr)
        avg_method <- match.arg(avg_method)
        plot_region(
            x = x,
            chr = chr,
            start = start,
            end = end,
            anno_regions = anno_regions,
            binary_threshold = binary_threshold,
            avg_method = avg_method,
            spaghetti = spaghetti,
            heatmap = heatmap,
            heatmap_subsample = heatmap_subsample,
            smoothing_window = smoothing_window,
            window_prop = window_prop,
            palette = palette,
            line_size = line_size,
            mod_scale = mod_scale
        )
    }
)

#' @rdname plot_region
#' @export
setMethod("plot_region",
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
        anno_regions = NULL,
        binary_threshold = NULL,
        avg_method = c("mean", "median"),
        spaghetti = FALSE,
        heatmap = FALSE,
        heatmap_subsample = 50,
        smoothing_window = 2000,
        window_prop = 0,
        palette = ggplot2::scale_colour_brewer(palette = "Set1"),
        line_size = 1,
        mod_scale = c(0, 1),
        span = NULL
    ) {
        if (!missing("span")) {
            warning("the 'span' argument has been deprecated, please use 'smoothing_window' instead")
        }
        chr <- as.character(chr)
        avg_method <- match.arg(avg_method)
        plot_region(
            x = x,
            chr = chr,
            start = start,
            end = end,
            anno_regions = anno_regions,
            binary_threshold = binary_threshold,
            avg_method = avg_method,
            spaghetti = spaghetti,
            heatmap = heatmap,
            heatmap_subsample = heatmap_subsample,
            smoothing_window = smoothing_window,
            window_prop = window_prop,
            palette = palette,
            line_size = line_size,
            mod_scale = mod_scale
        )
    }
)

plot_region_impl <- function(
    x,
    chr,
    start,
    end,
    anno_regions,
    binary_threshold,
    avg_method,
    spaghetti,
    heatmap,
    heatmap_subsample,
    smoothing_window,
    window_prop,
    palette,
    line_size,
    mod_scale
) {
    sample_anno <- samples(x)

    if (length(window_prop) == 1) {
        window_prop <- c(window_prop, window_prop)
    }

    feature_width <- end - start
    window_left <- feature_width * window_prop[1]
    window_right <- feature_width * window_prop[2]

    # query data
    methy_data <- query_methy(
        x,
        chr,
        floor(start - window_left * 1.1),
        ceiling(end + window_right * 1.1),
        simplify = TRUE
    )

    if (nrow(methy_data) == 0) {
        warning("no methylation data in region, returning empty plot")
        return(ggplot() + theme_void())
    }

    methy_data <- methy_data %>%
        dplyr::select(-"strand") %>%
        tibble::as_tibble()


    # setup base plot
    title <- glue::glue("{chr}:{start}-{end}")
    xlim <- round(c(start - window_left, end + window_right))
    p1 <- plot_methylation_data(
        methy_data = methy_data,
        start = start,
        end = end,
        chr = chr,
        title = title,
        anno_regions = anno_regions,
        binary_threshold = binary_threshold,
        avg_method = avg_method,
        spaghetti = spaghetti,
        sample_anno = sample_anno,
        smoothing_window = smoothing_window,
        palette_col = palette,
        line_size = line_size,
        mod_scale = mod_scale
    ) +
        ggplot2::coord_cartesian(
            xlim = xlim,
            expand = FALSE
        )

    p_out <- p1

    # if exon anno exists, append it to plot
    if (nrow(exons(x)) != 0) {
        exons_anno <- query_exons_region(exons(x), chr = chr, start = start, end = end)

        p2 <- plot_gene_annotation(exons_anno, xlim[1], xlim[2]) +
            ggplot2::coord_cartesian(xlim = xlim, expand = FALSE) +
            ggplot2::scale_x_continuous(labels = scales::label_number(scale_cut = scales::cut_si("b")))

        anno_height <- attr(p2, "plot_height")
        heights <- c(1, 0.075 * anno_height)

        p_out <- p1 / p2 + patchwork::plot_layout(heights = heights)
    }

    # if heatmap requested, append it to plot
    if (heatmap) {
        p_heatmap <- plot_region_heatmap(x, chr, start, end, window_prop = window_prop, subsample = heatmap_subsample) +
            ggplot2::coord_cartesian(
                xlim = xlim
            )

        p_out <- stack_plots(p_out, ggrastr::rasterise(p_heatmap, dpi = 300))
    }

    p_out
}
