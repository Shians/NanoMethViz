#' Plot region
#'
#' @param x the NanoMethResult object.
#' @param chr the chromosome to plot.
#' @param start the start of the plotting region.
#' @param end the end of the plotting region.
#' @param ... additional arguments.
#'
#' @return a patchwork plot containing the methylation profile in the specified
#'   region.
#'
#' @examples
#' nmr <- load_example_nanomethresult()
#' plot_region(nmr, "chr7", 6703892, 6730431)
#'
#' @importFrom ggrastr rasterise
#' @export
setGeneric("plot_region", function(x, chr, start, end, ...) {
    standardGeneric("plot_region")
})

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
#' @param span the span for loess smoothing.
#' @param window_prop the size of flanking region to plot. Can be a vector of two
#'   values for left and right window size. Values indicate proportion of gene
#'   length.
#' @param palette the ggplot colour palette used for groups.
#' @param line_size the size of the lines.
#'
#' @return a patchwork plot containing the methylation profile in the specified
#'   region.
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
        span = NULL,
        window_prop = 0,
        palette = ggplot2::scale_colour_brewer(palette = "Set1"),
        line_size = 2
    ) {
        avg_method = match.arg(avg_method)

        .plot_region(
            x = x,
            chr = chr,
            start = start,
            end = end,
            anno_regions = anno_regions,
            binary_threshold = binary_threshold,
            spaghetti = spaghetti,
            avg_method = avg_method,
            heatmap = heatmap,
            span = span,
            window_prop = window_prop,
            palette = palette,
            line_size = line_size
        )
    }
)

#' @rdname plot_region
#'
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
        span = NULL,
        window_prop = 0,
        palette = ggplot2::scale_colour_brewer(palette = "Set1"),
        line_size = 2
    ) {
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
            span = span,
            window_prop = window_prop,
            palette = palette,
            line_size = line_size
        )
    }
)

.plot_region <- function(
    x,
    chr,
    start,
    end,
    anno_regions,
    binary_threshold,
    avg_method,
    spaghetti,
    heatmap,
    span,
    window_prop,
    palette,
    line_size
) {
    sample_anno <- samples(x)
    exons_anno <- query_exons_region(
        exons(x),
        chr = chr,
        start = start,
        end = end
    )

    if (length(window_prop) == 1) {
        window_prop <- c(window_prop, window_prop)
    }

    feature_width <- end - start
    window_left <- feature_width * window_prop[1]
    window_right <- feature_width * window_prop[2]
    xlim <- round(c(start - window_left, end + window_right))

    methy_data <-
        query_methy(
            x,
            chr,
            floor(start - window_left * 1.1),
            ceiling(end + window_right * 1.1),
            simplify = TRUE) %>%
        dplyr::select(-"strand") %>%
        tibble::as_tibble()

    if (nrow(methy_data) == 0) {
        warning("no methylation data in region")
        return(ggplot() + theme_void())
    }

    title <- glue::glue("{chr}:{start}-{end}")
    p1 <- plot_methylation_internal(
        methy_data = methy_data,
        start = start,
        end = end,
        chr = chr,
        xlim = xlim,
        title = title,
        anno_regions = anno_regions,
        binary_threshold = binary_threshold,
        avg_method = avg_method,
        spaghetti = spaghetti,
        sample_anno = sample_anno,
        span = span,
        palette_col = palette,
        line_size = line_size
    ) +
        ggplot2::scale_x_continuous(
            limits = xlim,
            expand = ggplot2::expansion(),
            labels = scales::label_number(scale_cut = scales::cut_si("b"))
        )

    p2 <- plot_gene_annotation(exons_anno, xlim[1], xlim[2]) +
        ggplot2::scale_x_continuous(
            limits = xlim,
            expand = ggplot2::expansion(),
            labels = scales::label_number(scale_cut = scales::cut_si("b"))
        )

    anno_height <- attr(p2, "plot_height")

    heights <- c(1, 0.075 * anno_height)
    p_out <- p1 / p2 + patchwork::plot_layout(heights = heights)

    if (heatmap) {
        p_heatmap <- plot_region_heatmap(x, chr, start, end, window_prop = window_prop) +
            ggplot2::scale_x_continuous(
                limits = xlim,
                expand = ggplot2::expansion(),
                labels = scales::label_number(scale_cut = scales::cut_si("b"))
            )

        p_out <- stack_plots(p_out, ggrastr::rasterise(p_heatmap, dpi = 300))
    }

    p_out
}
