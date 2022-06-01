#' Plot GRanges
#'
#' @param x the NanoMethResult object.
#' @param grange the GRanges object with one entry.
#' @inheritParams plot_region
#' @inherit plot_region return
#'
#' @examples
#' nmr <- load_example_nanomethresult()
#' plot_grange(nmr, GRanges("chr7:6703892-6730431"))
#'
#' @export
plot_grange <- function(
    x,
    grange,
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
    avg_method <- match.arg(avg_method)

    assert_that(
        is(grange, "GRanges"),
        length(grange) == 1
    )

    chr <- as.character(GenomicRanges::seqnames(grange))
    start <- GenomicRanges::start(grange)
    end <- GenomicRanges::end(grange)

    plot_region(
        x,
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

#' Plot GRanges heatmap
#'
#' @param x the NanoMethResult object.
#' @param grange the GRanges object with one entry.
#' @param pos_style the style for plotting the base positions along the x-axis.
#'   Defaults to "to_scale", plotting (potentially) overlapping squares
#'   along the genomic position to scale. The "compact" options plots only the
#'   positions with measured modification.
#' @param window_prop the size of flanking region to plot. Can be a vector of two
#'   values for left and right window size. Values indicate proportion of gene
#'   length.
#'
#' @return a ggplot plot containing the heatmap.
#'
#' @examples
#' nmr <- load_example_nanomethresult()
#' gr <- GenomicRanges::GRanges(data.frame(chr = "chr7", start = 6703892, end = 6730431))
#' plot_grange_heatmap(nmr, gr[1, ])
#'
#' @export
plot_grange_heatmap <- function(
    x,
    grange,
    pos_style = c("to_scale", "compact"),
    window_prop = 0
) {

    assert_that(
        is(grange, "GRanges"),
        length(grange) == 1
    )

    chr <- as.character(GenomicRanges::seqnames(grange))
    start <- GenomicRanges::start(grange)
    end <- GenomicRanges::end(grange)

    plot_region_heatmap(
        x,
        chr = chr,
        start = start,
        end = end,
        pos_style = pos_style,
        window_prop = window_prop
    )
}
