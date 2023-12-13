#' Plot GRanges
#'
#' @param x the NanoMethResult object.
#' @param grange the GRanges object with one entry.
#' @inheritParams plot_region
#' @inherit plot_region return
#'
#' @examples
#' nmr <- load_example_nanomethresult()
#' plot_grange(nmr, GenomicRanges::GRanges("chr7:6703892-6730431"))
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
    heatmap_subsample = 50,
    smoothing_window = 500,
    window_prop = 0,
    palette = ggplot2::scale_colour_brewer(palette = "Set1"),
    line_size = 1,
        span = NULL
) {
    if (!missing("span")) {
        warning("the 'span' argument has been deprecated, please use 'smoothing_window' instead")
    }
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
        heatmap_subsample = heatmap_subsample,
        smoothing_window = smoothing_window,
        window_prop = window_prop,
        palette = palette,
        line_size = line_size
    )
}

#' Plot GRanges heatmap
#'
#' @param x the NanoMethResult object.
#' @param grange the GRanges object with one entry.
#' @inheritParams plot_region_heatmap
#'
#' @return a ggplot plot containing the heatmap.
#'
#' @examples
#' nmr <- load_example_nanomethresult()
#' plot_grange_heatmap(nmr, GenomicRanges::GRanges("chr7:6703892-6730431"))
#'
#' @export
plot_grange_heatmap <- function(
    x,
    grange,
    pos_style = c("to_scale", "compact"),
    window_prop = 0,
    subsample = 50
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
        window_prop = window_prop,
        subsample = subsample
    )
}
