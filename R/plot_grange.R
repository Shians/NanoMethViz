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
    spaghetti = FALSE,
    span = NULL
) {

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
        spaghetti = spaghetti,
        span = span
    )
}
