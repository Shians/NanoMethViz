#' Plot region
#'
#' @param x the NanoMethResults object
#' @param chr the chromosome to plot
#' @param start the start of the plotting region
#' @param end the end of the plotting region
#' @param ... additional arguments
#'
#' @return None
#' @export
setGeneric("plot_region", function(x, chr, start, end, ...) {
  standardGeneric("plot_region")
})

#' @rdname plot_region
#' @param anno_regions the data.frame of regions to be annotated
#' @param spaghetti whether or not individual reads should be shown.
#'
#' @export
setMethod("plot_region", signature(x = "NanoMethResult", chr = "character", start = "numeric", end = "numeric"),
    function(x, chr, start, end, anno_regions, spaghetti) {
      .plot_region(
        x = x,
        chr = chr,
        start = start,
        end = end,
        anno_regions = anno_regions,
        spaghetti = spaghetti
      )
    }
)


.plot_region <- function(
    x,
    chr,
    start,
    end,
    anno_regions = NULL,
    spaghetti = FALSE
    ) {
    sample_anno <- samples(x)
    exons_anno <- query_exons_region(exons(x), chr = chr, start = start, end = end)

    feature <- list()
    feature$chr <- unique(exons_anno$chr)
    feature$start <- min(exons_anno$start)
    feature$end <- max(exons_anno$end)

    p1 <- with(exons_anno,
         plot_feature(
            feature,
            title = glue::glue("{chr}:{start}-{end}"),
            methy = methy(x),
            sample_anno = sample_anno,
            anno_regions = anno_regions,
            spaghetti = spaghetti
        )
    )

    xlim <- .get_ggplot_range_x(p1)

    p2 <- plot_gene_annotation(exons_anno, xlim[1], xlim[2])

    n_unique <- function(x) { length(unique(x)) }

    heights <- c(1, 0.075 * n_unique(exons_anno$transcript_id))
    p1 / p2 + patchwork::plot_layout(heights = heights)
}
