setGeneric("plot_region", function(x, gene, ...) {
  standardGeneric("plot_region")
})

setMethod("plot_region", signature(x = "NanoMethResult", gene = "character"),
    function(x, gene, ...) {
        plot_region(x, gene, ...)
    }
)

#' Plot region
#'
#' @param x the NanoMethResults object
#' @param chr the chromosome to plot
#' @param start the start of the plotting region
#' @param end the end of the plotting region
#' @param anno_regions the data.frame of regions to be annotated
#' @param spaghetti whether spaghettis should be drawn
#'
#' @return
#' @export
plot_region <- function(
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
