#' Plot gene
#'
#' @param x the NanoMethResults object.
#' @param gene the gene symbol for the gene to plot.
#' @param ... additional arguments
#'
#' @return None
#' @export
setGeneric("plot_gene", function(x, gene, ...) {
  standardGeneric("plot_gene")
})

#' @rdname plot_gene
#'
#' @param anno_regions the data.frame of regions to annotate.
#' @param spaghetti whether or not individual reads should be shown.
#' @param span the span for loess smoothing.
#'
#' @export
setMethod("plot_gene", signature(x = "NanoMethResult", gene = "character"),
    function(x, gene, anno_regions = NULL, spaghetti = FALSE, span = NULL) {
        .plot_gene(x, gene, anno_regions = anno_regions, spaghetti = spaghetti, span = span)
    }
)

.plot_gene <- function(
    x,
    gene,
    anno_regions = NULL,
    spaghetti = FALSE,
    span = NULL
    ) {
    assertthat::assert_that(
      nrow(exons(x)) > 0,
      msg = "exons(x) is empty, gene cannot be queried"
    )

    sample_anno <- samples(x)
    exons_anno <- query_exons_symbol(exons(x), symbol = gene)

    feature <- list()
    feature$chr <- unique(exons_anno$chr)
    feature$start <- min(exons_anno$start)
    feature$end <- max(exons_anno$end)

    p1 <- with(exons_anno,
         plot_feature(
            feature,
            title = gene,
            methy = methy(x),
            sample_anno = sample_anno,
            anno_regions = anno_regions,
            spaghetti = spaghetti,
            span = span
        )
    )

    xlim <- .get_ggplot_range_x(p1)

    p2 <- plot_gene_annotation(exons_anno, xlim[1], xlim[2])

    n_unique <- function(x) { length(unique(x)) }

    heights <- c(1, 0.075 * n_unique(exons_anno$transcript_id))
    p1 / p2 + patchwork::plot_layout(heights = heights)
}
