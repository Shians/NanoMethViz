#' Plot gene
#'
#' @param x the NanoMethResult object.
#' @param gene the gene symbol for the gene to plot.
#' @param ... additional arguments
#'
#' @return a patchwork plot containing the methylation profile in the specified
#'   region.
#'
#' @examples
#' nmr <- load_example_nanomethresult()
#' plot_gene(nmr, "Peg3")
#'
#' @export
setGeneric("plot_gene", function(x, gene, ...) {
    standardGeneric("plot_gene")
})

#' @rdname plot_gene
#'
#' @param window_prop the size of flanking region to plot. Can be a vector of two
#'   values for left and right window size. Values indicate proportion of gene
#'   length.
#' @param anno_regions the data.frame of regions to annotate.
#' @param binary_threshold the modification probability such that calls with
#'   modification probability above the threshold are set to 1 and probabilities
#'   equalt to or below the threshold are set to 0.
#' @param spaghetti whether or not individual reads should be shown.
#' @param span the span for loess smoothing.
#' @param gene_anno whether or not gene annotation tracks are plotted.
#'
#' @return a patchwork plot containing the methylation profile in the specified
#'   region.
#'
#' @examples
#' nmr <- load_example_nanomethresult()
#' plot_gene(nmr, "Peg3")
#'
#' @export
setMethod("plot_gene", signature(x = "NanoMethResult", gene = "character"),
    function(
        x,
        gene,
        window_prop = 0.3,
        anno_regions = NULL,
        binary_threshold = NULL,
        spaghetti = FALSE,
        span = NULL,
        gene_anno = TRUE
    ) {
        .plot_gene(
            x,
            gene,
            window_prop = window_prop,
            anno_regions = anno_regions,
            binary_threshold = binary_threshold,
            spaghetti = spaghetti,
            span = span,
            gene_anno = gene_anno
        )
    }
)

.plot_gene <- function(
    x,
    gene,
    window_prop,
    anno_regions,
    binary_threshold,
    spaghetti,
    span,
    gene_anno
) {
    assertthat::assert_that(
        nrow(exons(x)) > 0,
        msg = "exons(x) is empty, gene cannot be queried"
    )

    if (length(window_prop) == 1) {
        # convert to two sided window_prop
        window_prop <- c(window_prop, window_prop)
    }

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
            window_prop = window_prop,
            sample_anno = sample_anno,
            anno_regions = anno_regions,
            binary_threshold = binary_threshold,
            spaghetti = spaghetti,
            span = span
        )
    )

    if (gene_anno) {
        # if gene annotation is needed
        xlim <- .get_ggplot_range_x(p1)

        p2 <- plot_gene_annotation(exons_anno, xlim[1], xlim[2])

        n_unique <- function(x) { length(unique(x)) }

        heights <- c(1, 0.075 * n_unique(exons_anno$transcript_id))
        p_out <- p1 / p2 + patchwork::plot_layout(heights = heights)
    } else {
        # if no gene annotation
        p_out <- p1
    }

    p_out
}
