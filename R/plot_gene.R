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
#'   equal to or below the threshold are set to 0.
#' @param avg_method the average method for pre-smoothing at each genomic position.
#'   Data is pre-smoothed at each genomic position before the smoothed aggregate line
#'   is generated for performance reasons. The default is "mean" which corresponds
#'   to the average methylation fraction. The alternative "median" option is
#'   closer to an average within the more common methylation state.
#' @param spaghetti whether or not individual reads should be shown.
#' @param heatmap whether or not read-methylation heatmap should be shown.
#' @param span the span for loess smoothing.
#' @param gene_anno whether or not gene annotation tracks are plotted.
#' @param palette the ggplot colour palette used for groups.
#' @param line_size the size of the lines.
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
        avg_method = c("mean", "median"),
        spaghetti = FALSE,
        heatmap = FALSE,
        span = NULL,
        gene_anno = TRUE,
        palette = ggplot2::scale_colour_brewer(palette = "Set1"),
        line_size = 2
    ) {
        avg_method = match.arg(avg_method)
        .plot_gene(
            x,
            gene,
            window_prop = window_prop,
            anno_regions = anno_regions,
            binary_threshold = binary_threshold,
            avg_method = avg_method,
            spaghetti = spaghetti,
            heatmap = heatmap,
            span = span,
            gene_anno = gene_anno,
            palette = palette,
            line_size = line_size
        )
    }
)

.plot_gene <- function(
    x,
    gene,
    window_prop,
    anno_regions,
    binary_threshold,
    avg_method,
    spaghetti,
    heatmap,
    span,
    gene_anno,
    palette,
    line_size
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
    window_size <- round((feature$end - feature$start) * window_prop)

    plot_left <- feature$start - window_size[1]
    plot_right <- feature$end + window_size[2]

    p1 <- with(exons_anno,
        plot_feature(
            feature,
            title = gene,
            methy = methy(x),
            window_size = window_size,
            sample_anno = sample_anno,
            anno_regions = anno_regions,
            binary_threshold = binary_threshold,
            avg_method = avg_method,
            spaghetti = spaghetti,
            span = span,
            palette = palette,
            line_size = line_size
        )
    )
    p1 <- p1 + ggplot2::scale_x_continuous(
            limits = c(plot_left, plot_right),
            expand = ggplot2::expansion()
        )


    if (gene_anno) {
        # if gene annotation is needed
        xlim <- .get_ggplot_range_x(p1)

        p2 <- plot_gene_annotation(exons_anno, plot_left, plot_right) +
            ggplot2::scale_x_continuous(
                limits = c(plot_left, plot_right),
                expand = ggplot2::expansion()
            )

        n_unique <- function(x) { length(unique(x)) }

        heights <- c(1, 0.075 * n_unique(exons_anno$transcript_id))
        p_out <- p1 / p2 + patchwork::plot_layout(heights = heights)
    } else {
        # if no gene annotation
        p_out <- p1
    }

    if (heatmap) {
        p_heatmap <- plot_gene_heatmap(
            x,
            gene,
            window_prop
        ) +
            ggplot2::scale_x_continuous(
                limits = c(plot_left, plot_right),
                expand = ggplot2::expansion())

        p_out <- stack_plots(p_out, p_heatmap)
    }

    p_out
}
