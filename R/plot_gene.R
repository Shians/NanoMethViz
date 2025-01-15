plot_gene_impl <- function(
        x,
        gene,
        window_prop = 0.3,
        anno_regions = NULL,
        binary_threshold = NULL,
        avg_method = c("mean", "median"),
        spaghetti = FALSE,
        heatmap = TRUE,
        heatmap_subsample = 50,
        smoothing_window = 2000,
        gene_anno = TRUE,
        palette = ggplot2::scale_colour_brewer(palette = "Set1"),
        line_size = 1,
        mod_scale = c(0, 1),
        span = NULL
) {
    if (!missing("span")) {
        warning("the 'span' argument has been deprecated, please use 'smoothing_window' instead")
    }

    avg_method <- match.arg(avg_method)

    assertthat::assert_that(
        nrow(exons(x)) > 0,
        msg = "exons(x) is empty, gene cannot be queried"
    )

    if (length(window_prop) == 1) {
        # convert to two sided window_prop
        window_prop <- c(window_prop, window_prop)
    }

    exons_anno <- query_exons_symbol(x, symbol = gene)

    feature_chr <- unique(exons_anno$chr)
    feature_start <- min(exons_anno$start)
    feature_end <- max(exons_anno$end)

    plot_region(
        x = x,
        chr = feature_chr,
        start = feature_start,
        end = feature_end,
        anno_regions = anno_regions,
        binary_threshold = binary_threshold,
        avg_method = avg_method,
        spaghetti = spaghetti,
        heatmap = heatmap,
        heatmap_subsample = heatmap_subsample,
        smoothing_window = smoothing_window,
        gene_anno = gene_anno,
        window_prop = window_prop,
        palette = palette,
        line_size = line_size,
        mod_scale = mod_scale
    )
}

#' @rdname plot_gene
#'
#' @inheritParams plot_region
#'
#' @details
#' This function plots the methylation data for a given gene. The main trendline plot shows the average methylation
#' probability across the gene. The heatmap plot shows the methylation probability for each read across the gene. The
#' gene annotation plot shows the exons of the gene. In the heatmap, each row represents one or more non-overlapping reads
#' where the coloured segments represent the methylation probability at each position. Data along a read is connected by
#' a grey line. The gene annotation plot shows the isoforms and exons of the gene, with arrows indicating the direction
#' of transcription.
#'
#' Since V3.0.0 NanoMethViz has changed the smoothing strategy from a loess smoothing to a weighted moving average. This
#' is because the loess smoothing was too computationally expensive for large datasets and had a span parameter that was
#' difficult to tune. The new smoothing strategy is controlled by the smoothing_window argument.
#'
#' @examples
#' nmr <- load_example_nanomethresult()
#' plot_gene(nmr, "Peg3")
#'
#' @export
setMethod("plot_gene", signature(x = "NanoMethResult", gene = "character"),
    plot_gene_impl
)

#' @describeIn plot_gene S4 method for ModBamResult
setMethod("plot_gene", signature(x = "ModBamResult", gene = "character"),
    plot_gene_impl
)
