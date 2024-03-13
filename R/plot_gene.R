#' @rdname plot_gene
#'
#' @inheritParams plot_region
#' @param gene_anno whether to show gene annotation.
#'
#' @details This function plots the methylation data for a given gene. Since
#' V3.0.0 NanoMethViz has changed the smoothing strategy from a loess smoothing
#' to a weighted moving average. This is because the loess smoothing was too
#' computationally expensive for large datasets and had a span parameter that
#' was difficult to tune. The new smoothing strategy is controlled by the
#' smoothing_window argument.
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

        plot_gene_impl(
            x,
            gene,
            window_prop = window_prop,
            anno_regions = anno_regions,
            binary_threshold = binary_threshold,
            avg_method = avg_method,
            spaghetti = spaghetti,
            heatmap = heatmap,
            heatmap_subsample = heatmap_subsample,
            smoothing_window = smoothing_window,
            gene_anno = gene_anno,
            palette = palette,
            line_size = line_size,
            mod_scale = mod_scale
        )
    }
)

#' @describeIn plot_gene S4 method for ModBamResult
setMethod("plot_gene", signature(x = "ModBamResult", gene = "character"),
    function(
        x,
        gene,
        window_prop = 0.3,
        anno_regions = NULL,
        binary_threshold = NULL,
        avg_method = c("mean", "median"),
        spaghetti = FALSE,
        heatmap = FALSE,
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

        plot_gene_impl(
            x,
            gene,
            window_prop = window_prop,
            anno_regions = anno_regions,
            binary_threshold = binary_threshold,
            avg_method = avg_method,
            spaghetti = spaghetti,
            heatmap = heatmap,
            heatmap_subsample = heatmap_subsample,
            smoothing_window = smoothing_window,
            gene_anno = gene_anno,
            palette = palette,
            line_size = line_size,
            mod_scale = mod_scale
        )
    }
)

plot_gene_impl <- function(
    x,
    gene,
    window_prop,
    anno_regions,
    binary_threshold,
    avg_method,
    spaghetti,
    heatmap,
    heatmap_subsample,
    smoothing_window,
    gene_anno,
    palette,
    line_size,
    mod_scale = c(0, 1)
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
    exons_anno <- query_exons_symbol(x, symbol = gene)

    feature <- list()
    feature_chr <- unique(exons_anno$chr)
    feature_start <- min(exons_anno$start)
    feature_end <- max(exons_anno$end)
    window_size <- round((feature_end - feature_start) * window_prop)

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
