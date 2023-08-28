#' @rdname plot_gene
#'
#' @inheritParams plot_region
#' @param gene_anno whether to show gene annotation.
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
        span = NULL,
        gene_anno = TRUE,
        palette = ggplot2::scale_colour_brewer(palette = "Set1"),
        line_size = 1,
        mod_scale = c(0, 1)
    ) {
        avg_method <- match.arg(avg_method)
        .plot_gene(
            x,
            gene,
            window_prop = window_prop,
            anno_regions = anno_regions,
            binary_threshold = binary_threshold,
            avg_method = avg_method,
            spaghetti = spaghetti,
            heatmap = heatmap,
            heatmap_subsample = heatmap_subsample,
            span = span,
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
        span = NULL,
        gene_anno = TRUE,
        palette = ggplot2::scale_colour_brewer(palette = "Set1"),
        line_size = 1,
        mod_scale = c(0, 1)
    ) {
        avg_method <- match.arg(avg_method)
        .plot_gene(
            x,
            gene,
            window_prop = window_prop,
            anno_regions = anno_regions,
            binary_threshold = binary_threshold,
            avg_method = avg_method,
            spaghetti = spaghetti,
            heatmap = heatmap,
            heatmap_subsample = heatmap_subsample,
            span = span,
            gene_anno = gene_anno,
            palette = palette,
            line_size = line_size,
            mod_scale = mod_scale
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
    heatmap_subsample,
    span,
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
            methy = x,
            window_size = window_size,
            sample_anno = sample_anno,
            anno_regions = anno_regions,
            binary_threshold = binary_threshold,
            avg_method = avg_method,
            spaghetti = spaghetti,
            span = span,
            palette = palette,
            line_size = line_size,
            mod_scale = mod_scale
        )
    )
    p1 <- p1 + ggplot2::coord_cartesian(
        xlim = c(plot_left, plot_right),
        expand = FALSE
    )


    if (gene_anno) {
        # if gene annotation is needed
        xlim <- .get_ggplot_range_x(p1)

        p2 <- plot_gene_annotation(exons_anno, plot_left, plot_right) +
            ggplot2::coord_cartesian(
                xlim = c(plot_left, plot_right),
                expand = FALSE
            ) +
                ggplot2::scale_x_continuous(labels = scales::label_number(scale_cut = scales::cut_si("b")))

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
            window_prop,
            subsample = heatmap_subsample
        ) +
            ggplot2::coord_cartesian(
                xlim = c(plot_left, plot_right)
            )

        p_out <- stack_plots(p_out, ggrastr::rasterise(p_heatmap))
    }

    p_out
}
