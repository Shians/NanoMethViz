#' Plot region
#'
#' @param x the NanoMethResult object.
#' @param chr the chromosome to plot.
#' @param start the start of the plotting region.
#' @param end the end of the plotting region.
#' @param ... additional arguments.
#'
#' @return a patchwork plot containing the methylation profile in the specified
#'   region.
#'
#' @examples
#' nmr <- load_example_nanomethresult()
#' plot_region(nmr, "chr7", 6703892, 6730431)
#'
#' @export
setGeneric("plot_region", function(x, chr, start, end, ...) {
    standardGeneric("plot_region")
})

#' @rdname plot_region
#'
#' @param anno_regions the data.frame of regions to be annotated.
#' @param spaghetti whether or not individual reads should be shown.
#' @param span the span for loess smoothing.
#' @param window_prop the size of flanking region to plot. Can be a vector of two
#'   values for left and right window size. Values indicate proportion of gene
#'   length.
#'
#' @return a patchwork plot containing the methylation profile in the specified
#'   region.
#'
#' @examples
#' nmr <- load_example_nanomethresult()
#' plot_region(nmr, "chr7", 6703892, 6730431)
#'
#' @export
setMethod("plot_region",
    signature(
        x = "NanoMethResult",
        chr = "character",
        start = "numeric",
        end = "numeric"),
    function(
        x,
        chr,
        start,
        end,
        anno_regions = NULL,
        spaghetti = FALSE,
        span = NULL,
        window_prop = 0
    ) {
        .plot_region(
            x = x,
            chr = chr,
            start = start,
            end = end,
            anno_regions = anno_regions,
            spaghetti = spaghetti,
            span = span,
            window_prop = window_prop
        )
    }
)

#' @rdname plot_region
#'
#' @export
setMethod("plot_region",
    signature(
        x = "NanoMethResult",
        chr = "factor",
        start = "numeric",
        end = "numeric"),

    function(
        x,
        chr,
        start,
        end,
        anno_regions = NULL,
        spaghetti = FALSE,
        span = NULL,
        window_prop = 0
    ) {
        chr <- as.character(chr)
        plot_region(
            x = x,
            chr = chr,
            start = start,
            end = end,
            anno_regions = anno_regions,
            spaghetti = spaghetti,
            span = span,
            window_prop = window_prop
        )
    }
)


.plot_region <- function(
    x,
    chr,
    start,
    end,
    anno_regions = NULL,
    spaghetti = FALSE,
    span = NULL,
    window_prop = 0
) {
    sample_anno <- samples(x)
    exons_anno <- query_exons_region(
        exons(x),
        chr = chr,
        start = start,
        end = end
    )

    if (length(window_prop) == 1) {
        window_prop <- c(window_prop, window_prop)
    }

    feature_width <- end - start
    window_left <- feature_width * window_prop[1]
    window_right <- feature_width * window_prop[2]
    xlim <- round(c(start - window_left, end + window_right))

    methy_data <-
        query_methy(
            x,
            chr,
            floor(start - window_left * 1.1),
            ceiling(end + window_right * 1.1),
            simplify = TRUE) %>%
        dplyr::select(-"strand") %>%
        tibble::as_tibble()

    if (nrow(methy_data) == 0) {
        warning("no methylation data in region")
        return(ggplot() + theme_void())
    }

    title <- glue::glue("{chr}:{start}-{end}")
    p1 <- plot_methylation_internal(
        methy_data = methy_data,
        start = start,
        end = end,
        chr = chr,
        xlim = xlim,
        title = title,
        anno_regions = anno_regions,
        spaghetti = spaghetti,
        sample_anno = sample_anno,
        span = span
    )

    xlim <- .get_ggplot_range_x(p1)

    p2 <- plot_gene_annotation(exons_anno, xlim[1], xlim[2])

    anno_height <- attr(p2, "plot_height")

    heights <- c(1, 0.075 * anno_height)
    p1 / p2 + patchwork::plot_layout(heights = heights)
}
