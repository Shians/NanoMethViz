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
#'
#' @param anno_regions the data.frame of regions to be annotated
#' @param spaghetti whether or not individual reads should be shown.
#' @param span the span for loess smoothing.
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
        span = NULL
        ) {
        .plot_region(
            x = x,
            chr = chr,
            start = start,
            end = end,
            anno_regions = anno_regions,
            spaghetti = spaghetti,
            span = span
        )
    })

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
        span = NULL
        ) {
        chr <- as.character(chr)
        plot_region(
            x = x,
            chr = chr,
            start = start,
            end = end,
            anno_regions = anno_regions,
            spaghetti = spaghetti,
            span = span
        )
    })


.plot_region <- function(
    x,
    chr,
    start,
    end,
    anno_regions = NULL,
    spaghetti = FALSE,
    span = NULL,
    window_prop = c(0.3, 0.3)
    ) {
    sample_anno <- samples(x)
    exons_anno <- query_exons_region(
        exons(x),
        chr = chr,
        start = start,
        end = end
    )

    window_left <- (end - start) * window_prop[1]
    window_right <- (end - start) * window_prop[2]

    methy_data <-
        query_methy(methy(x), chr, start - window_left, end + window_right) %>%
        dplyr::bind_rows() %>%
        dplyr::select(-"strand", -"modified") %>%
        tibble::as_tibble()

    if (nrow(methy_data) == 0) {
        warning("no methylation data in region")
        return(ggplot() + theme_void())
    }

    p1 <- with(
        exons_anno,
        plot_methylation_internal(
            methy_data = methy_data,
            start = start,
            end = end,
            chr = chr,
            title = glue::glue("{chr}:{start}-{end}"),
            anno_regions = anno_regions,
            spaghetti = spaghetti,
            sample_anno = sample_anno,
            span = span
        )
    )

    xlim <- .get_ggplot_range_x(p1)

    p2 <- plot_gene_annotation(exons_anno, xlim[1], xlim[2])

    anno_height <- attr(p2, "plot_height")

    heights <- c(1, 0.075 * anno_height)
    p1 / p2 + patchwork::plot_layout(heights = heights)
}
