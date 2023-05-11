#' @rdname plot_gene_heatmap
#'
#' @param window_prop the size of flanking region to plot. Can be a vector of two
#'   values for left and right window size. Values indicate proportion of gene
#'   length.
#' @param pos_style the style for plotting the base positions along the x-axis.
#'   Defaults to "to_scale", plotting (potentially) overlapping squares
#'   along the genomic position to scale. The "compact" options plots only the
#'   positions with measured modification.
#' @param subsample the number of read of packed read rows to subsample to.
#'
#' @return a ggplot plot containing the heatmap.
#'
#' @examples
#' nmr <- load_example_nanomethresult()
#' plot_gene_heatmap(nmr, "Peg3")
#'
#' @importFrom scico scale_colour_scico
#' @export
setMethod(
    "plot_gene_heatmap",
    signature(x = "NanoMethResult", gene = "character"),
    function(
        x,
        gene,
        window_prop = 0.3,
        pos_style = c("to_scale", "compact"),
        subsample = 50
    ) {
        pos_style <- match.arg(pos_style)

        .plot_gene_heatmap(
            x = x,
            gene = gene,
            window_prop = window_prop,
            pos_style = pos_style,
            subsample = subsample
        )
    }
)

#' @rdname plot_gene_heatmap
#'
#' @export
setMethod(
    "plot_gene_heatmap",
    signature(x = "ModBamResult", gene = "character"),
    function(
        x,
        gene,
        window_prop = 0.3,
        pos_style = c("to_scale", "compact"),
        subsample = 50
    ) {
        pos_style <- match.arg(pos_style)

        .plot_gene_heatmap(
            x = x,
            gene = gene,
            window_prop = window_prop,
            pos_style = pos_style,
            subsample = subsample
        )
    }
)

.plot_gene_heatmap <- function(
    x,
    gene,
    window_prop,
    pos_style,
    subsample
) {
    assertthat::assert_that(
        nrow(exons(x)) > 0,
        msg = "exons(x) is empty, gene cannot be queried"
    )

    if (length(window_prop) == 1) {
        # convert to two sided window
        window_prop <- c(window_prop, window_prop)
    }

    # query_gene_methy
    if (!gene %in% exons(x)$symbol) {
        stop(glue::glue("gene {gene} not found in exon annotation"))
    }

    pos_range <- gene_pos_range(x, gene)

    chr <- exons(x) %>%
        dplyr::filter(.data$symbol == gene) %>%
        dplyr::slice(1) %>%
        dplyr::pull(chr)

    .plot_region_heatmap(
        x = x,
        chr = chr,
        start = pos_range[1],
        end = pos_range[2],
        window_prop = window_prop,
        pos_style = pos_style,
        subsample = subsample
    )
}
