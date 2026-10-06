gene_to_region <- function(x, gene) {
    validate_gene_symbol(x, gene)

    pos_range <- gene_pos_range(x, gene)

    chr <- exons(x) %>%
        dplyr::filter(.data$symbol == gene) %>%
        dplyr::slice(1) %>%
        dplyr::pull(chr)

    list(chr = chr, start = pos_range[1], end = pos_range[2])
}

plot_gene_heatmap_impl <- function(
    x,
    gene,
    window_prop = 0.3,
    pos_style = c("to_scale", "compact"),
    subsample = 50
) {
    pos_style <- match.arg(pos_style)
    region <- gene_to_region(x, gene)

    plot_region_heatmap_impl(
        x = x,
        chr = region$chr,
        start = region$start,
        end = region$end,
        window_prop = window_prop,
        pos_style = pos_style,
        subsample = subsample
    )
}

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
#' @details
#' This function creates a heatmap visualisation of methylation data for a specific gene.
#' Each row in the heatmap represents one or more packed reads, where colored segments
#' indicate methylation probability at each genomic position.
#'
#' @examples
#' nmr <- load_example_nanomethresult()
#' plot_gene_heatmap(nmr, "Peg3")
#'
#' @export
setMethod(
    "plot_gene_heatmap",
    signature(x = "NanoMethResult", gene = "character"),
    plot_gene_heatmap_impl
)

#' @rdname plot_gene_heatmap
#'
#' @export
setMethod(
    "plot_gene_heatmap",
    signature(x = "ModBamResult", gene = "character"),
    plot_gene_heatmap_impl
)
