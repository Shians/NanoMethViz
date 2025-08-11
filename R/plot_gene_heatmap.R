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
#' This function creates a heatmap visualization of methylation data for a specific gene.
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
    function(
        x,
        gene,
        window_prop = 0.3,
        pos_style = c("to_scale", "compact"),
        subsample = 50
    ) {
        pos_style <- match.arg(pos_style)

        plot_gene_heatmap_impl(
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

        plot_gene_heatmap_impl(
            x = x,
            gene = gene,
            window_prop = window_prop,
            pos_style = pos_style,
            subsample = subsample
        )
    }
)

plot_gene_heatmap_impl <- function(
    x,
    gene,
    window_prop,
    pos_style,
    subsample
) {
    # validate genes
    if (is.na(gene) || gene == "") {
        stop("Gene symbol cannot be empty or NA. Please provide a valid gene symbol.")
    }

    if (nrow(exons(x)) == 0) {
        stop(glue::glue(
            "No exon annotation found in the data object.\\n",
            "Gene '{gene}' cannot be plotted without exon information.\\n",
            "Please add exon annotation using: exons(your_object) <- your_exons"
        ))
    }

    if (length(window_prop) == 1) {
        # convert to two sided window
        window_prop <- c(window_prop, window_prop)
    }

    # query_gene_methy
    if (!gene %in% exons(x)$symbol) {
        available_genes <- unique(exons(x)$symbol)
        similar_genes <- available_genes[grepl(paste0("^", substr(gene, 1, 3)), available_genes, ignore.case = TRUE)]

        error_msg <- glue::glue(
            "Gene '{gene}' not found in exon annotation.\\n",
            "Please check the gene symbol spelling."
        )

        if (length(similar_genes) > 0 && length(similar_genes) <= 10) {
            error_msg <- paste0(
                error_msg, "\\n",
                glue::glue("Similar genes found: {paste(similar_genes, collapse = ', ')}")
            )
        } else if (length(available_genes) <= 20) {
            error_msg <- paste0(
                error_msg, "\\n",
                glue::glue("Available genes: {paste(available_genes, collapse = ', ')}")
            )
        } else {
            error_msg <- paste0(
                error_msg, "\\n",
                glue::glue("Use unique(exons(your_object)$symbol) to see all {length(available_genes)} available genes.")
            )
        }

        stop(error_msg)
    }

    pos_range <- gene_pos_range(x, gene)

    chr <- exons(x) %>%
        dplyr::filter(.data$symbol == gene) %>%
        dplyr::slice(1) %>%
        dplyr::pull(chr)

    plot_region_heatmap_impl(
        x = x,
        chr = chr,
        start = pos_range[1],
        end = pos_range[2],
        window_prop = window_prop,
        pos_style = pos_style,
        subsample = subsample
    )
}
