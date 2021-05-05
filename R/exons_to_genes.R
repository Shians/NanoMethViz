#' Convert exon annotation to genes
#'
#' @param x the exon level annotation containing columns "gene_id", "chr",
#'   "strand" and "symbol".
#'
#' @return the gene level annotation where each gene is taken to span the
#'   earliest start position and latest end position of its exons.
#' @export
#'
#' @examples
#' nmr <- load_example_nanomethresult()
#' exons_to_genes(NanoMethViz::exons(nmr))
exons_to_genes <- function(x) {
    x %>%
        dplyr::group_by(.data$gene_id, .data$chr, .data$strand, .data$symbol) %>%
        dplyr::summarise(start = min(.data$start), end = max(.data$end), .groups = "drop")
}
