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
        dplyr::group_by(gene_id, chr, strand, symbol) %>%
        dplyr::summarise(start = min(start), end = max(end), .groups = "drop")
}
