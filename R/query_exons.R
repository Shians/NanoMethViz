#' Query exons
#'
#' Query a data.frame of exons for a subset.
#'
#' @return data.frame of queried exons.
#'
#' @name query_exons
NULL

#' @param exons the data.frame of exons.
#' @param chr the chromosome to query.
#' @param start the start of the query region.
#' @param end the end of the query region.
#'
#' @export
#' @describeIn query_exons Query region.
query_exons_region <- function(exons, chr, start, end) {
    genes_df <- exons %>%
        dplyr::group_by(gene_id, chr) %>%
        dplyr::summarise(
            start = min(start),
            end = max(end)
        )

    genes <- genes_df %>%
        dplyr::filter(
            start <= !!end & end >= !!start,
            chr == !!chr
        ) %>%
        dplyr::select(gene_id)

    exons %>%
        dplyr::inner_join(genes, by = "gene_id")
}

#' @param exons the data.frame of exons.
#' @param gene_id the gene_id to query.
#'
#' @export
#' @describeIn query_exons Query gene ID.
query_exons_gene_id <- function(exons, gene_id) {
    out <- exons %>%
        dplyr::filter(gene_id %in% !!gene_id)

    if (nrow(out) == 0) {
        stop(glue::glue("gene {gene_id} not found in exon annotation"))
    }

    out
}

#' @param exons the data.frame of exons.
#' @param symbol the gene_id to query.
#'
#' @export
#' @describeIn query_exons Query gene symbol.
query_exons_symbol <- function(exons, symbol) {
    out <- exons %>%
        dplyr::filter(symbol %in% !!symbol)

    if (nrow(out) == 0) {
        stop(glue::glue("gene {symbol} not found in exon annotation"))
    }

    out
}

