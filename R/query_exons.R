#' Query exons
#'
#' Query a data.frame, NanoMethResult or ModBamResult for exon annotation.
#'
#' @param x the object to query.
#' 
#' @return data.frame of queried exons.
#'
#' @name query_exons
NULL

#' @param chr the chromosome to query.
#' @param start the start of the query region.
#' @param end the end of the query region.
#'
#' @export
#' @describeIn query_exons Query region.
query_exons_region <- function(x, chr, start, end) {
    exons <- get_exons(x)

    genes_df <- exons %>%
        dplyr::group_by(.data$gene_id, .data$chr) %>%
        dplyr::summarise(
            start = min(start),
            end = max(end)
        )

    genes <- genes_df %>%
        dplyr::filter(
            start <= !!end & end >= !!start,
            chr == !!chr
        ) %>%
        dplyr::select("gene_id")

    exons %>%
        dplyr::inner_join(genes, by = "gene_id", multiple = "all")
}

#' @param gene_id the gene_id to query.
#'
#' @export
#' @describeIn query_exons Query gene ID.
query_exons_gene_id <- function(x, gene_id) {
    exons <- get_exons(x)

    out <- exons %>%
        dplyr::filter(gene_id %in% !!gene_id)

    if (nrow(out) == 0) {
        stop(glue::glue("gene {gene_id} not found in exon annotation"))
    }

    out
}

#' @param symbol the gene_id to query.
#'
#' @export
#' @describeIn query_exons Query gene symbol.
query_exons_symbol <- function(x, symbol) {
    exons <- get_exons(x)

    out <- exons %>%
        dplyr::filter(symbol %in% !!symbol)

    if (nrow(out) == 0) {
        stop(glue::glue("gene {symbol} not found in exon annotation"))
    }

    out
}

get_exons <- function(x) {
    if (is(x, "NanoMethResult")) {
        exons(x)
    } else if (is(x, "ModBamResult")) {
        exons(x)
    } else if (is(x, "data.frame")) {
        x
    } else {
        stop("x must be a NanoMethResult, ModBamResult or data.frame")
    }
}