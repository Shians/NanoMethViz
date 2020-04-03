query_exons_region <- function(exons, chr, start, end) {
    message(glue::glue("querying {chr}:{scales::comma(start)}-{scales::comma(end)}"))
    exons %>%
        dplyr::filter(
            start >= !!start & end <= !!end,
            chr == !!chr
        )
}

query_exons_gene_id <- function(exons, gene_id) {
    exons %>%
        dplyr::filter(gene_id == !!gene_id)
}

query_exons_symbol <- function(exons, symbol) {
    exons %>%
        dplyr::filter(symbol == !!symbol)
}

