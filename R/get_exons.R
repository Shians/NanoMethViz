#' Get exon annotations for mus musculus
#'
#' @return data.frame containing exons
#'
#' @importFrom readr read_tsv
#' @export
get_exons_mus_musculus <- function() {
    url <- "https://github.com/Shians/exon_annotations/raw/master/mus_musculus_exons.tsv.gz"

    readr::read_tsv(
        url,
        col_types = readr::cols(
            gene_id = readr::col_double(),
            chr = readr::col_character(),
            strand = readr::col_character(),
            start = readr::col_double(),
            end = readr::col_double(),
            transcript_id = readr::col_double(),
            symbol = readr::col_character()
        )
    )
}

#' Get exon annotations for homo sapiens
#'
#' @return data.frame containing exons
#'
#' @importFrom readr read_tsv
#' @export
get_exons_homo_sapiens <- function() {
    url <- "https://github.com/Shians/exon_annotations/raw/master/homo_sapiens_exons.tsv.gz"

    readr::read_tsv(
        url,
        col_types = readr::cols(
            gene_id = readr::col_double(),
            chr = readr::col_character(),
            strand = readr::col_character(),
            start = readr::col_double(),
            end = readr::col_double(),
            transcript_id = readr::col_double(),
            symbol = readr::col_character()
        )
    )
}
