#' Get exon annotations for mus musculus
#'
#' @return data.frame containing exons
#'
#' @examples
#' m_musculus_exons <- get_exons_mus_musculus()
#'
#' @importFrom readr read_tsv
#' @importFrom tibble as_tibble
#' @importFrom dplyr rename
#' @importFrom utils install.packages packageVersion
#' @export
get_exons_mus_musculus <- function() {
    genes <-  AnnotationDbi::keys(Mus.musculus::Mus.musculus, "GENEID")
    exon_data <- suppressMessages(AnnotationDbi::select(
        Mus.musculus::Mus.musculus,
        keys = genes,
        keytype = "GENEID",
        columns = c(
            "GENEID",
            "TXID",
            "EXONCHROM",
            "EXONSTRAND",
            "EXONSTART",
            "EXONEND",
            "SYMBOL"
        )
    ))

    tibble::as_tibble(exon_data) %>%
        dplyr::rename(
            gene_id = "GENEID",
            chr = "EXONCHROM",
            strand = "EXONSTRAND",
            start = "EXONSTART",
            end = "EXONEND",
            transcript_id = "TXID",
            symbol = "SYMBOL"
        )
}

#' Get example exon annotations for mus musculus
#'
#' This is a small subset of the exons returned by
#' \code{get_exons_mus_musculus()} for demonstrative purposes. It contains
#' the exons for the genes Brca1, Brca2, Impact, Meg3, Peg3 and Xist.
#'
#' @return data.frame containing exons
#'
#' @examples
#' example_exons <- get_example_exons_mus_musculus()
#'
#' @export
get_example_exons_mus_musculus <- function() {
    genes <-  c("12189", "12190", "16210", "17263", "18616", "213742")
    exon_data <- suppressMessages(AnnotationDbi::select(
        Mus.musculus::Mus.musculus,
        keys = genes,
        keytype = "GENEID",
        columns = c(
            "GENEID",
            "TXID",
            "EXONCHROM",
            "EXONSTRAND",
            "EXONSTART",
            "EXONEND",
            "SYMBOL"
        )
    ))

    tibble::as_tibble(exon_data) %>%
        dplyr::rename(
            gene_id = "GENEID",
            chr = "EXONCHROM",
            strand = "EXONSTRAND",
            start = "EXONSTART",
            end = "EXONEND",
            transcript_id = "TXID",
            symbol = "SYMBOL"
        )
}

#' Get exon annotations for homo sapiens
#'
#' @return data.frame containing exons
#'
#' @examples
#' h_sapiens_exons <- get_exons_homo_sapiens()
#'
#' @importFrom readr read_tsv
#' @importFrom tibble as_tibble
#' @importFrom dplyr rename
#' @importFrom utils install.packages packageVersion
#' @export
get_exons_homo_sapiens <- function() {
    genes <-  AnnotationDbi::keys(Homo.sapiens::Homo.sapiens, "GENEID")
    exon_data <- suppressMessages(AnnotationDbi::select(
        Homo.sapiens::Homo.sapiens,
        keys = genes,
        keytype = "GENEID",
        columns = c(
            "GENEID",
            "TXID",
            "EXONCHROM",
            "EXONSTRAND",
            "EXONSTART",
            "EXONEND",
            "SYMBOL"
        )
    ))

    tibble::as_tibble(exon_data) %>%
        dplyr::rename(
            gene_id = "GENEID",
            chr = "EXONCHROM",
            strand = "EXONSTRAND",
            start = "EXONSTART",
            end = "EXONEND",
            transcript_id = "TXID",
            symbol = "SYMBOL"
        )

}
