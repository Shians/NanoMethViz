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
#' @importFrom utils install.packages
#' @export
get_exons_mus_musculus <- function() {
    if ("Mus.musculus" %in% utils::installed.packages()) {
        library(Mus.musculus)
    } else {
        stop("package 'Mus.musculus' is not installed, please install using BiocManager::install('Mus.musculus')")
    }

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
#' @importFrom utils install.packages
#' @export
get_exons_homo_sapiens <- function() {
    if ("Homo.sapiens" %in% utils::installed.packages()) {
        library(Homo.sapiens)
    } else {
        stop("package 'Homo.sapiens' is not installed, please install using BiocManager::install('Homo.sapiens')")
    }

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
