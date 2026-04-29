#' Get exon annotations
#'
#' Helper functions are provided for obtaining exon annotations from relevant
#' TxDb packages on Bioconductor for the construction of NanoMethResult
#' objects.
#'
#' @name get_exons
#' @rdname get_exons
#'
#' @return tibble (data.frame) object containing exon annotation.
#'
#' @importFrom tibble as_tibble
#' @importFrom dplyr rename
NULL

#' Get exon annotations for Mus musculus (mm10)
#'
#' @return data.frame containing exons
#'
#' @examples
#' m_musculus_exons <- get_exons_mus_musculus()
#'
#' @importFrom tibble as_tibble
#' @importFrom dplyr rename
#' @export
get_exons_mus_musculus <- function() {
    package_check("Mus.musculus", "1.3.1")
    get_exons_from_organism_db(Mus.musculus::Mus.musculus)
}

#' @rdname get_exons
#'
#' @examples
#' mm10_exons <- get_exons_mm10()
#'
#' @export
get_exons_mm10 <- function() {
    package_check(
        c("TxDb.Mmusculus.UCSC.mm10.knownGene", "org.Mm.eg.db"),
        c("3.10.0", "3.15.0")
    )
    get_exons_from_txdb_orgdb(
        TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene,
        org.Mm.eg.db::org.Mm.eg.db
    )
}

#' @rdname get_exons
#'
#' @examples
#' grcm39_exons <- get_exons_grcm39()
#'
#' @export
get_exons_grcm39 <- function() {
    package_check(
        c("TxDb.Mmusculus.UCSC.mm39.refGene", "org.Mm.eg.db"),
        c("3.10.0", "3.15.0")
    )
    get_exons_from_txdb_orgdb(
        TxDb.Mmusculus.UCSC.mm39.refGene::TxDb.Mmusculus.UCSC.mm39.refGene,
        org.Mm.eg.db::org.Mm.eg.db
    )
}

#' Get example exon annotations for mus musculus (mm10)
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
    package_check("Mus.musculus", "1.3.1")
    get_exons_from_organism_db(
        Mus.musculus::Mus.musculus,
        keys = c("12189", "12190", "16210", "17263", "18616", "213742")
    )
}

#' Get exon annotations for Homo sapiens (hg19)
#'
#' @return data.frame containing exons
#'
#' @examples
#' h_sapiens_exons <- get_exons_homo_sapiens()
#'
#' @importFrom tibble as_tibble
#' @importFrom dplyr rename
#' @export
get_exons_homo_sapiens <- function() {
    package_check("Homo.sapiens", "1.3.1")
    get_exons_from_organism_db(Homo.sapiens::Homo.sapiens)
}

#' @rdname get_exons
#'
#' @examples
#' hg19_exons <- get_exons_hg19()
#'
#' @export
get_exons_hg19 <- function() {
    package_check(
        c("org.Hs.eg.db", "TxDb.Hsapiens.UCSC.hg19.knownGene"),
        c("3.15.0", "3.2.2")
    )
    get_exons_from_txdb_orgdb(
        TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene,
        org.Hs.eg.db::org.Hs.eg.db
    )
}

#' @rdname get_exons
#'
#' @examples
#' hg38_exons <- get_exons_hg38()
#'
#' @export
get_exons_hg38 <- function() {
    package_check(
        c("org.Hs.eg.db", "TxDb.Hsapiens.UCSC.hg38.knownGene"),
        c("3.15.0", "3.15.0")
    )
    get_exons_from_txdb_orgdb(
        TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene,
        org.Hs.eg.db::org.Hs.eg.db
    )
}

#' @rdname get_exons
#'
#' @examples
#' t2t_exons <- get_exons_t2t()
#'
#' @export
get_exons_t2t <- function() {
    anno_file <- system.file(
        "exons_t2t.rds",
        package = "NanoMethViz",
        mustWork = TRUE
    )
    readRDS(anno_file) %>%
        filter_primary_chr()
}

get_exons_from_organism_db <- function(db, keys = NULL) {
    if (is.null(keys)) keys <- AnnotationDbi::keys(db, "GENEID")
    suppressMessages(AnnotationDbi::select(
        db,
        keys = keys,
        keytype = "GENEID",
        columns = c("GENEID", "TXID", "EXONCHROM", "EXONSTRAND", "EXONSTART", "EXONEND", "SYMBOL")
    )) %>%
        tibble::as_tibble() %>%
        standardise_exon_col_names() %>%
        filter_primary_chr()
}

get_exons_from_txdb_orgdb <- function(txdb, orgdb) {
    genes <- AnnotationDbi::keys(txdb, "GENEID")
    exon_data <- suppressMessages(AnnotationDbi::select(
        txdb,
        keys = genes,
        keytype = "GENEID",
        columns = c("GENEID", "TXID", "EXONCHROM", "EXONSTRAND", "EXONSTART", "EXONEND")
    ))
    symbols_data <- suppressMessages(AnnotationDbi::select(
        orgdb,
        keys = genes,
        keytype = "ENTREZID",
        columns = "SYMBOL"
    )) %>%
        dplyr::rename(GENEID = "ENTREZID")

    dplyr::left_join(exon_data, symbols_data, by = "GENEID", multiple = "all") %>%
        tibble::as_tibble() %>%
        standardise_exon_col_names() %>%
        filter_primary_chr()
}

standardise_exon_col_names <- function(exon_data) {
    dplyr::rename(exon_data,
        gene_id = "GENEID",
        chr = "EXONCHROM",
        strand = "EXONSTRAND",
        start = "EXONSTART",
        end = "EXONEND",
        transcript_id = "TXID",
        symbol = "SYMBOL"
    )
}

filter_primary_chr <- function(exon_data) {
    dplyr::filter(exon_data, grepl("^chr([0-9]+|[XYM])$", .data$chr))
}

#' @importFrom utils installed.packages packageVersion
package_check <- function(packages, req_versions) {
    assertthat::assert_that(length(packages) == length(req_versions))

    any_missing <- FALSE
    for (i in seq_along(packages)) {
        package <- packages[i]
        req_version <- req_versions[i]
        if (!package %in% utils::installed.packages() || packageVersion(package) < req_version) {
            any_missing <- TRUE
            message(glue::glue(
                "required package '{package} (>= {req_version})' is not installed, please install using BiocManager::install(\"{package}\")"
            ))
        }
    }

    if (any_missing) {
        stop("please install required packages before continuing")
    }
}
