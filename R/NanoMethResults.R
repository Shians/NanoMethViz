#' Nanopore Methylation Result
#'
#' @slot methy the path to methylation database (sql or indexed tabix).
#' @slot samples the data.frame of sample annotation containg at least columns
#'   sample and group.
#' @slot exons the data.frame of exon information containing at least columns
#'   gene_id, chr, strand, start, end, transcript_id and symbol.
#'
#' @return a NanoMethResult object to be used with plotting functions
#' @export
setClass(
    "NanoMethResult",
    representation(
        methy = "character",
        samples = "data.frame",
        exons = "data.frame"
    )
)

#' Nanopore Methylation Result Constructor
#'
#' @param methy the path to methylation database (sql or indexed tabix).
#' @param samples the data.frame of sample annotation containg at least columns
#'   sample and group.
#' @param exons the data.frame of exon information containing at least columns
#'   gene_id, chr, strand, start, end, transcript_id and symbol.
#'
#' @return a NanoMethResult object to be used with plotting functions
#' @export
NanoMethResult <- function(methy, samples, exons = NULL) {
    if (is.null(exons)) {
        exons <- tibble::tibble(
            gene_id = character(),
            chr = character(),
            strand = character(),
            start = integer(),
            end = integer(),
            transcript_id = character(),
            symbol = character()
        )
    }

    new(
        "NanoMethResult",
        methy = methy,
        samples = tibble::as_tibble(samples),
        exons = tibble::as_tibble(exons)
    )
}


#' Get methylation data
#'
#' @param object
#'
#' @return
#' @export
setGeneric("methy", valueClass = "character", function(object) {
    standardGeneric("methy")
})
setMethod("methy", signature("NanoMethResult"), function(object)
{
    object@methy
})

#' Get sample annotation
#'
#' @param object
#'
#' @return
#' @export
setGeneric("samples", valueClass = "data.frame", function(object) {
    standardGeneric("samples")
})
setMethod("samples", signature("NanoMethResult"), function(object)
{
    object@samples
})

#' Get exon annotation
#'
#' @param object
#'
#' @return
#' @export
setGeneric("exons", valueClass = "data.frame", function(object) {
    standardGeneric("exons")
})
setMethod("exons", signature("NanoMethResult"), function(object)
{
    object@exons
})