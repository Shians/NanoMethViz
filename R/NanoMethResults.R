#' Nanopore Methylation Result
#'
#' A NanoMethResult object stores data used for NanoMethViz visualisation. It
#' contains stores a path to the methylation data, sample information and
#' optional exon information. The object is constructed using the
#' NanoMethResult() constructor function described in "Usage".
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

#' @describeIn NanoMethResult-class Constructor
#'
#' @param methy the path to methylation database (sql or indexed tabix).
#' @param samples the data.frame of sample annotation containg at least columns
#'   sample and group.
#' @param exons (optional) the data.frame of exon information containing at least columns
#'   gene_id, chr, strand, start, end, transcript_id and symbol.
#'
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

    assertthat::is.readable(methy)

    methods::new(
        "NanoMethResult",
        methy = methy,
        samples = tibble::as_tibble(samples),
        exons = tibble::as_tibble(exons)
    )
}


#' Get methylation data
#' @keywords internal
#'
#' @param object the object.
#'
#' @export
setGeneric("methy", valueClass = "character", function(object) {
    standardGeneric("methy")
})

#' @describeIn NanoMethResult-class methyation data path getter.
#'
#' @param object the NanoMethResult object.
#'
#' @export
setMethod("methy", signature("NanoMethResult"), function(object)
{
    object@methy
})

#' Get sample annotation
#'
#' @param object the object.
#'
#' @export
#'
#' @keywords internal
setGeneric("samples", valueClass = "data.frame", function(object) {
    standardGeneric("samples")
})

#' @describeIn NanoMethResult-class sample annotation getter.
#'
#' @param object the NanoMethResult object.
#'
#' @export
setMethod("samples", signature("NanoMethResult"), function(object)
{
    object@samples
})

#' Get exon annotation
#' @keywords internal
#'
#' @param object the object.
#'
#' @export
setGeneric("exons", valueClass = "data.frame", function(object) {
    standardGeneric("exons")
})

#' @describeIn NanoMethResult-class exon annotation getter.
#'
#' @param object the NanoMethResult object.
#'
#' @export
setMethod("exons", signature("NanoMethResult"), function(object)
{
    object@exons
})
