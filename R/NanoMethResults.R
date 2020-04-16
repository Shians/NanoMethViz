#' Nanopore Methylation Result
#'
#' @slot methy the path to methylation database (sql or indexed tabix).
#' @slot samples the data.frame of sample annotation containg at least columns
#'   sample and group.
#' @slot exons the data.frame of exon information containing at least columns
#'   gene_id, chr, strand, start, end, transcript_id and symbol.
#'
#' @seealso [NanoMethResult()], [methy()], [samples()], [exons()].
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

    methods::new(
        "NanoMethResult",
        methy = methy,
        samples = tibble::as_tibble(samples),
        exons = tibble::as_tibble(exons)
    )
}


#' Get methylation data
#'
#' @param object the object.
#'
#' @return path to methylation data.
#' @export
setGeneric("methy", valueClass = "character", function(object) {
    standardGeneric("methy")
})

#' Get methylation data
#'
#' @param object the NanoMethResult object.
#'
#' @return path to methylation data.
#' @describeIn NanoMethResult-class methyation data path getter.
#' @export
setMethod("methy", signature("NanoMethResult"), function(object)
{
    object@methy
})

#' Get sample annotation
#'
#' @param object the object.
#'
#' @return data.frame of sample annotation.
#' @export
setGeneric("samples", valueClass = "data.frame", function(object) {
    standardGeneric("samples")
})

#' Get sample annotation
#'
#' @param object the NanoMethResult object.
#'
#' @return data.frame of sample annotation.
#' @describeIn NanoMethResult-class sample annotation getter.
#' @export
setMethod("samples", signature("NanoMethResult"), function(object)
{
    object@samples
})

#' Get exon annotation
#'
#' @param object the object.
#'
#' @return data.frame of exon annotation.
#' @export
setGeneric("exons", valueClass = "data.frame", function(object) {
    standardGeneric("exons")
})

#' Get exon annotation
#'
#' @param object the NanoMethResult object.
#'
#' @return data.frame of exon annotation.
#' @describeIn NanoMethResult-class exon annotation getter.
#' @export
setMethod("exons", signature("NanoMethResult"), function(object)
{
    object@exons
})