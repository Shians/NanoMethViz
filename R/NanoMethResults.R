#' Nanopore Methylation Result
#'
#' A NanoMethResult object stores data used for NanoMethViz visualisation. It
#' contains stores a path to the methylation data, sample information and
#' optional exon information. The object is constructed using the
#' NanoMethResult() constructor function described in "Usage".
#'
#' @slot methy the path to the methylation tabix file.
#' @slot samples the data.frame of sample annotation containing at least columns
#'   sample and group.
#' @slot exons the data.frame of exon information containing at least columns
#'   gene_id, chr, strand, start, end, transcript_id and symbol.
#'
#' @return a NanoMethResult object to be used with plotting functions
#'
#' @examples
#' methy <- system.file(package = "NanoMethViz", "methy_subset.tsv.bgz", mustWork = FALSE)
#' sample <- c(
#'     "B6Cast_Prom_1_bl6",
#'     "B6Cast_Prom_1_cast",
#'     "B6Cast_Prom_2_bl6",
#'     "B6Cast_Prom_2_cast",
#'     "B6Cast_Prom_3_bl6",
#'     "B6Cast_Prom_3_cast"
#' )
#' group <- c(
#'     "bl6",
#'     "cast",
#'     "bl6",
#'     "cast",
#'     "bl6",
#'     "cast"
#' )
#' sample_anno <- data.frame(sample, group, stringsAsFactors = FALSE)
#' exon_tibble <- get_example_exons_mus_musculus()
#' NanoMethResult(methy, sample_anno, exon_tibble)
#'
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
#' @param methy the path to the methylation tabix file.
#' @param samples the data.frame of sample annotation containing at least
#'   columns sample and group.
#' @param exons (optional) the data.frame of exon information containing at
#'   least columns gene_id, chr, strand, start, end, transcript_id and symbol.
#'
#' @export
NanoMethResult <- function(methy, samples, exons = NULL) {
    # validate input
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

    # validate inputs
    assert_readable(methy)
    assert_valid_samples(samples, context = "NanoMethResult constructor")
    assert_valid_exons(exons)

    # convert group to factor
    samples$group <- as.factor(samples$group)

    # ensure samples match
    assert_valid_methy_samples(methy, samples)

    methods::new(
        "NanoMethResult",
        methy = methy,
        samples = tibble::as_tibble(samples),
        exons = tibble::as_tibble(exons)
    )
}

# methy ----

#' @describeIn NanoMethResult-class methylation data path getter.
#'
#' @param object the NanoMethResult object.
#'
#' @return the path to the methylation data.
#'
#' @examples
#' x <- load_example_nanomethresult()
#' methy(x)
#'
#' @export
setMethod("methy", signature("NanoMethResult"), function(object) {
    object@methy
})

#' @describeIn NanoMethResult-class methylation data path setter.
#'
#' @param object the NanoMethResult object.
#' @param value the path to the methylation data.
#'
#' @export
setMethod("methy<-", signature("NanoMethResult"), function(object, value) {
    assert_readable(value)

    object@methy <- value
    object
})

# sample ----
#' @describeIn NanoMethResult-class sample annotation getter.
#'
#' @param object the NanoMethResult object.
#'
#' @return the sample annotation.
#'
#' @export
setMethod("samples", signature("NanoMethResult"), function(object) {
    object@samples
})

#' @describeIn NanoMethResult-class sample annotation setter.
#'
#' @param object the NanoMethResult object.
#' @param value the data.frame of sample annotation containing at least columns
#'   sample and group.
#'
#' @export
setMethod("samples<-", signature("NanoMethResult", "data.frame"), function(object, value) {
    assert_valid_samples(value, context = "NanoMethResult samples setter")

    # Convert group to factor
    value$group <- as.factor(value$group)

    object@samples <- tibble::as_tibble(value)
    object
})

# exons ----
#' @describeIn NanoMethResult-class exon annotation getter.
#'
#' @param object the NanoMethResult object.
#'
#' @return the exon annotation.
#'
#' @export
setMethod("exons", signature("NanoMethResult"), function(object) {
    object@exons
})

#' @describeIn NanoMethResult-class exon annotation setter.
#'
#' @param object the NanoMethResult object.
#' @param value the exon annotation.
#'
#' @export
setMethod("exons<-", signature("NanoMethResult", "data.frame"), function(object, value) {
    assert_valid_exons(value)

    object@exons <- tibble::as_tibble(value)
    object
})
