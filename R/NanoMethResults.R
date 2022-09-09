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
#' methy <- system.file(package = "NanoMethViz", "methy_subset.tsv.bgz")
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
    assert_has_columns(
        exons,
        c("gene_id", "chr", "strand", "start", "end", "transcript_id", "symbol")
    )
    assert_has_columns(samples, c("sample", "group"))
    samples$group <- as.factor(samples$group)

    # Check in first 1000 entries that samples are in sample annotation
    head_values <- read.table(
        gzfile(methy),
        col.names = methy_col_names(),
        nrows = 1000
    )
    if (length(intersect(head_values$sample, samples$sample)) == 0) {
        stop("in first 1000 entries, no sample names matched samples from annotation")
    }
    if (!all(head_values$sample %in% samples$sample)) {
        warning(glue::glue(
            "in first 1000 entires, the following samples were not in annotation: {missing_cols}",
            missing_cols = paste(setdiff(head_values$sample, samples$sample), collapse = ", ")
        ))
    }

    assertthat::assert_that(any(head_values$sample %in% samples$sample))

    methods::new(
        "NanoMethResult",
        methy = methy,
        samples = tibble::as_tibble(samples),
        exons = tibble::as_tibble(exons)
    )
}

# methy ----

#' Get methylation data
#' @keywords internal
#'
#' @param object the object.
#'
#' @return the path to the methylation data.
#'
#' @examples
#' showMethods("methy")
#'
#' @export
setGeneric("methy", valueClass = "character", function(object) {
    standardGeneric("methy")
})

#' Set methylation data
#' @keywords internal
#' @export
setGeneric("methy<-", function(object, value) {
    standardGeneric("methy<-")
})

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
    assertthat::is.readable(value)

    object@methy <- value
    object
})

# sample ----

#' Get sample annotation
#'
#' @param object the object.
#'
#' @return the sample annotation.
#'
#' @examples
#' showMethods("samples")
#'
#' @export
#'
#' @keywords internal
setGeneric("samples", valueClass = "data.frame", function(object) {
    standardGeneric("samples")
})

#' Set sample annotation
#' @keywords internal
#' @export
setGeneric("samples<-", function(object, value) {
    standardGeneric("samples<-")
})

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
#' @param value the data.frame of sample annotation containg at least columns
#'   sample and group.
#'
#' @export
setMethod("samples<-", signature("NanoMethResult", "data.frame"), function(object, value) {
    assert_has_columns(value, c("sample", "group"))

    object@samples <- value
    object
})

# exons ----

#' Get exon annotation
#' @keywords internal
#'
#' @param object the object.
#'
#' @return the exon annotation.
#'
#' @examples
#' showMethods("exons")
#'
#' @export
setGeneric("exons", valueClass = "data.frame", function(object) {
    standardGeneric("exons")
})

#' Set exon annotation
#' @keywords internal
#' @export
setGeneric("exons<-", function(object, value) {
    standardGeneric("exons<-")
})

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

#' @describeIn NanoMethResult-class exon annotation getter.
#'
#' @param object the NanoMethResult object.
#' @param value the exon annotation.
#'
#' @export
setMethod("exons<-", signature("NanoMethResult", "data.frame"), function(object, value) {
    assert_has_columns(
        value,
        c("gene_id", "chr", "strand", "start", "end", "transcript_id", "symbol")
    )
    object@exons <- value
    object
})
