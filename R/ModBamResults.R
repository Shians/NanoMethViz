
setClass(
    "ModBamFiles",
    contains = "data.frame"
)

ModBamFiles <- function(samples, paths) {
    assert_readable(paths)

    x <- data.frame(
        sample = samples,
        path = paths
    )

    new("ModBamFiles", x)
}

setMethod("show", signature("ModBamFiles"), function(object) {
    print(glue::glue("A ModBamFiles object containing {nrow(object)} samples:"))
    print(object)
})

#' Modbam methylation results
#'
#' (Experimental) A ModBamResult object stores data used for NanoMethViz
#' #'visualisation. It contains stores a path to the methylation data, sample
#' information and optional exon information. The object is constructed using
#' the NanoMethResult() constructor function described in "Usage".
#'
#' @slot methy a ModBamFiles data.frame specifying the samples and paths to bam
#'   files.
#' @slot samples the data.frame of sample annotation containing at least columns
#'   sample and group.
#' @slot exons the data.frame of exon information containing at least columns
#'   gene_id, chr, strand, start, end, transcript_id and symbol.
#'
#' @return a NanoMethResult object to be used with plotting functions
#'
#' @export
setClass(
    "ModBamResult",
    representation(
        methy = "ModBamFiles",
        samples = "data.frame",
        exons = "data.frame"
    )
)

setMethod(
    "methy",
    signature(object = "ModBamResult"),
    definition = function(object) {
        object@methy
    }
)

setMethod(
    "methy<-",
    signature(object = "ModBamResult", value = "ModBamFiles"),
    definition = function(object, value) {
        object@methy <- value
        object
    }
)

setMethod(
    "samples",
    signature(object = "ModBamResult"),
    definition = function(object) {
        object@samples
    }
)

setMethod(
    "samples<-",
    signature(object = "ModBamResult", value = "data.frame"),
    definition = function(object, value) {
        object@samples <- value
        object
    }
)

setMethod(
    "exons",
    signature(object = "ModBamResult"),
    definition = function(object) {
        object@exons
    }
)

setMethod(
    "exons<-",
    signature(object = "ModBamResult", value = "data.frame"),
    definition = function(object, value) {
        object@exons <- value
        object
    }
)

ModBamResult <- function(methy, samples, exons = NULL) {
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

    assert_has_columns(
        exons,
        c("gene_id", "chr", "strand", "start", "end", "transcript_id", "symbol")
    )
    assert_has_columns(samples, c("sample", "group"))
    samples$group <- as.factor(samples$group)

    methods::new(
        "ModBamResult",
        methy = methy,
        samples = tibble::as_tibble(samples),
        exons = tibble::as_tibble(exons)
    )
}
