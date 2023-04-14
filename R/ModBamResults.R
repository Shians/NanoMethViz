#' ModBamFiles class
#'
#' This is a class for holding information about modbam files. It is a
#' data.frame containing information about samples and paths to modbam files.
#'
#' @export
setClass(
    "ModBamFiles",
    contains = "data.frame"
)

#' Constructor for a ModBamFiles object
#'
#' This function creates a ModBamFiles object containing information about the
#' samples and file paths.
#'
#' @param samples a character vector with the names of the samples.
#' @param paths a character vector with the file paths for the BAM files.
#'
#' @return A ModBamFiles object with the sample and path information.
#'
#' @export
ModBamFiles <- function(samples, paths) {
    assert_readable(paths)

    x <- data.frame(
        sample = samples,
        path = paths
    )

    new("ModBamFiles", x)
}

#' @docType methods
#' @rdname ModBamFiles
#'
#' @param object a ModBamFiles object.
#' @export
setMethod("show", signature("ModBamFiles"), function(object) {
    print(glue::glue(
        "A ModBamFiles object containing {sample_n} samples:",
        sample_n = nrow(object)
    ))
    print(object)
})

#' Modbam methylation results
#'
#' (Experimental) A ModBamResult object stores modbam data used for NanoMethViz
#' visualisation.
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

#' @describeIn ModBamResult-class modbam information getter.
#'
#' @param object the ModBamResult object.
#'
#' @return a ModBamFiles data.frame.
#'
#' @export
setMethod(
    "methy",
    signature(object = "ModBamResult"),
    definition = function(object) {
        object@methy
    }
)

#' @describeIn ModBamResult-class modbam information setter.
#'
#' @param object the ModBamResult object.
#' @param value the path to the methylation data.
#'
#' @export
setMethod(
    "methy<-",
    signature(object = "ModBamResult", value = "ModBamFiles"),
    definition = function(object, value) {
        object@methy <- value
        object
    }
)

#' @describeIn ModBamResult-class sample annotation getter.
#'
#' @param object the ModBamResult object.
#'
#' @return the sample annotation.
#'
#' @export
setMethod(
    "samples",
    signature(object = "ModBamResult"),
    definition = function(object) {
        object@samples
    }
)

#' @describeIn ModBamResult-class sample annotation setter.
#'
#' @param object the ModBamResult object.
#' @param value the data.frame of sample annotation containing at least columns
#'   sample and group.
#'
#' @export
setMethod(
    "samples<-",
    signature(object = "ModBamResult", value = "data.frame"),
    definition = function(object, value) {
        object@samples <- value
        object
    }
)

#' @describeIn ModBamResult-class exon annotation getter.
#'
#' @param object the ModBamResult object.
#'
#' @return the exon annotation.
#'
#' @export
setMethod(
    "exons",
    signature(object = "ModBamResult"),
    definition = function(object) {
        object@exons
    }
)

#' @describeIn ModBamResult-class exon annotation setter.
#'
#' @param object the ModBamResult object.
#' @param value the exon annotation.
#'
#' @return the exon annotation.
#'
#' @export
setMethod(
    "exons<-",
    signature(object = "ModBamResult", value = "data.frame"),
    definition = function(object, value) {
        object@exons <- value
        object
    }
)

#' @describeIn ModBamResult-class Constructor
#'
#' @param methy a ModBamFiles object.
#' @param samples the data.frame of sample annotation containing at least
#'   columns sample and group.
#' @param exons (optional) the data.frame of exon information containing at
#'   least columns gene_id, chr, strand, start, end, transcript_id and symbol.
#'
#' @export
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
