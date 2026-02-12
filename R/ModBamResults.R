#' ModBamFiles class
#'
#' This is a class for holding information about modBAM files. It is a
#' data.frame containing information about samples and paths to modBAM files.
#'
#' @export
setClass(
    "ModBamFiles",
    contains = "data.frame"
)

#' Constructor for a ModBamFiles object
#'
#' This function creates a ModBamFiles object containing information about the
#' samples and file paths. This constructor checks that the files are readable
#' and have an index.
#'
#' @param samples a character vector with the names of the samples.
#' @param paths a character vector with the file paths for the BAM files.
#'
#' @return A ModBamFiles object with the sample and path information.
#'
#' @importFrom tibble tibble
#' @export
ModBamFiles <- function(samples, paths) {
    # validate inputs
    if (length(samples) != length(paths)) {
        stop(glue::glue(
            "Length of samples ({length(samples)}) must equal length of paths ({length(paths)}).\n",
            "Please provide matching samples and paths."
        ))
    }

    if (length(samples) == 0) {
        stop("At least one sample and path must be provided.")
    }

    if (any(is.na(samples) | samples == "")) {
        stop("Sample names cannot be empty or NA. Please provide valid sample identifiers.")
    }

    if (any(duplicated(samples))) {
        duplicates <- samples[duplicated(samples)]
        stop(glue::glue(
            "Found duplicate sample names: {paste(unique(duplicates), collapse = ', ')}\n",
            "Each sample must have a unique identifier."
        ))
    }

    assert_readable(paths)
    assert_has_index(paths)

    x <- tibble::tibble(
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

#' modBAM methylation results
#'
#' A ModBamResult object stores modBAM data used for NanoMethViz
#' visualisation. It contains stores a ModBamFiles object, sample information
#' and optional exon information. The object is constructed using the
#' ModBamResult() constructor function described in "Usage".
#'
#' @slot methy a ModBamFiles data.frame specifying the samples and paths to bam
#'   files.
#' @slot samples the data.frame of sample annotation containing at least columns
#'   sample and group.
#' @slot exons the data.frame of exon information containing at least columns
#'   gene_id, chr, strand, start, end, transcript_id and symbol.
#' @slot mod_code the modification code of interest.
#'
#' @return a ModBamResult object to be used with plotting functions
#'
#' @export
setClass(
    "ModBamResult",
    representation(
        methy = "ModBamFiles",
        samples = "data.frame",
        exons = "data.frame",
        mod_code = "character"
    )
)

#' @describeIn ModBamResult-class modBAM information getter.
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

#' @describeIn ModBamResult-class modBAM information setter.
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
        assert_valid_samples(value, context = "ModBamResult samples setter")

        # Convert group to factor
        value$group <- as.factor(value$group)

        object@samples <- tibble::as_tibble(value)
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
#' @export
setMethod(
    "exons<-",
    signature(object = "ModBamResult", value = "data.frame"),
    definition = function(object, value) {
        assert_valid_exons(value)

        object@exons <- tibble::as_tibble(value)
        object
    }
)

#' @describeIn ModBamResult-class mod code getter.
#'
#' @param object the ModBamResult object.
#'
#' @return the mod code.
#'
#' @export
setMethod(
    "mod_code",
    signature(object = "ModBamResult"),
    definition = function(object) {
        object@mod_code
    }
)

#' @describeIn ModBamResult-class mod code setter.
#'
#' @param object the ModBamResult object.
#' @param value the mod code.
#'
#' @export
setMethod(
    "mod_code<-",
    signature(object = "ModBamResult", value = "character"),
    definition = function(object, value) {
        # validate mod_code
        if (is.na(value) || value == "") {
            stop("Modification code cannot be empty or NA. Common codes are 'm' for 5mC, 'h' for 5hmC.")
        }

        value <- as.character(value)
        if (nchar(value) != 1) {
            stop(glue::glue(
                "Modification code must be a single character. Got: '{value}'\n",
                "Common codes are 'm' for 5mC, 'h' for 5hmC, 'a' for 6mA."
            ))
        }

        object@mod_code <- value
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
#' @param mod_code a character with the mod code of interest. Defaults to "m"
#'  for 5mC. See details for other options.
#'
#' @details
#' The possible tags for mod_code can be found at
#' \url{https://samtools.github.io/hts-specs/SAMtags.pdf} under the
#' 'Base modifications' section.
#'
#' @export
ModBamResult <- function(methy, samples, exons = NULL, mod_code = "m") {
    # validate inputs
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

    # validate ModBamFiles object
    if (!is(methy, "ModBamFiles")) {
        stop(glue::glue(
            "The 'methy' argument must be a ModBamFiles object.\n",
            "Got object of class: {class(methy)[1]}\n",
            "Create a ModBamFiles object using: ModBamFiles(samples = ..., paths = ...)"
        ))
    }

    # validate samples and exons
    assert_valid_samples(samples, context = "ModBamResult constructor")
    assert_valid_exons(exons)

    # validate mod_code
    if (is.na(mod_code) || mod_code == "") {
        stop("Modification code cannot be empty or NA. Common codes are 'm' for 5mC, 'h' for 5hmC.")
    }

    mod_code <- as.character(mod_code)
    if (nchar(mod_code) != 1) {
        stop(glue::glue(
            "Modification code must be a single character. Got: '{mod_code}'\n",
            "Common codes are 'm' for 5mC, 'h' for 5hmC, 'a' for 6mA."
        ))
    }

    # Check that samples in annotation match ModBamFiles samples
    modbam_samples <- methy$sample
    anno_samples <- samples$sample

    matched_samples <- intersect(modbam_samples, anno_samples)
    if (length(matched_samples) == 0) {
        stop(glue::glue(
            "No sample names match between ModBamFiles and sample annotation.\n",
            "ModBamFiles samples: {paste(modbam_samples, collapse = ', ')}\n",
            "Annotation samples: {paste(anno_samples, collapse = ', ')}\n",
            "Please ensure sample names match between your BAM files and annotation."
        ))
    }

    # warn about unmatched samples
    unmatched_bam <- setdiff(modbam_samples, anno_samples)
    unmatched_anno <- setdiff(anno_samples, modbam_samples)

    if (length(unmatched_bam) > 0) {
        warning(glue::glue(
            "Found {length(unmatched_bam)} samples in ModBamFiles not in annotation: {paste(unmatched_bam, collapse = ', ')}\n",
            "These samples will be ignored in analysis."
        ))
    }

    if (length(unmatched_anno) > 0) {
        warning(glue::glue(
            "Found {length(unmatched_anno)} samples in annotation not in ModBamFiles: {paste(unmatched_anno, collapse = ', ')}\n",
            "These samples will have no data for analysis."
        ))
    }

    # convert group to factor
    samples$group <- as.factor(samples$group)

    message(glue::glue(
        "Successfully created ModBamResult with {length(matched_samples)} matched samples."
    ))

    methods::new(
        "ModBamResult",
        methy = methy,
        samples = tibble::as_tibble(samples),
        exons = tibble::as_tibble(exons),
        mod_code = mod_code
    )
}
