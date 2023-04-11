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
