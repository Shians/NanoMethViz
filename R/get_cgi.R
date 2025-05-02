#' Get CpG islands annotations
#'
#' Helper functions are provided for obtaining CpG islands annotations from UCSC table browser.
#'
#' @keywords internal
#'
#' @return tibble (data.frame) object containing CpG islands annotation.
#'
get_cgi <- function(genome) {
    available_genomes <- c("hg19", "hg38", "t2t", "mm10", "grcm39")
    assertthat::assert_that(
        genome %in% available_genomes,
        msg = sprintf("genome must be one of %s", paste(available_genomes, collapse = ", "))
    )

    rds_path <- switch(
        genome,
        hg19 = system.file("cgi_hg19.rds", package = "NanoMethViz", mustWork = FALSE),
        hg38 = system.file("cgi_hg38.rds", package = "NanoMethViz", mustWork = FALSE),
        t2t = system.file("cgi_t2t.rds", package = "NanoMethViz", mustWork = FALSE),
        mm10 = system.file("cgi_mm10.rds", package = "NanoMethViz", mustWork = FALSE),
        grcm39 = system.file("cgi_GRCm39.rds", package = "NanoMethViz", mustWork = FALSE)
    )

    readRDS(rds_path)
}

#' @rdname get_exons
#'
#' @examples
#' cgi_mm10 <- get_cgi_mm10()
#'
#' @export
get_cgi_mm10 <- function() {
    get_cgi("mm10")
}

#' @rdname get_exons
#'
#' @examples
#' cgi_GRCm39 <- get_cgi_grcm39()
#'
#' @export
get_cgi_grcm39 <- function() {
    get_cgi("grcm39")
}

#' @rdname get_exons
#'
#' @examples
#' cgi_t2t <- get_cgi_t2t()
#'
#'# @export
get_cgi_t2t <- function() {
    get_cgi("t2t")
}

#' @rdname get_exons
#'
#' @examples
#' cgi_hg19 <- get_cgi_hg19()
#'
#' @export
get_cgi_hg19 <- function() {
    get_cgi("hg19")
}

#' @rdname get_exons
#'
#' @examples
#' cgi_hg38 <- get_cgi_hg38()
#'
#' @export
get_cgi_hg38 <- function() {
    get_cgi("hg38")
}
