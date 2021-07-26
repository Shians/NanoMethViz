#' Convert BSseq object to edgeR methylation matrix
#'
#' @param bsseq the BSseq object.
#'
#' @return a matrix compatible with the edgeR differential methylation pipeline
#' @export
#'
#' @examples
#' methy <- system.file("methy_subset.tsv.bgz", package = "NanoMethViz")
#' bsseq <- methy_to_bsseq(methy)
#' edger_mat <- bsseq_to_edger(bsseq)
bsseq_to_edger <- function(bsseq) {
    edger_col_names <- .get_edger_col_names(bsseq)
    edger_row_names <- .get_edger_row_names(bsseq)

    # construct matrix
    methylated <- .get_me_mat(bsseq)
    unmethylated <- .get_un_mat(bsseq)

    edger_mat <- matrix(
        0,
        ncol = 2*ncol(methylated),
        nrow = nrow(methylated),
        dimnames = list(edger_row_names, edger_col_names)
    )

    for (i in 0:(ncol(methylated) - 1)) {
        edger_mat[, 2*i + 1] <- methylated[, i + 1]
        edger_mat[, 2*i + 2] <- unmethylated[, i + 1]
    }

    edger_mat
}

#' Convert BSseq object to log-methylation-ratio matrix
#'
#' Creates a log-methylation-ratio matrix from a BSseq object that is useful for
#' dimensionality reduction plots.
#'
#' @param bsseq the BSseq object.
#' @param regions the regions to calculate log-methylation ratios over. If left NULL, ratios will be calculated per site.
#' @param prior_count the prior count added to avoid taking log of 0.
#'
#' @return a matrix containing log-methylation-ratios.
#' @export
#'
#' @examples
#' nmr <- load_example_nanomethresult()
#' bsseq <- methy_to_bsseq(nmr)
#' regions <- exons_to_genes(NanoMethViz::exons(nmr))
#' log_m_ratio <- bsseq_to_log_methy_ratio(bsseq, regions)
bsseq_to_log_methy_ratio <- function(bsseq, regions = NULL, prior_count = 2) {
    if (prior_count < 1) {
        warning("prior_count of 1 or higher is recommended")
    }

    if (!is.null(regions)) {
        assertthat::assert_that(is(regions, "data.frame"))
        regions <- GenomicRanges::GRanges(regions)
        row_names <- regions$gene_id
    } else {
        row_names <- .get_edger_row_names(bsseq)
    }

    col_names <- SummarizedExperiment::colData(bsseq)$sample
    methylated <- .get_me_mat(bsseq, regions)
    unmethylated <- .get_un_mat(bsseq, regions)

    log_mat <- log2(methylated + prior_count) - log2(unmethylated + prior_count)

    dimnames(log_mat) <- list(row_names, col_names)

    log_mat
}

.get_me_mat <- function(bsseq, regions) {
    if (!is.null(regions)) {
        bsseq::getCoverage(bsseq, regions, type = "M", what = "perRegionTotal")
    } else {
        bsseq::getCoverage(bsseq, type = "M")
    }
}

.get_un_mat <- function(bsseq, regions) {
    if (!is.null(regions)) {
        cov <- bsseq::getCoverage(bsseq, regions, type = "Cov", what = "perRegionTotal")
    } else {
        cov <- bsseq::getCoverage(bsseq, type = "Cov")
    }
    me_mat <- .get_me_mat(bsseq, regions)
    cov - me_mat
}

.get_edger_col_names <- function(bsseq) {
    samples <- SummarizedExperiment::colData(bsseq)$sample
    me_names <- paste0(samples, "_Me")
    un_names <- paste0(samples, "_Un")
    edger_col_names <- character(2 * length(samples))
    for (i in 0:(length(samples) - 1)) {
        edger_col_names[2*i + 1] <- me_names[i + 1]
        edger_col_names[2*i + 2] <- un_names[i + 1]
    }

    edger_col_names
}

.get_edger_row_names <- function(bsseq) {
    gr <- bsseq::getBSseq(bsseq, type = "gr")
    seq <- as.character(SummarizedExperiment::seqnames(gr))
    pos <- as.integer(SummarizedExperiment::start(gr))
    edger_row_names <- paste(seq, pos, sep = "-")

    edger_row_names
}
