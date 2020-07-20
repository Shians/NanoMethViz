#' Create BSSeq object from methylation tabix file
#'
#' @param methy the path to the methylation tabix file.
#' @param out_folder the folder to store intermediate files. One file is created
#'   for each sample and contains columns "chr", "pos", "total" and
#'   "methylated".
#'
#' @return a BSSeq object.
#' @export
#'
#' @examples
#' methy <- system.file("methy_subset.tsv.bgz", package = "NanoMethViz")
#' bsseq <- methy_to_bsseq(methy)
methy_to_bsseq <- function(
    methy,
    out_folder = tempdir()
) {
    timed_log("creating intermediate files...")
    files <- convert_methy_to_dss(methy, out_folder)

    timed_log("creating bsseq object...")
    out <- create_bsseq_from_files(files$file_path, files$sample)
    timed_log("done")

    out
}

convert_methy_to_dss <- function(
    methy,
    out_folder
) {
    assert_that(
        fs::file_exists(methy)
    )

    if (!fs::dir_exists(out_folder)) {
        fs::dir_create(out_folder, recurse = TRUE)
    }

    samples <- convert_methy_to_dss_cpp(
        fs::path_expand(methy),
        fs::path_expand(out_folder)
    )

    data.frame(
        sample = samples,
        file_path = path(out_folder, paste0(samples, ".txt"))
    )
}

create_bsseq_from_files <- function(paths, samples) {
    read_dss <- purrr::partial(
        read_tsv,
        col_types = cols(
            chr = col_character(),
            pos = col_double(),
            total = col_double(),
            methylated = col_double()
        ),
        progress = FALSE
    )

    timed_log("reading in parsed data...")
    # read in data
    dat <- purrr::map(paths, read_dss)

    timed_log("constructing matrices...")
    # get unique positions
    combine_distinct_gpos <- function(x, y) {
        rbind(
            dplyr::select(x, chr, pos),
            dplyr::select(y, chr, pos)
        ) %>%
            dplyr::distinct()
    }

    unique_pos_df <- purrr::reduce(
            dat,
            combine_distinct_gpos) %>%
        dplyr::arrange(chr, pos) %>%
        dplyr::mutate(id = paste(chr, pos))

    # create methylation matrix
    M_mat <- matrix(
        0,
        nrow = nrow(unique_pos_df),
        ncol = length(dat),
        dimnames = list(NULL, samples))

    for (i in 1:length(dat)) {
        row_inds <- with(dat[[i]], paste(chr, pos)) %>%
            factor(levels = unique_pos_df$id)

        M_mat[row_inds, i] <- dat[[i]]$methylated
    }

    # create coverage matrix
    cov_mat <- matrix(
        0,
        nrow = nrow(unique_pos_df),
        ncol = length(dat),
        dimnames = list(NULL, samples))

    for (i in 1:length(dat)) {
        row_inds <- with(dat[[i]], paste(chr, pos)) %>%
            factor(levels = unique_pos_df$id)

        cov_mat[row_inds, i] <- dat[[i]]$total
    }

    # create BSseq object
    result <- bsseq::BSseq(
        chr = unique_pos_df$chr,
        pos = unique_pos_df$pos,
        M = M_mat,
        Cov = cov_mat,
        sampleNames = samples
    )

    SummarizedExperiment::colData(result) <- S4Vectors::DataFrame(
        sample = samples
    )

    result
}
