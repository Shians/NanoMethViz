#' @title Convert BAM with modifications to tabix format
#'
#' @description The `modbam_to_tabix` function takes a ModBamResult object and
#'   converts it into a tabix file format, which is efficient for indexing and
#'   querying large datasets.
#'
#' @param x the `ModBamResult` object.
#' @param out_file the path of the output tabix.
#'
#' @return invisibly returns the name of the created tabix file.
#'
#' @examples
#' out_file <- tempfile()
#' mbr <- ModBamResult(
#'     methy = ModBamFiles(
#'         samples = "sample1",
#'         paths = "inst/peg3.bam"
#'     ),
#'     samples = data.frame(
#'         sample = "sample1",
#'         group = "group1"
#'     )
#' )
#'
#' modbam_to_tabix(mbr, out_file)
#'
#' @export
modbam_to_tabix <- function(x, out_file) {
    assertthat::assert_that(is(x, "ModBamResult"))

    bam_info <- dplyr::inner_join(samples(x), methy(x), by = dplyr::join_by(sample))

    if (fs::file_exists(out_file)) {
        cli::cli_progress_step(paste0("Output file exists, overwriting ", out_file))
        fs::file_delete(out_file)
    }

    parse_and_select <- function(x) {
        parse_modbam(x[[1]], sample) %>%
            select("sample", "chr", "pos", "strand", "statistic", "read_name")
    }

    total <- get_bam_total_reads(path)
    fname <- fs::path_file(path)
    n_files <- nrow(bam_info)
    for (i in seq_len(n_files)) {
        path <- bam_info$path[i]
        sample <- bam_info$sample[i]
        bam_file <- Rsamtools::BamFile(path, yieldSize = 15000)

        cli::cli_progress_bar(
            glue::glue("Converting file {i}/{n_files}: {fname}"),
            total = total,
            format_done = paste0(
                "{.alert-success Data converted: ", fname, " {.timestamp {cli::pb_elapsed}}}"),
            format_failed = paste0(
                "{.alert-danger Data conversion failed: ", fname, " {.timestamp {cli::pb_elapsed}}}"),
            clear = FALSE
        )

        open(bam_file)
        while (Rsamtools::isIncomplete(bam_file)) {
            reads <- read_modbam(bam_file)
            data <- parse_and_select(reads)

            readr::write_tsv(data, out_file, append = TRUE, progress = FALSE)
            cli::cli_progress_update(length(reads[[1]][[1]]))
        }
        close(bam_file)
    }

    cli::cli_progress_step("Sorting data")

    f <- sort_methy_file(out_file)

    cli::cli_progress_step("Compressing data")
    out <- tabix_compress(f)
    fs::file_delete(f)

    cli::cli_progress_step(paste0("Tabix file created: ", out))
    invisible(out)
}
