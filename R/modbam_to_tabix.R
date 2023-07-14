#' @title Convert BAM with modifications to tabix format
#'
#' @description The `modbam_to_tabix` function takes a ModBamResult object and
#'   converts it into a tabix file format, which is efficient for indexing and
#'   querying large datasets.
#'
#' @param x the `ModBamResult` object.
#' @param out_file the path of the output tabix.
#' @param mod_code the modification code to use, defaults to 'm' for 5mC
#'  methylation.
#'
#' @details
#' The possible tags for mod_code can be found at
#'  \url{https://samtools.github.io/hts-specs/SAMtags.pdf} under the
#'  'Base modifications' section.
#'
#' @return invisibly returns the name of the created tabix file.
#'
#' @examples
#' out_file <- tempfile()
#' mbr <- ModBamResult(
#'     methy = ModBamFiles(
#'         samples = "sample1",
#'         paths = system.file("peg3.bam", package = "NanoMethViz")
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
modbam_to_tabix <- function(x, out_file, mod_code = NanoMethViz::mod_code(x)) {
    assertthat::assert_that(is(x, "ModBamResult"))

    bam_info <- dplyr::inner_join(samples(x), methy(x), by = dplyr::join_by(sample))

    if (fs::file_exists(out_file)) {
        cli::cli_progress_step(paste0("Output file exists, overwriting ", out_file))
        fs::file_delete(out_file)
    }

    parse_read_chunk <- function(x) {
        parse_modbam(x[[1]], sample, mod_code = mod_code) %>%
            select("sample", "chr", "pos", "strand", "statistic", "read_name")
    }

    n_files <- nrow(bam_info)
    for (i in seq_len(n_files)) {
        path <- bam_info$path[i]
        sample <- bam_info$sample[i]
        bam_file <- Rsamtools::BamFile(path, yieldSize = 15000)

        total <- get_bam_total_reads(path)
        fname <- fs::path_file(path)
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
            reads <- read_bam(bam_file)
            if (!is.null(reads[[1]])) {
                # parse if valid data exists
                data <- parse_read_chunk(reads)
                readr::write_tsv(data, out_file, append = TRUE, progress = FALSE)
                cli::cli_progress_update(length(reads[[1]][[1]]))
            }
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
