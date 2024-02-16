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
#' out_file <- paste0(tempfile(), ".tsv.bgz")
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

    if (fs::file_exists(out_file)) {
        cli::cli_progress_step(paste0("Output file exists, overwriting ", out_file))
        fs::file_delete(out_file)
    }

    cli::cli_progress_step("Converting data to TSV")
    tsv_file <- convert_modbam_to_tsv(x, out_file, mod_code)

    cli::cli_progress_step("Sorting data")
    f <- sort_methy_file(tsv_file)

    cli::cli_progress_step("Compressing data")
    out <- tabix_compress(f)
    fs::file_delete(f)

    cli::cli_progress_step(paste0("Tabix file created: ", out))
    invisible(out)
}

convert_modbam_to_tsv <- function(x, out_file, mod_code) {
    # if .gz at end of output name then trim it so final output
    # doesn't end with .bgz.bgz
    if (stringr::str_detect(out_file, ".bgz$")) {
        out_file <- out_file %>%
            stringr::str_remove(".bgz$")
    }

    bam_info <- dplyr::inner_join(samples(x), methy(x), by = dplyr::join_by(sample))

    parse_read_chunk <- function(x) {
        parse_modbam(x[[1]], sample, mod_code = mod_code) %>%
            select("sample", "chr", "pos", "strand", "statistic", "read_name")
    }

    n_files <- nrow(bam_info)
    for (i in seq_len(n_files)) {
        path <- bam_info$path[i]
        sample <- bam_info$sample[i]

        total <- get_bam_total_reads(path)
        fname <- fs::path_file(path)

        prog_bar_id <- cli::cli_progress_bar(
            glue::glue("Converting file {i}/{n_files}: {fname}"),
            total = total,
            format_done = paste0("{.alert-success Data converted: ", fname, " {.timestamp {cli::pb_elapsed}}}"),
            format_failed = paste0("{.alert-danger Data conversion failed: ", fname, " {.timestamp {cli::pb_elapsed}}}"),
            clear = FALSE
        )

        bam_file <- Rsamtools::BamFile(path, yieldSize = 15000)
        open(bam_file)
        while (Rsamtools::isIncomplete(bam_file)) {
            reads <- read_bam(bam_file)

            n_reads <- length(reads[[1]][[1]])
            cli::cli_progress_update(n_reads, id = prog_bar_id)

            if (!is.null(reads[[1]]) && length(reads[[1]]$qname) > 0) {
                # parse if valid data exists
                data <- parse_read_chunk(reads)
                readr::write_tsv(data, out_file, append = TRUE, progress = FALSE)
            }
        }

    }

    close(bam_file)

    out_file
}
