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
#'         paths = system.file("peg3.bam", package = "NanoMethViz",
#'         mustWork = FALSE)
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
    assert_that(is(x, "ModBamResult"))

    assert_that(
        tools::file_ext(out_file) == "bgz",
        msg = "output_file must end with .bgz extension."
    )

    assert_that(
        !(fs::file_exists(out_file) && !fs::is_file(out_file)),
        msg = "output_file exists and is not a file"
    )

    if (fs::file_exists(out_file)) {
        cli::cli_alert_info(paste0("Output file exists, will overwrite ", out_file))
    }

    tmp_tsv_path <- base::tempfile(fileext = ".tsv")
    cli::cli_progress_step(paste0("Writing data to temporary file: ", tmp_tsv_path))

    cli::cli_progress_step("Converting data to TSV")
    tsv_file <- run_modbam_to_tsv_converter(x, tmp_tsv_path, mod_code)

    cli::cli_progress_step("Sorting data")
    sorted <- sort_methy_file(tsv_file)

    cli::cli_progress_step("Compressing data")
    tmp_out <- tabix_compress(sorted)

    cli::cli_alert_info(paste0("Moving data to final location: ", out_file))

    output_dir <- fs::path_dir(out_file)
    if (output_dir != "" && fs::file_exists(output_dir)) {
        fs::dir_create(fs::path_dir(out_file))
    }

    fs::file_move(tmp_out, out_file)

    cli::cli_progress_step(paste0("Tabix file created: ", out_file))
    invisible(out_file)
}

run_modbam_to_tsv_converter <- function(x, out_file, mod_code) {
    # if .bgz at end of output name then trim it so final output
    # doesn't end with .bgz.bgz
    if (stringr::str_detect(out_file, ".bgz$")) {
        out_file <- out_file %>%
            stringr::str_remove(".bgz$")
    }

    bam_info <- dplyr::inner_join(samples(x), methy(x), by = dplyr::join_by(sample))

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
        modbam_file_to_tsv(path, out_file, sample, mod_code, prog_bar_id)
    }

    out_file
}

modbam_file_to_tsv <- function(path, out_file, sample, mod_code, prog_bar_id = NULL) {
    parse_read_chunk <- function(x) {
        parse_modbam(x[[1]], sample, mod_code = mod_code) %>%
            select("sample", "chr", "pos", "strand", "statistic", "read_name")
    }

    bam_file <- Rsamtools::BamFile(path, yieldSize = 2000)
    open(bam_file)
    while (Rsamtools::isIncomplete(bam_file)) {
        reads <- read_bam(bam_file)

        if (!is.null(reads[[1]]) && length(reads[[1]]$qname) > 0) {
            # parse if valid data exists
            data <- parse_read_chunk(reads)
            readr::write_tsv(data, out_file, append = TRUE, progress = FALSE)
        }

        if (!is.null(prog_bar_id)) {
            n_reads <- length(reads[[1]][[1]])
            cli::cli_progress_update(n_reads, id = prog_bar_id)
        }
    }
    close(bam_file)
}
