#' Sort methylation file
#' @keywords internal
#'
#' @param x the path to the methylation file to sort
#'
#' @return invisibly returns path of sorted file
sort_methy_file <- function(x) {
    assert_that(is.readable(x))

    if (.Platform$OS.type == "windows") {
        methy_df <- readr::read_tsv(
            x,
            col_names = methy_col_names(),
            col_types = methy_col_types())
        methy_df <- dplyr::arrange(methy_df, .data$chr, .data$pos)
        readr::write_tsv(methy_df, x, col_names = FALSE)
    } else {
        cmd <- glue::glue("sort -k2,3V {x} -o {x}")
        system(cmd)
    }

    invisible(x)
}

tabix_compress <- function(x, index = TRUE) {
    assert_that(is.readable(x))

    f <- Rsamtools::bgzip(x, overwrite = TRUE)
    if (index) {
        tabix_index(f)
    }

    f
}

tabix_index <- function(x) {
    assert_that(is.readable(x))

    Rsamtools::indexTabix(x, seq = 2, start = 3, end = 3)
}

#' Convert methylation file to tabix format
#' @keywords internal
#'
#' @param x the path to the sorted methylation file
#'
#' @return invisibly returns the path to the tabix file
raw_methy_to_tabix <- function(x) {
    assert_that(is.readable(x))

    bgz_name <- tabix_compress(x)
    tabix_index(bgz_name)

    invisible(bgz_name)
}

#' Create a tabix file using methylation calls
#'
#' @param input_files the files to convert
#' @param output_file the output file to write results to (must end in .bgz)
#' @param samples the names of samples corresponding to each file
#' @param verbose TRUE if progress messages are to be printed
#'
#' @return invisibly returns the output file path, creates a tabix file (.bgz)
#'   and its index (.bgz.tbi)
#' @export
#'
#' @examples
#' methy_calls <- system.file(package = "NanoMethViz",
#'     c("sample1_nanopolish.tsv.gz", "sample2_nanopolish.tsv.gz"))
#' temp_file <- paste0(tempfile(), ".tsv.bgz")
#'
#' create_tabix_file(methy_calls, temp_file)
create_tabix_file <- function(
    input_files,
    output_file,
    samples = extract_file_names(input_files),
    verbose = TRUE
) {
    assert_that(
        tools::file_ext(output_file) == "bgz",
        "output_file must end with .bgz extension."
    )
    assert_that(
        assertthat::not_empty(input_files),
        is.character(input_files),
        is.string(output_file),
        is.character(samples),
        length(input_files) == length(samples)
    )

    if (.Platform$OS.type == "windows") {
        timed_log("WARNING: creating tabix file on windows requires at least twice as much memory as total size of methylation data")
    }

    if (packageVersion("cpp11") < package_version("0.2.5")) {
        warning("cpp11 versions < 0.2.5 may crash when reading in large tables")
    }

    temp_file <- tempfile()
    if (verbose) {
        timed_log("creating methylation table")
    }
    convert_methy_format(input_files, temp_file, samples, verbose = verbose)

    if (verbose) {
        timed_log("sorting methylation table")
    }
    sort_methy_file(temp_file)

    if (verbose) {
        timed_log("compressing methylation table to tabix with index")
    }
    raw_methy_to_tabix(temp_file)

    fs::file_move(paste0(temp_file, ".bgz"), output_file)
    fs::file_move(paste0(temp_file, ".bgz.tbi"), paste0(output_file, ".tbi"))
    fs::file_delete(temp_file)

    invisible(output_file)
}
