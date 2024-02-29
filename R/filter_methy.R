#' Create filtered methylation file
#'
#' Create a filtered methylation file from an existing one.
#'
#' @param x the path to the methylation file or a NanoMethResult object.
#' @param output_file the output file to write results to (must end in .bgz).
#' @param ... filtering criteria given in dplyr syntax. Use methy_col_names()
#'   to get available column names.
#'
#' @return
#' invisibly returns 'output_file' if x is a file path, otherwise returns
#' NanoMethResult object with methy(x) replaced with filtered value.
#'
#' @export
#'
#' @examples
#' nmr <- load_example_nanomethresult()
#' output_file <- paste0(tempfile(), ".tsv.bgz")
#' filter_methy(nmr, output_file = output_file, chr == "chrX")
#' filter_methy(methy(nmr), output_file = output_file, chr == "chrX")
filter_methy <- function(x, output_file, ...) {
    if (is(x, "NanoMethResult")) {
        return_nmr <- TRUE
        nmr <- x
        input_file <- methy(x)
    } else {
        return_nmr <- FALSE
        input_file <- x
    }

    assertthat::assert_that(fs::file_exists(input_file))
    assertthat::assert_that(
        input_file != output_file,
        msg = "target file name must differ from original methylation file."
    )

    assertthat::assert_that(
        tools::file_ext(output_file) == "bgz",
        msg = "output_file must end with .bgz extension."
    )

    output_tsv <- stringr::str_remove(output_file, ".bgz")

    assert_that(is.dir(fs::path_dir(output_tsv)))
    file.create(path.expand(output_tsv))
    assert_that(is.writeable(output_tsv))

    lines_read <- 0
    lines_kept <- 0

    filter_fn <- function(x) {
        dplyr::filter(x, ...)
    }
    writer_fn <- function(x, i) {
        lines_read <<- lines_read + nrow(x)
        filtered_x <- filter_fn(x)
        lines_kept <<- lines_kept + nrow(filtered_x)
        data.table::fwrite(
            filtered_x,
            file = output_tsv,
            sep = "\t",
            append = TRUE,
            scipen = 999L
        )
    }
    readr::local_edition(1) # temporary fix for vroom bad value
    readr::read_tsv_chunked(
        input_file,
        callback = readr::SideEffectChunkCallback$new(writer_fn),
        col_names = methy_col_names(),
        col_types = methy_col_types()
    )

    # assume filter does not change order of rows, therefore data should
    # still be sorted after filtering
    tabix_compress(x = output_tsv, index = TRUE)
    file.remove(output_tsv)

    kept_pct <- scales::percent(lines_kept / lines_read, accuracy = 0.01)
    lines_kept <- scales::comma(lines_kept)
    lines_read <- scales::comma(lines_read)

    message(glue::glue(
        "{lines_kept} of {lines_read} ({kept_pct}) entries kept after filtering"
    ))

    index_name <- paste0(output_file, ".tbi")
    message(glue::glue("results written to '{output_file}' along with index file '{index_name}'"))

    if (return_nmr) {
        methy(nmr) <- output_file
        invisible(nmr)
    } else {
        invisible(output_file)
    }
}
