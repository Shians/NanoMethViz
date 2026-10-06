#' Sort methylation file
#' @keywords internal
#'
#' @param x the path to the methylation file to sort
#'
#' @return invisibly returns path of sorted file
sort_methy_file <- function(x) {
    assert_readable(x)

    if (.Platform$OS.type == "windows") {
        methy_df <- readr::read_tsv(
            x,
            col_names = methy_col_names(),
            col_types = methy_col_types()
        )
        methy_df <- dplyr::arrange(methy_df, .data$chr, .data$pos)
        readr::write_tsv(methy_df, x, col_names = FALSE, progress = FALSE)
    } else {
        args <- c(
            "--compress-program=gzip",
            "-S", "4G",
            "-k2,3V",
            shQuote(x),
            "-o", shQuote(x)
        )
        status <- system2("sort", args)
        if (status != 0) {
            stop(glue::glue(
                "Failed to sort methylation file '{x}'.\n",
                "The 'sort' command exited with status {status}.\n",
                "Please check that 'sort' is available and the file is readable."
            ))
        }
    }

    invisible(x)
}

tabix_compress <- function(x, index = TRUE) {
    assert_readable(x)

    f <- Rsamtools::bgzip(x, overwrite = TRUE)
    if (index) {
        tabix_index(f)
    }

    f
}

tabix_index <- function(x) {
    assert_readable(x)

    Rsamtools::indexTabix(x, seq = 2, start = 3, end = 3)
}

#' Convert methylation file to tabix format
#' @keywords internal
#'
#' @param x the path to the sorted methylation file
#'
#' @return invisibly returns the path to the tabix file
raw_methy_to_tabix <- function(x) {
    assert_readable(x)

    bgz_name <- tabix_compress(x)
    tabix_index(bgz_name)

    invisible(bgz_name)
}
