#' Sort methylation file
#' @keywords internal
#'
#' @param x the path to the methylation file to sort
#'
#' @return invisibly returns path of sorted file
sort_methy_file <- function(x) {
    assert_that(is.readable(x))

    if (.Platform$OS.type == "windows") {
        methy_df <- data.table::fread(
            x,
            header = FALSE,
            col.names = methy_col_names(),
            showProgress = FALSE,
            data.table = FALSE)
        methy_df <- dplyr::arrange(methy_df, .data$chr, .data$pos)
        readr::write_tsv(methy_df, x, col_names = FALSE)
    } else {
        cmd <- glue::glue("sort --compress-program=gzip -k2,3V {x} -o {x}")
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
