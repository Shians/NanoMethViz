# Skeleton for merge function for two methylation files, currently Unix only
merge_methy_files <- function(inputs, output) {
    if (.Platform$OS.type != "unix") {
        stop(glue::glue(
            "merge_methy_files is only supported on Unix systems, as it relies ",
            "on the 'sort -m' command which is unavailable on this platform.\n",
            "Please run the merge on a Unix system."
        ))
    }

    assertthat::assert_that(
        all(fs::is_file(inputs)),
        all(fs::file_exists(inputs))
    )

    assertthat::assert_that(
        stringr::str_detect(output, "[.]bgz$"),
        msg = "'output' must end in .bgz"
    )

    output <- stringr::str_remove(output, ".bgz")

    temp_files <- purrr::map_chr(seq_along(inputs), ~tempfile())

    temp_merged <- tempfile()

    purrr::walk2(
        inputs,
        temp_files,
        ~R.utils::gunzip(.x, destname = .y, remove = FALSE)
    )

    args <- c(
        "-m", "-k2,3V",
        "-o", shQuote(temp_merged),
        shQuote(temp_files)
    )
    status <- system2("sort", args)
    if (status != 0) {
        stop(glue::glue(
            "Failed to merge methylation files: {paste(inputs, collapse = ', ')}.\n",
            "The 'sort' command exited with status {status}.\n",
            "Please check that 'sort' is available and the files share the same sort key convention."
        ))
    }

    fs::file_copy(temp_merged, output, overwrite = TRUE)
    if (fs::file_exists(paste0(output, ".bgz.tbi"))) {
        fs::file_delete(paste0(output, ".bgz.tbi"))
    }

    tabix_compress(output)
    fs::file_delete(output)
}
