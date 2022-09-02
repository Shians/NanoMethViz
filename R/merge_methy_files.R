# Skeleton for merge function for two methylation files, currently Unix only
merge_methy_files <- function(inputs, output) {
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

    files_str <- paste(temp_files, collapse = " ")
    cmd <- glue::glue("sort -m -k2,3V -o {temp_merged} {files_str}")
    system(cmd)

    fs::file_copy(temp_merged, output, overwrite = TRUE)
    if (fs::file_exists(paste0(output, ".bgz.tbi"))) {
        fs::file_delete(paste0(output, ".bgz.tbi"))
    }

    tabix_compress(output)
    fs::file_delete(output)
}

