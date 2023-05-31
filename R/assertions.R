#' @importFrom purrr map_lgl
assert_readable <- function(x) {
    readable <- purrr::map_lgl(x, fs::file_exists)

    if (any(!readable)) {
        unreadables <- x[!readable]
        if (length(unreadables) == 1) {
            stop(glue::glue("Path '{unreadables}' does not exist"))
        } else {
            unreadables <- paste(glue::glue("'{unreadables}'"), collapse = ", ")
            stop(glue::glue("Paths {unreadables} do not exist"))
        }
    }
}

#' @importFrom purrr map_lgl
assert_has_index <- function(x) {
    has_index <- purrr::map_lgl(paste0(x, ".bai"), fs::file_exists)

    if (any(!has_index)) {
        no_index <- x[!has_index]
        if (length(no_index) == 1) {
            stop(glue::glue("files '{no_index}' does not have bam index"))
        } else {
            no_index <- paste(glue::glue("'{no_index}'"), collapse = ", ")
            stop(glue::glue("files {no_index} do not have bam index"))
        }
    }
}
