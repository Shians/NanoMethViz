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
