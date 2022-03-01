expect_nrow <- function(object, n) {
    # 1. Capture object and label
    act <- testthat::quasi_label(rlang::enquo(object), arg = "object")

    # 2. Call expect()
    act$n <- nrow(act$val)
    testthat::expect(
        act$n == n,
        sprintf("%s has nrow %i, not nrow %i.", act$lab, act$n, n)
    )

    # 3. Invisibly return the value
    invisible(act$val)
}

expect_ncol <- function(object, n) {
    # 1. Capture object and label
    act <- testthat::quasi_label(rlang::enquo(object), arg = "object")

    # 2. Call expect()
    act$n <- ncol(act$val)
    testthat::expect(
        act$n == n,
        sprintf("%s has ncol %i, not ncol %i.", act$lab, act$n, n)
    )

    # 3. Invisibly return the value
    invisible(act$val)
}
