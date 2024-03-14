# nocov start
# code adapted from zzz.R of readr
.onLoad <- function(libname, pkgname) {
    opt <- options()
    new_opt <- list(
        NanoMethViz.site_filter = 1L,
        NanoMethViz.highlight_col = "grey50"
    )

    to_set <- ! names(new_opt) %in% names(opt)

    if (any(to_set)) options(new_opt[to_set])
    invisible()
}
# nocov end
