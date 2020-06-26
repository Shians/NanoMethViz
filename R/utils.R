.get_ggplot_range_x <- function(x) {
    # get x-axis range from a ggplot object
    # returns c(x_min, x_max)
    ggplot2::ggplot_build(x)$layout$panel_scales_x[[1]]$range$range
}

# create a list where the nth element contains the nth values of the original
# vectors
vec_zip <- function(..., .names = NULL) {
    x <- do.call(data.frame, list(..., stringsAsFactors = FALSE))
    setNames(split(x, 1:nrow(x)), .names)
}

extract_file_names <- function(x) {
    fs::path_ext_remove(fs::path_file(x))
}

logit <- function(p) {
    log(p / (1-p))
}

assert_has_columns <- function(x, cols) {
    if (!all(cols %in% colnames(x))) {
        stop(glue::glue(
            "columns missing from {input}: {missing_cols}",
            input = deparse(substitute(x)),
            missing_cols = paste(setdiff(cols, colnames(x)), collapse = ", ")
        ))
    }
}
