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
