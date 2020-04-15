# get x-axis range from a ggplot object
.get_ggplot_range_x <- function(x) {
    # returns c(x_min, x_max)
    ggplot2::ggplot_build(x)$layout$panel_scales_x[[1]]$range$range
}