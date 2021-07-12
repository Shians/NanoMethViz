StatLowess <- ggplot2::ggproto("StatLowess", ggplot2::Stat,
    required_aes = c("x", "y"),

    compute_group = function(
        data,
        scales,
        params,
        span = 2/3
    ) {
        lowess_fit <- lowess(data$x, data$y, f = span)
        as.data.frame(lowess_fit)
    }
)

stat_lowess <- function(
    mapping = NULL,
    data = NULL,
    geom = "line",
    position = "identity",
    na.rm = TRUE,
    show.legend = NA,
    inherit.aes = TRUE,
    span = 2/3,
    ...
) {
    ggplot2::layer(
        stat = StatLowess,
        data = data,
        mapping = mapping,
        geom = geom,
        position = position,
        show.legend = show.legend,
        inherit.aes = inherit.aes,
        params = list(
            span = span,
            na.rm = na.rm,
            ...
        )
    )
}
