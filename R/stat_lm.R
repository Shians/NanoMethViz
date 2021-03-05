StatLm <- ggplot2::ggproto("StatLm", Stat,
    required_aes = c("x", "y"),

    compute_group = function(data, scales, params, n = 20) {
        if (nrow(data) <= 1) {
            return(data.frame(x = NULL, y = NULL))
        }

        poly_deg <- min(6, nrow(data) - 1)
        rng <- range(data$x, na.rm = TRUE)
        grid <- data.frame(x = seq(rng[1], rng[2], length = n))

        mod <- lm(y ~ poly(x, poly_deg), data = data)
        grid$y <- predict(mod, newdata = grid)

        grid
    }
)

stat_lm <- function(
    mapping = NULL, data = NULL, geom = "line",
    position = "identity", na.rm = FALSE, show.legend = NA,
    inherit.aes = TRUE, n = 50,
    ...
) {
    ggplot2::layer(
        stat = StatLm, data = data, mapping = mapping, geom = geom,
        position = position, show.legend = show.legend, inherit.aes = inherit.aes,
        params = list(n = n, na.rm = na.rm, ...)
    )
}
