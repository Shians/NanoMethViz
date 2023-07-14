plot_methy_data_heatmap <- function(
        methy_data,
        pos_style,
        subsample,
        group_col = "group",
        mod_scale = c(0, 1)
) {
    # subsample reads down to a certain number of read groups
    subsample_groups <- function(x, key, subsample) {
        if (nrow(x) > subsample) {
            x <- dplyr::sample_n(x, subsample)
        }
        x
    }

    methy_data <- methy_data %>%
        dplyr::nest_by(group = .data[[group_col]], read_group = .data$read_group) %>%
        dplyr::group_by(.data$group) %>%
        dplyr::group_modify(subsample_groups, subsample = subsample)
    methy_data <- tidyr::unnest(methy_data, "data")

    read_data <- methy_data %>%
        dplyr::group_by(.data$read_name) %>%
        dplyr::summarise(start = min(.data$pos), end = max(.data$pos)) %>%
        dplyr::inner_join(
            dplyr::select(methy_data, "read_name", "read_group", "group"),
            by = "read_name"
        )

    if (pos_style == "compact") {
        # only plots sites with measured modification, evenly spaced
        p <- ggplot2::ggplot(
            methy_data,
            aes(
                x = factor(.data$pos),
                y = .data$read_group,
                fill = .data$mod_prob
            )
        ) +
            ggplot2::geom_raster() +
            heatmap_fill_scale +
            ggplot2::theme(
                axis.ticks.x = ggplot2::element_blank(),
                axis.text.x = ggplot2::element_blank(),
                axis.ticks.y = ggplot2::element_blank(),
                axis.title.y = ggplot2::element_blank(),
                axis.text.y = ggplot2::element_blank()
            ) +
            ggplot2::xlab("Site")
    } else if (pos_style == "to_scale") {
        # plots all sites in range, evenly spaced with square geoms
        # data will overlap
        p <- ggplot2::ggplot(methy_data, aes(y = .data$read_group)) +
            ggplot2::geom_errorbarh(
                ggplot2::aes(xmin = .data$start, xmax = .data$end),
                data = read_data,
                alpha = 1,
                color = "darkgray",
                linewidth = 0.8,
                height = 0
            ) +
            ggplot2::geom_point(aes(x = .data$pos, col = .data$mod_prob), alpha = 1, shape = 15) +
            ggplot2::scale_x_continuous(
                labels = scales::label_number(scale_cut = scales::cut_si("b")),
                expand = ggplot2::expansion(0, 0)
            ) +
            heatmap_col_scale +
            theme_methy_heatmap +
            ggplot2::xlab("Position")
    }

    p + ggplot2::facet_wrap(~group, scales = "free_y", ncol = 1, strip.position = "right")
}
