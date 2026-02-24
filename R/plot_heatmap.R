subsample_read_groups <- function(methy_data, subsample, group_col = "group") {
    subsample_fn <- function(x, key, subsample) {
        if (nrow(x) > subsample) dplyr::sample_n(x, subsample) else x
    }

    methy_data %>%
        dplyr::nest_by(group = .data[[group_col]], read_group = .data$read_group) %>%
        dplyr::group_by(.data$group) %>%
        dplyr::group_modify(subsample_fn, subsample = subsample) %>%
        tidyr::unnest("data")
}

build_heatmap_compact <- function(methy_data) {
    ggplot2::ggplot(
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
            axis.text.y = ggplot2::element_blank(),
            panel.grid.minor = ggplot2::element_blank(),
            panel.grid.major = ggplot2::element_blank()
        ) +
        ggplot2::xlab("Site")
}

build_heatmap_to_scale <- function(methy_data, read_data) {
    ggplot2::ggplot(methy_data, aes(y = .data$read_group)) +
        ggplot2::geom_errorbarh(
            ggplot2::aes(xmin = .data$start, xmax = .data$end),
            data = read_data,
            alpha = 1,
            color = "darkgray",
            linewidth = 0.8,
            width = 0
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

sort_read_groups_by_methy <- function(methy_data) {
    read_group_levels <- methy_data %>%
        dplyr::group_by(.data$group, .data$read_group) %>%
        dplyr::summarise(mean_mod = mean(.data$mod_prob, na.rm = TRUE), .groups = "drop") %>%
        dplyr::arrange(.data$group, .data$mean_mod) %>%
        dplyr::pull(.data$read_group)

    methy_data$read_group <- factor(methy_data$read_group, levels = read_group_levels)
    methy_data
}

plot_methy_data_heatmap <- function(
    methy_data,
    pos_style,
    subsample,
    group_col = "group"
) {
    methy_data <- subsample_read_groups(methy_data, subsample, group_col)
    methy_data <- sort_read_groups_by_methy(methy_data)

    if (pos_style == "compact") {
        p <- build_heatmap_compact(methy_data)
    } else {
        read_data <- methy_data %>%
            dplyr::group_by(.data$read_name) %>%
            dplyr::summarise(start = min(.data$pos), end = max(.data$pos)) %>%
            dplyr::inner_join(
                dplyr::select(methy_data, "read_name", "read_group", "group"),
                by = "read_name"
            )
        p <- build_heatmap_to_scale(methy_data, read_data)
    }

    p + ggplot2::facet_wrap(~group, scales = "free_y", ncol = 1, strip.position = "right")
}
