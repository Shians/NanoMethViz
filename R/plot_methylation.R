plot_methylation_internal <- function(
    methy_data,
    sample_anno,
    chr,
    start,
    end,
    title,
    anno_regions = NULL,
    spaghetti = FALSE,
    prop = TRUE
    ) {
    if (!is.null(anno_regions)) {
        anno_regions <- anno_regions %>%
            dplyr::filter(
                .data$chr == unique(methy_data$chr),
                .data$end >= min(methy_data$pos),
                .data$start <= max(methy_data$pos)
            )
    }

    if (prop) {
        # plot proportions
        plot_data_prop <- methy_data %>%
            dplyr::inner_join(sample_anno, by = "sample") %>%
            dplyr::mutate(
                mod_prob = e1071::sigmoid(.data$statistic)
            ) %>%
            dplyr::group_by(
                .data$group, .data$pos
            ) %>%
            dplyr::summarise(
                mod_prop = sum(.data$mod_prob > 0.5) / dplyr::n()
            )
    }

    # extract group information and convert probabilities
    plot_data <- methy_data %>%
        dplyr::inner_join(sample_anno, by = "sample") %>%
        dplyr::mutate(
            mod_prob = e1071::sigmoid(.data$statistic)
        )

    smooth_span <- min(8000 / (end - start), 0.4)
    if (prop) {
        p <- ggplot(plot_data, aes(x = .data$pos, col = .data$group)) +
            ggplot2::stat_smooth(
                data = plot_data_prop,
                aes(y = .data$mod_prop),
                geom = "smooth",
                method = "loess",
                na.rm = TRUE,
                size = 3,
                span = smooth_span,
                formula = y ~ x
            )
    } else {
        p <- ggplot(plot_data, aes(x = .data$pos, col = .data$group)) +
            ggplot2::stat_smooth(
                aes(y = .data$mod_prob),
                geom = "smooth",
                method = "loess",
                na.rm = TRUE,
                size = 3,
                span = smooth_span,
                formula = y ~ x
            )
    }

    if (spaghetti) {
        p <- p +
            ggplot2::stat_smooth(
                aes(y = .data$mod_prob, group = .data$read_name),
                alpha = 0.25,
                geom = "line",
                method = "loess",
                na.rm = TRUE,
                se = FALSE,
                span = 1,
                formula = y ~ x
            )
    }

    if (!is.null(anno_regions)) {
        for (i in seq_len(nrow(anno_regions))) {
            region <- anno_regions[i,]
            p <- p +
                ggplot2::annotate(
                    "rect",
                    xmin = region$start,
                    xmax = region$end,
                    ymin = 0,
                    ymax = 1,
                    alpha = 0.2,
                    col = "green"
                )
        }
    }

    x_min <- max(min(plot_data$pos), start)
    x_max <- min(max(plot_data$pos), end)

    p +
        ggplot2::geom_rug(aes(col = NULL), sides = "b") +
        ggplot2::ggtitle(title) +
        ggplot2::ylim(0, 1) +
        ggplot2::xlab(chr) +
        ggplot2::scale_x_continuous(
            breaks = c(x_min, x_max),
            labels = scales::comma(c(x_min, x_max))) +
        ggplot2::scale_color_brewer(palette = "Set1") +
        ggplot2::scale_fill_brewer(palette = "Set1") +
        ggthemes::theme_tufte()
}

plot_agg <- function(
    chr,
    start,
    end,
    title,
    methy,
    sample_anno,
    anno_regions = NULL,
    window_prop = 0.3,
    prop = TRUE
    ) {
    window <- (end - start) * window_prop

    methy_data <- query_methy(methy, chr, start - window, end + window) %>%
        dplyr::bind_rows() %>%
        dplyr::select(-"strand", -"modified") %>%
        tibble::as_tibble()

    plot_methylation_internal(
        methy_data = methy_data,
        start = start,
        end = end,
        chr = chr,
        title = title,
        anno_regions = anno_regions,
        spaghetti = FALSE,
        sample_anno = sample_anno,
        prop = prop
    )
}

plot_spaghetti_and_agg <- function(
    chr,
    start,
    end,
    title,
    methy,
    sample_anno,
    anno_regions = NULL,
    window_prop = 0.3
    ) {
    window <- (end - start) * window_prop

    methy_data <- query_methy(methy, chr, start - window, end + window) %>%
        dplyr::bind_rows() %>%
        dplyr::select(-"strand", -"modified") %>%
        tibble::as_tibble()

    plot_methylation_internal(
        methy_data = methy_data,
        start = start,
        end = end,
        chr = chr,
        title = title,
        spaghetti = TRUE,
        sample_anno = sample_anno,
        anno_regions = anno_regions
    )
}

plot_feature <- function(
    feature,
    title = "",
    methy,
    sample_anno,
    anno_regions = NULL,
    window_prop = 0.3,
    spaghetti = TRUE
    ) {
    if (spaghetti) {
        plot_spaghetti_and_agg(
            methy = methy,
            anno_regions = anno_regions,
            feature$chr,
            feature$start,
            feature$end,
            title = title,
            sample_anno = sample_anno,
            window_prop = window_prop
        )
    } else {
        plot_agg(
            methy = methy,
            sample_anno = sample_anno,
            feature$chr,
            feature$start,
            feature$end,
            title = title,
            anno_regions = anno_regions,
            window_prop = window_prop
        )
    }
}
