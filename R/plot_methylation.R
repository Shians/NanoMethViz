plot_methylation_internal <- function(
    methy_data,
    sample_anno,
    chr,
    start,
    end,
    title,
    anno_regions = NULL,
    spaghetti = FALSE,
    span = NULL
    ) {

    if (!is.null(anno_regions)) {
        anno_regions <- anno_regions %>%
            dplyr::filter(
                .data$chr == unique(methy_data$chr),
                .data$end >= min(methy_data$pos),
                .data$start <= max(methy_data$pos)
            )
    }

    if (is.null(span)) {
        span <- min(8000 / (end - start), 0.4)
    }

    # extract group information and convert probabilities
    plot_data <- methy_data %>%
        dplyr::inner_join(sample_anno, by = "sample") %>%
        dplyr::mutate(
            mod_prob = e1071::sigmoid(.data$statistic)
        )

    # set up plot
    p <- ggplot(plot_data, aes(x = .data$pos, col = .data$group))

    # add annotated regions
    if (!is.null(anno_regions)) {
        for (i in seq_len(nrow(anno_regions))) {
            region <- anno_regions[i,]
            p <- p +
                ggplot2::annotate(
                    "rect",
                    xmin = region$start,
                    xmax = region$end,,
                    ymin = -Inf,
                    ymax = Inf,
                    alpha = 0.2
                )
        }
    }

    # add spaghetti
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

    # add smoothed line
    plot_data_smooth <- plot_data %>%
        dplyr::group_by(.data$group, .data$pos) %>%
        dplyr::summarise(mod_prob = mean(.data$mod_prob))

    p <- p +
        ggplot2::stat_smooth(
            aes(y = .data$mod_prob, fill = .data$group),
            data = plot_data_smooth,
            geom = "smooth",
            method = "loess",
            span = span,
            na.rm = TRUE,
            size = 3,
            formula = y ~ x
        )

    # add auxiliary elements and style
    x_min <- max(min(plot_data$pos), start)
    x_max <- min(max(plot_data$pos), end)
    p +
        ggplot2::geom_rug(aes(col = NULL), sides = "b") +
        ggplot2::ggtitle(title) +
        ggplot2::xlab(chr) +
        ggplot2::coord_cartesian(ylim = c(0, 1), clip = "on") +
        ggplot2::scale_x_continuous(
            breaks = c(x_min, x_max),
            labels = scales::comma(c(x_min, x_max))) +
        ggplot2::scale_color_brewer(palette = "Set1") +
        ggplot2::scale_fill_brewer(palette = "Set1") +
        ggthemes::theme_tufte()
}

plot_feature <- function(
    feature,
    title = "",
    methy,
    sample_anno,
    anno_regions = NULL,
    window_prop = c(0.3, 0.3),
    spaghetti = TRUE,
    span = NULL
    ) {

    chr <- feature$chr
    start <- feature$start
    end <- feature$end
    window_left <- (end - start) * window_prop[1]
    window_right <- (end - start) * window_prop[2]

    methy_data <-
        query_methy(methy, chr, start - window_left, end + window_right) %>%
        dplyr::bind_rows() %>%
        dplyr::select(-"strand") %>%
        tibble::as_tibble()


    if (nrow(methy_data) == 0) {
        warning("no methylation data in region")
        return(ggplot() + theme_void())
    }

    plot_methylation_internal(
        methy_data = methy_data,
        start = start,
        end = end,
        chr = chr,
        title = title,
        anno_regions = anno_regions,
        spaghetti = spaghetti,
        sample_anno = sample_anno,
        span = span
    )
}
