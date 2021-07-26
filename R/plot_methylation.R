plot_methylation_internal <- function(
    methy_data,
    sample_anno,
    chr,
    start,
    end,
    xlim,
    title,
    palette_col = ggplot2::scale_colour_brewer(palette = "Set1"),
    anno_regions = NULL,
    spaghetti = FALSE,
    span = NULL
) {
    if (!is.null(anno_regions)) {
        # filter annotation regions to be within plotting region
        anno_regions <- anno_regions %>%
            dplyr::filter(
                .data$chr == unique(methy_data$chr),
                .data$end >= min(methy_data$pos),
                .data$start <= max(methy_data$pos)
            )
    }

    if (is.null(span)) {
        span <- min(3000 / (end - start), 0.4)
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
                    xmax = region$end,
                    ymin = -Inf,
                    ymax = Inf,
                    alpha = 0.2
                )
        }
    }

    # add spaghetti
    if (spaghetti) {
        p <- p +
            stat_lm(
                aes(y = .data$mod_prob, group = .data$read_name),
                alpha = 0.25,
                na.rm = TRUE
            )
    }

    # add smoothed line
    plot_data_smooth <- plot_data %>%
        dplyr::group_by(.data$group, .data$pos) %>%
        dplyr::summarise(mod_prob = mean(.data$mod_prob))

    p <- p +
        stat_lowess(
            aes(y = .data$mod_prob),
            data = plot_data_smooth,
            span = span,
            na.rm = TRUE,
            size = 2
        )

    # add auxiliary elements and style
    x_min <- max(min(plot_data$pos), start)
    x_max <- min(max(plot_data$pos), end)
    breaks <- c(x_min, round((x_min + x_max) / 2), x_max)
    p +
        ggplot2::geom_rug(aes(col = NULL), sides = "b") +
        ggplot2::ggtitle(title) +
        ggplot2::xlab(chr) +
        ggplot2::coord_cartesian(ylim = c(0, 1), clip = "on") +
        ggplot2::scale_x_continuous(
            limits = xlim,
            breaks = breaks,
            labels = scales::comma(breaks)) +
        palette_col +
        theme_bw()
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

    feature_width <- end - start
    window_left <- feature_width * window_prop[1]
    window_right <- feature_width * window_prop[2]
    xlim <- c(start - window_left, end + window_right)

    methy_data <-
        query_methy(
            methy,
            chr,
            floor(start - window_left * 1.1),
            ceiling(end + window_right * 1.1),
            simplify = TRUE) %>%
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
        xlim = xlim,
        title = title,
        anno_regions = anno_regions,
        spaghetti = spaghetti,
        sample_anno = sample_anno,
        span = span
    )
}
