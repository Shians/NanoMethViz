# plot methylation data queried using NanoMethViz
plot_methylation_data <- function(
    methy_data,
    sample_anno,
    chr,
    start,
    end,
    title,
    read_anno = NULL,
    group_col = "group",
    palette_col = ggplot2::scale_colour_brewer(palette = "Set1"),
    anno_regions = NULL,
    binary_threshold = NULL,
    avg_method = c("mean", "median"),
    spaghetti = FALSE,
    points = FALSE,
    smoothing_window = 2000,
    highlight_col = getOption("NanoMethViz.highlight_col", "grey50"),
    line_size = 1,
    mod_scale = c(0, 1),
    span = NULL
) {
    if (!missing("span")) {
        warning("the 'span' argument has been deprecated, please use 'smoothing_window' instead")
    }
    # assign averaging method
    avg_method <- match.arg(avg_method)
    avg_func <- switch(
        avg_method,
        mean = confidence_weighted_mean,
        median = median
    )

    # if there are regions to annotate
    if (!is.null(anno_regions)) {
        # filter annotation regions to be within plotting region
        anno_regions <- anno_regions %>%
            dplyr::filter(
                .data$chr == unique(methy_data$chr),
                .data$end >= min(methy_data$pos),
                .data$start <= max(methy_data$pos)
            )
    }

    plot_data <- methy_data %>%
        binarise_statistics(binary_threshold = binary_threshold) %>%
        join_samples_anno(sample_anno) %>%
        add_read_anno(read_anno)

    # set up plot
    # rug first so it appears on the bottom layer
    p <- ggplot(plot_data, aes(x = .data$pos, col = .data[[group_col]])) +
        ggplot2::geom_rug(aes(col = NULL), sides = "b")

    # add annotated regions
    if (!is.null(anno_regions)) {
        for (i in seq_len(nrow(anno_regions))) {
            region <- anno_regions[i, ]
            p <- p +
                ggplot2::annotate(
                    "rect",
                    xmin = region$start,
                    xmax = region$end,
                    ymin = -Inf,
                    ymax = Inf,
                    alpha = 0.2,
                    fill = highlight_col
                )
        }
    }

    # add points
    if (points) {
        p <- p +
            ggplot2::geom_point(
                aes(y = .data$mod_prob),
                alpha = 0.75
            )
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
        dplyr::group_by(.data[[group_col]], .data$pos) %>%
        dplyr::summarise(
            mod_prob = avg_func(.data$mod_prob)
        ) %>%
        dplyr::mutate(
            mod_prob = rolling_average(
                .data$mod_prob,
                .data$pos,
                smoothing_window = smoothing_window
            )
        )

    p <- p +
        ggplot2::geom_line(
            aes(y = .data$mod_prob),
            data = plot_data_smooth,
            na.rm = TRUE,
            linewidth = line_size
        )

    # add auxiliary elements and style
    if (is.numeric(mod_scale) && length(mod_scale) == 2) {
        p <- p + ggplot2::ylim(mod_scale[1], mod_scale[2])
    }
    p +
        ggplot2::ggtitle(title) +
        ggplot2::xlab(chr) +
        palette_col +
        ggplot2::theme_bw() +
        ggplot2::theme(panel.grid.minor = ggplot2::element_blank()) +
        ggplot2::scale_x_continuous(labels = scales::label_number(scale_cut = scales::cut_si("b")))
}

binarise_statistics <- function(plot_data, binary_threshold) {
    if (!is.null(binary_threshold)) {
        # if binary threshold is provided, convert probabilities to binary values
        plot_data$mod_prob <- as.numeric(
            sigmoid(plot_data$statistic) > binary_threshold
        )
    }
    plot_data
}

join_samples_anno <- function(plot_data, sample_anno) {
    plot_data %>%
        dplyr::inner_join(sample_anno, by = "sample", multiple = "all")
}

add_read_anno <- function(plot_data, read_anno) {
    # incorporate read annotation if available
    if (!is.null(read_anno)) {
        # remove any duplicate columns and use plot_data as reference
        plot_data_names <- setdiff(names(plot_data), "read_name")
        read_anno <- dplyr::select(read_anno, -dplyr::any_of(plot_data_names))
        plot_data <- dplyr::inner_join(plot_data, read_anno, by = "read_name") %>%
            tidyr::drop_na()
    }

    plot_data
}

confidence_weighted_mean <- function(x, threshold = 0.5) {
    weights <- abs(x - threshold)^2
    weights <- weights / sum(weights)
    sum(x * weights)
}

rolling_average <- function(y, x, smoothing_window = 2000, neighbours = smoothing_window / 10) {
    tricube_kern <- function(x) {
        ifelse(abs(x) < 1, (1 - abs(x))^3, 0)
    }

    n <- length(y)
    y_smooth <- rep(NA, n)
    for (i in 1:n) {
        start <- max(1, i - neighbours)
        end <- min(n, i + neighbours)
        dists <- abs(x[start:end] - x[i])
        weights <- tricube_kern(dists / smoothing_window)
        weights <- weights / sum(weights)

        y_smooth[i] <- sum(y[start:end] * weights)
    }

    y_smooth
}
