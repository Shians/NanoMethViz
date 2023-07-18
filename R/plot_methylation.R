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
    span = NULL,
    line_size = 2,
    mod_scale = c(0, 1)
) {
    # assign averaging method
    avg_method <- match.arg(avg_method)
    avg_func <- switch(
        avg_method,
        mean = mean,
        median = median
    )

    if (!is.null(anno_regions)) {
        # filter annotation regions to be within plotting region
        anno_regions <- anno_regions %>%
            dplyr::filter(
                .data$chr == unique(methy_data$chr),
                .data$end >= min(methy_data$pos),
                .data$start <= max(methy_data$pos)
            )
    }

    # assign default span, roughly equal to 3000 bp with a maximum of 0.4
    if (is.null(span)) {
        span <- min(3000 / (end - start), 0.4)
    }

    # extract group information and convert probabilities
    plot_data <- methy_data
    if (!is.null(binary_threshold)) {
        # if binary threshold is provided, convert probabilities to binary values
        plot_data$mod_prob <- as.numeric(
            sigmoid(methy_data$statistic) > binary_threshold
        )
    }

    plot_data <- plot_data %>%
        dplyr::inner_join(sample_anno, by = "sample", multiple = "all")

    # incorporate read annotation if available
    if (!is.null(read_anno)) {
        # remove any duplicate columns and use plot_data as reference
        plot_data_names <- setdiff(names(plot_data), "read_name")
        read_anno <- dplyr::select(read_anno, -dplyr::any_of(plot_data_names))
        plot_data <- dplyr::inner_join(plot_data, read_anno, by = "read_name") %>%
            tidyr::drop_na()
    }

    # set up plot
    # rug first so it appears on the bottom layer
    p <- ggplot(plot_data, aes(x = .data$pos, col = .data[[group_col]])) +
        ggplot2::geom_rug(aes(col = NULL), sides = "b")

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
        dplyr::summarise(mod_prob = avg_func(.data$mod_prob))

    p <- p +
        stat_lowess(
            aes(y = .data$mod_prob),
            data = plot_data_smooth,
            span = span,
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
        ggplot2::scale_x_continuous(labels = scales::label_number(scale_cut = scales::cut_si("b")))
}

# wrapper function for plotting a single feature
# responsible for setting up the query and calling plot_methylation_data
plot_feature <- function(
    feature,
    title = "",
    methy,
    sample_anno,
    anno_regions = NULL,
    window_size = c(0, 0),
    binary_threshold = NULL,
    avg_method = c("mean", "median"),
    spaghetti = FALSE,
    span = NULL,
    palette = ggplot2::scale_colour_brewer(palette = "Set1"),
    line_size = 2,
    mod_scale = c(0, 1)
) {
    avg_method <- match.arg(avg_method)

    chr <- feature$chr
    start <- feature$start
    end <- feature$end

    feature_width <- end - start
    window_left <- window_size[1]
    window_right <- window_size[2]
    xlim <- c(start - window_left, end + window_right)

    methy_data <-
        query_methy(
            methy,
            chr,
            floor(start - window_left * 1.1),
            ceiling(end + window_right * 1.1),
            simplify = TRUE)

    if (nrow(methy_data) == 0) {
        warning("no methylation data in region")
        return(ggplot() + theme_void())
    }

    methy_data <- methy_data %>%
        dplyr::select(-"strand") %>%
        tibble::as_tibble()

    plot_methylation_data(
        methy_data = methy_data,
        start = start,
        end = end,
        chr = chr,
        title = title,
        anno_regions = anno_regions,
        binary_threshold = binary_threshold,
        avg_method = avg_method,
        spaghetti = spaghetti,
        sample_anno = sample_anno,
        span = span,
        palette_col = palette,
        line_size = line_size,
        mod_scale = mod_scale
    )
}
