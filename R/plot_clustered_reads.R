plot_clustered_reads <- function(x, chr, start, end, min_pts = 5, title = glue::glue("{chr}:{start}-{end}")) {
    methy_data <- query_methy(x, chr, start, end) %>%
        dplyr::filter(.data$pos > start & .data$pos < end) %>%
        dplyr::inner_join(samples(x), by = "sample")

    cluster_res <- cluster_reads(x, chr, start, end, min_pts = min_pts)

    append_read_group <- function(x, k) {
        x$read_group <- paste0(k, stacked_interval_inds(x))
        x
    }

    read_data <- methy_data %>%
        dplyr::group_by(.data$read_name) %>%
        dplyr::summarise(start = min(.data$pos), end = max(.data$pos))

    group_data <- methy_data %>%
        dplyr::group_by(.data$read_name) %>%
        dplyr::summarise(group = unique(.data$group))

    read_data <- dplyr::inner_join(
        read_data,
        group_data,
        by = "read_name",
        multiple = "all",
        suffix = c("_read", "_group")
    )

    read_data <- read_data %>%
        dplyr::group_by(.data$group) %>%
        dplyr::group_modify(append_read_group) %>%
        dplyr::ungroup()

    methy_data <- dplyr::inner_join(
        methy_data,
        read_data,
        by = dplyr::join_by("read_name", "group")
    )

    p1 <- plot_methylation_internal(
        methy_data,
        sample_anno = samples(x),
        read_anno = cluster_res,
        chr = chr,
        start = start,
        end = end,
        group_col = "cluster_id",
        title = title,
        points = TRUE
    ) +
        ggplot2::coord_cartesian(xlim = c(start, end), expand = FALSE) +
        ggplot2::scale_x_continuous(labels = scales::label_number(scale_cut = scales::cut_si("b")))

    heatmap_data <- dplyr::inner_join(
        methy_data,
        cluster_res,
        by = dplyr::join_by("read_name", "strand", "start", "end")
    )

    p2 <- plot_heatmap_internal(
        methy_data = heatmap_data,
        pos_style = "to_scale",
        subsample = 30,
        group_col = "cluster_id"
    ) +
        ggplot2::coord_cartesian(xlim = c(start, end), expand = TRUE) +
        ggplot2::scale_x_continuous(labels = scales::label_number(scale_cut = scales::cut_si("b")), expand = ggplot2::expansion(0, 0))


    p1 / p2
}
