plot_gene_annotation <- function(exons_df, plot_start, plot_end) {
    exons_df <- exons_df %>%
    dplyr::mutate(
        uid = factor(paste(gene_id, transcript_id, sep = ".")),
        y_offset = as.integer(uid) - 1
    )

    exons_count <- exons_df %>%
        dplyr::group_by(uid, y_offset) %>%
        dplyr::summarise(exons = n())

    gap <- exons_df %>%
        dplyr::inner_join(exons_count, by = c("uid", "y_offset")) %>%
        dplyr::filter(exons > 1) %>%
        dplyr::group_by(uid, y_offset, strand) %>%
        dplyr::summarise(
            gap_start = list(end[-length(end)]),
            gap_end = list(start[-1])
        ) %>%
        tidyr::unnest(cols = c(gap_start, gap_end))

    get_gaps <- function(gaps, strand) {
        gaps <- gaps %>%
            split(gap$strand)
        if (!is.null(gaps[[strand]])) {
            gaps[[strand]]
        } else {
            tibble::tibble(
                uid = character(0),
                y_offset = numeric(0),
                strand = character(0),
                gap_start = integer(0),
                gap_end = integer(0)
            )
        }
    }

    gap_pos <- get_gaps(gap, "+")
    gap_neg <- get_gaps(gap, "-")
    gap_none <- get_gaps(gap, "*")

    gene_middle <- exons_df %>%
        dplyr::group_by(gene_id, symbol, transcript_id, y_offset) %>%
        dplyr::summarise(gene_middle = (min(start) + max(end)) / 2)

    ggplot2::ggplot() +
        ggplot2::theme_void() +
        ggplot2::xlim(plot_start, plot_end) +
        ggplot2::geom_rect(ggplot2::aes(xmin = start, xmax = end, ymin = y_offset + 0.05, ymax = y_offset + 0.55), data = exons_df, fill = "#696969") +
        ggplot2::geom_segment(
            ggplot2::aes(
                x = gap_start,
                xend = gap_end,
                y = y_offset + 0.275,
                yend = y_offset + 0.275
            ),
            data = gap_pos,
            arrow = grid::arrow(
                type = "closed",
                length = unit(0.035, "npc")
            )
        ) +
        ggplot2::geom_segment(
            ggplot2::aes(
                x = gap_end,
                xend = gap_start,
                y = y_offset + 0.275,
                yend = y_offset + 0.275
            ),
            data = gap_neg,
            arrow = grid::arrow(
                type = "closed",
                length = unit(0.035, "npc")
            )
        ) +
        ggplot2::geom_segment(
            ggplot2::aes(
                x = gap_end,
                xend = gap_start,
                y = y_offset + 0.275,
                yend = y_offset + 0.275),
            data = gap_none) +
        ggplot2::geom_text(aes(x = gene_middle, y = y_offset + 0.8, label = symbol), data = gene_middle, hjust = "center", size = rel(3.5)) +
        ggplot2::ylim(0, 1 + max(exons_df$y_offset))
}
