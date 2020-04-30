plot_gene_annotation <- function(exons_df, plot_start, plot_end) {
    exons_df <- exons_df %>%
    dplyr::mutate(
        uid = factor(paste(.data$gene_id, .data$transcript_id, sep = ".")),
        y_offset = as.integer(.data$uid) - 1
    )

    exons_count <- exons_df %>%
        dplyr::group_by(.data$uid, .data$y_offset) %>%
        dplyr::summarise(exons = dplyr::n())

    gap <- exons_df %>%
        dplyr::inner_join(exons_count, by = c("uid", "y_offset")) %>%
        dplyr::filter(.data$exons > 1) %>%
        dplyr::group_by(.data$uid, .data$y_offset, .data$strand) %>%
        dplyr::summarise(
            gap_start = list(.data$end[-length(.data$end)]),
            gap_end = list(.data$start[-1])
        ) %>%
        tidyr::unnest(cols = c(.data$gap_start, .data$gap_end))

    .get_gaps <- function(gaps, strand = c("+", "-", "*")) {
        strand <- match.arg(strand)

        gaps <- gaps %>%
            split(gaps$strand)

        if (!is.null(gaps[[strand]])) {
            out <- gaps[[strand]]

            if (strand == "-") {
                temp <- out$gap_start
                out$gap_start <- out$gap_end
                out$gap_end <- temp
            }
        } else {
            out <- tibble::tibble(
                uid = character(0),
                y_offset = numeric(0),
                strand = character(0),
                gap_start = integer(0),
                gap_end = integer(0)
            )
        }

        dplyr::rename(
            out,
            start = "gap_start",
            end = "gap_end"
        )
    }

    gap_pos <- .get_gaps(gap, "+")
    gap_neg <- .get_gaps(gap, "-")
    gap_none <- .get_gaps(gap, "*")

    gene_labels <- exons_df %>%
        dplyr::group_by(.data$gene_id, .data$symbol, .data$transcript_id, .data$y_offset, .data$strand) %>%
        dplyr::summarise(gene_middle = (min(.data$start) + max(.data$end)) / 2)

    .exons <- function(exons_df) {
        ggplot2::geom_rect(
            ggplot2::aes(
                xmin = .data$start,
                xmax = .data$end,
                ymin = .data$y_offset + 0.05,
                ymax = .data$y_offset + 0.55
            ),
            data = exons_df,
            fill = "#696969"
        )
    }

    .connector_arrows <- function(gaps) {
        ggplot2::geom_segment(
            ggplot2::aes(
                x = .data$start,
                xend = .data$end,
                y = .data$y_offset + 0.275,
                yend = .data$y_offset + 0.275
            ),
            data = gaps,
            arrow = grid::arrow(
                type = "closed",
                length = ggplot2::unit(5, "points")
            )
        )
    }

    .connector_lines <- function(gaps) {
        ggplot2::geom_segment(
            ggplot2::aes(
                x = .data$end,
                xend = .data$start,
                y = .data$y_offset + 0.275,
                yend = .data$y_offset + 0.275
            ),
            data = gaps
        )
    }

    .gene_labels <- function(gene_labels) {
        gene_labels$symbol[gene_labels$strand == "+"] <- paste(
            gene_labels$symbol[gene_labels$strand == "+"],
            "»"
        )

        gene_labels$symbol[gene_labels$strand == "-"] <- paste(
            "«",
            gene_labels$symbol[gene_labels$strand == "-"]
        )

        ggplot2::geom_text(
            aes(x = .data$gene_middle, y = .data$y_offset + 0.8, label = .data$symbol),
            data = gene_labels,
            hjust = "center",
            size = ggplot2::rel(3.5)
        )
    }

    .truncate_region <- function(x, plot_start, plot_end, strand) {
        if (strand == "-") {
            x <- x %>%
                dplyr::filter(end <= plot_end, start >= plot_start)
            x$end[x$end < plot_start] <- plot_start
            x$start[x$start > plot_end] <- plot_end
        } else {
            x <- x %>%
                dplyr::filter(start <= plot_end, end >= plot_start)
            x$start[x$start < plot_start] <- plot_start
            x$end[x$end > plot_end] <- plot_end
        }

        x
    }

    exons_df <- .truncate_region(exons_df, plot_start, plot_end, "*")
    gap_pos <- .truncate_region(gap_pos, plot_start, plot_end, "+")
    gap_neg <- .truncate_region(gap_neg, plot_start, plot_end, "-")
    gap_none <- .truncate_region(gap_none, plot_start, plot_end, "*")
    gene_labels <- gene_labels %>%
        dplyr::filter(
            dplyr::between(.data$gene_middle, plot_start, plot_end)
        )

    ggplot2::ggplot() +
        ggplot2::theme_void() +
        .exons(exons_df) +
        .connector_arrows(gap_pos) +
        .connector_arrows(gap_neg) +
        .connector_lines(gap_none) +
        .gene_labels(gene_labels) +
        ggplot2::xlim(plot_start, plot_end) +
        ggplot2::ylim(0, 1 + max(exons_df$y_offset))
}
