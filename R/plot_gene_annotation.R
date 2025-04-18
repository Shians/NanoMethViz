#' Plot gene annotations with exons and introns
#'
#' @description
#' Creates a ggplot visualization of gene annotations showing exons and intron regions
#' within a specified genomic range.
#'
#' @param exons_df A data frame containing exon information with columns:
#'   \itemize{
#'     \item gene_id: Gene identifier
#'     \item transcript_id: Transcript identifier
#'     \item symbol: Gene symbol
#'     \item strand: Strand direction ("+", "-", or "*")
#'     \item start: Exon start position
#'     \item end: Exon end position
#'   }
#' @param plot_start Numeric, the start position of the plotting region
#' @param plot_end Numeric, the end position of the plotting region
#'
#' @return A ggplot object showing gene annotations. The plot has an attribute
#'   "plot_height" indicating the height needed to display all transcripts.
#'   Returns an empty plot with height 0 if no exons are in the region.
#'
#' @details
#' The function filters exons to the specified region, organizes them by transcript,
#' and visualizes them with connecting lines for introns. Gene symbols are placed
#' near the middle of each transcript. Arrows indicate the strand direction when available.
#'
#' @importFrom dplyr mutate group_by summarise filter inner_join arrange ungroup
#' @importFrom tidyr unnest
#' @importFrom tibble tibble
#' @importFrom ggplot2 ggplot theme_void ylim
#'
#' @noRd
plot_gene_annotation <- function(exons_df, plot_start, plot_end) {
    if (nrow(exons_df) == 0) {
        p <- ggplot2::ggplot() + ggplot2::theme_void()
        attr(p, "plot_height") <- 0
        return(p)
    }

    # remove transcripts outside plot area
    exons_df <- .anno_filter_regions(exons_df, plot_start, plot_end)

    exons_df <- exons_df %>%
        dplyr::mutate(
            uid = factor(paste(.data$gene_id, .data$transcript_id, sep = ".")),
            y_offset = as.integer(.data$uid) - 1
        )

    exons_count <- exons_df %>%
        dplyr::group_by(.data$uid) %>%
        dplyr::summarise(exons = dplyr::n())

    gap <- exons_df %>%
        dplyr::inner_join(exons_count, by = c("uid"), multiple = "all") %>%
        dplyr::filter(.data$exons > 1) %>%
        dplyr::group_by("transcript_id") %>%
        dplyr::arrange(.data$start) %>%
        dplyr::ungroup()

    if (nrow(gap) > 0) {
        gap <- gap %>%
            dplyr::group_by(.data$uid, .data$strand, .data$y_offset) %>%
            dplyr::summarise(
                gap_start = list(.data$end[-length(.data$end)]),
                gap_end = list(.data$start[-1])
            ) %>%
            tidyr::unnest(cols = c("gap_start", "gap_end"))
    } else {
        gap <- tibble::tibble(
            uid = character(),
            y_offset = numeric(),
            strand = character(),
            gap_start = numeric(),
            gap_end = numeric()
        )
    }

    gap_pos <- .anno_get_gaps(gap, "+")
    gap_neg <- .anno_get_gaps(gap, "-")
    gap_none <- .anno_get_gaps(gap, "*")

    region_width <- plot_end - plot_start
    gene_labels <- exons_df %>%
        dplyr::group_by(.data$gene_id, .data$symbol, .data$transcript_id, .data$y_offset, .data$strand) %>%
        dplyr::summarise(
            gene_middle = (min(.data$start) + max(.data$end)) / 2,
            label_pos = .data$gene_middle,
            label_pos = pmin(.data$label_pos, plot_end - region_width * 0.05),
            label_pos = pmax(.data$label_pos, plot_start +  region_width * 0.05)
        )

    exons_df <- .anno_truncate_region(exons_df, plot_start, plot_end, "*")

    gap_pos <- .anno_truncate_region(gap_pos, plot_start, plot_end, "+")
    gap_neg <- .anno_truncate_region(gap_neg, plot_start, plot_end, "-")
    gap_none <- .anno_truncate_region(gap_none, plot_start, plot_end, "*")
    gene_labels <- gene_labels %>%
        dplyr::filter(
            dplyr::between(.data$label_pos, plot_start, plot_end)
        )

    if (length(exons_df$y_offset) > 0) {
        plot_height <- 1 + max(exons_df$y_offset)
    } else {
        plot_height <- 0
    }

    p <- ggplot2::ggplot() +
        ggplot2::theme_void() +
        .anno_connector_arrows(gap_pos) +
        .anno_connector_arrows(gap_neg) +
        .anno_connector_lines(gap_none) +
        .anno_exons(exons_df) +
        .anno_gene_labels(gene_labels) +
        ggplot2::ylim(0, plot_height)

    attr(p, "plot_height") <- plot_height

    p
}

# helper functions ----
.anno_exons <- function(exons_df) {
    ggplot2::geom_rect(
        ggplot2::aes(
            xmin = .data$start,
            xmax = .data$end,
            ymin = .data$y_offset + 0.05,
            ymax = .data$y_offset + 0.55
        ),
        data = exons_df,
        fill = "#696969",
        colour = "black"
    )
}

.anno_connector_arrows <- function(gaps) {
    list(
        # first half with arrow
        ggplot2::geom_segment(
            ggplot2::aes(
                x = .data$start,
                xend = (.data$start + .data$end)/2,
                y = .data$y_offset + 0.275,
                yend = .data$y_offset + 0.275
            ),
            lineend = "butt",
            linejoin = "mitre",
            data = gaps,
            arrow = grid::arrow(
                type = "open",
                length = ggplot2::unit(5, "points")
            )
        ),
        # second half without arrow
        ggplot2::geom_segment(
            ggplot2::aes(
                x = (.data$start + .data$end)/2,
                xend = .data$end,
                y = .data$y_offset + 0.275,
                yend = .data$y_offset + 0.275
            ),
            lineend = "butt",
            linejoin = "mitre",
            data = gaps
        )
    )
}

.anno_connector_lines <- function(gaps) {
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

.anno_gene_labels <- function(gene_labels) {
    gene_labels$symbol[gene_labels$strand == "+"] <- paste(
        gene_labels$symbol[gene_labels$strand == "+"],
        ">"
    )

    gene_labels$symbol[gene_labels$strand == "-"] <- paste(
        "<",
        gene_labels$symbol[gene_labels$strand == "-"]
    )

    ggplot2::geom_text(
        ggplot2::aes(x = .data$label_pos, y = .data$y_offset + 0.8, label = .data$symbol),
        data = gene_labels,
        hjust = "center",
        size = ggplot2::rel(2.5)
    )
}

.anno_filter_regions <- function(exons_df, plot_start, plot_end) {
    transcripts <- exons_df %>%
        dplyr::summarise(
            .by = "transcript_id",
            start = min(.data$start),
            end = max(.data$end)
        )

    transcripts <- transcripts %>%
        dplyr::filter(
            .data$start <= plot_end & .data$end >= plot_start
        )

    exons_df %>%
        dplyr::filter(
            .data$transcript_id %in% transcripts$transcript_id
        )
}

.anno_truncate_region <- function(x, plot_start, plot_end, strand) {
    if (strand == "-") {
        x <- x %>%
            dplyr::filter(.data$end <= plot_end, .data$start >= plot_start)
        x$end[x$end < plot_start] <- plot_start
        x$start[x$start > plot_end] <- plot_end
    } else {
        x <- x %>%
            dplyr::filter(.data$start <= plot_end, .data$end >= plot_start)
        x$start[x$start < plot_start] <- plot_start
        x$end[x$end > plot_end] <- plot_end
    }

    x
}

.anno_get_gaps <- function(gaps, strand = c("+", "-", "*")) {
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
