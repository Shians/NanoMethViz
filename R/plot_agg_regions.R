#' Plot aggregate regions
#'
#' @param x the NanoMethResult object.
#' @param regions a table of regions containing at least columns chr, strand,
#'   start and end. Any additiona columns can be used for grouping.
#' @param group_col the column to group aggregated trends by. This column can
#'   be in from the regions table or samples(x).
#' @param flank the number of flanking bases to add to each side of each region.
#' @param stranded TRUE if negative strand features should have coordinates
#'   flipped to reflect features like transcription start sites.
#' @param span the span for loess smoothing.
#' @param palette the colour palette used for groups.
#'
#' @return a ggplot object containing the aggregate methylation trend.
#'
#' @examples
#' nmr <- load_example_nanomethresult()
#' gene_anno <- exons_to_genes(NanoMethViz::exons(nmr))
#' plot_agg_regions(nmr, gene_anno)
#' plot_agg_regions(nmr, gene_anno, group_col = "sample")
#' plot_agg_regions(nmr, gene_anno, group_col = "group")
#'
#' @export
plot_agg_regions <- function(
    x,
    regions,
    group_col = NULL,
    flank = 2000,
    stranded = TRUE,
    span = 0.05,
    palette = ggplot2::scale_colour_brewer(palette = "Set1")
) {
    if (!is.null(group_col)) {
        avail_columns <- c(colnames(samples(x)), colnames(regions))
        assertthat::assert_that(
            group_col %in% avail_columns,
            msg = glue::glue("'{group_col}' could not be found in columns of 'regions' or samples(x)")
        )
    }

    # query methylation data
    methy_data <- .get_methy_data(x, regions, flank)
    methy_data <- .scale_methy_data(methy_data, stranded, flank)
    methy_data <- .pos_avg(methy_data)
    # unnest and annotate after position averaging to reduce data burden
    methy_data <- .unnest_anno(methy_data, NanoMethViz::samples(x))
    methy_data <- .bin_avg(methy_data, group_col = group_col)

    # set up plot
    p <- ggplot2::ggplot() +
        ggplot2::ylim(c(0, 1)) +
        ggplot2::theme_minimal() +
        ggplot2::stat_smooth(
            ggplot2::aes_string(
                x = "binned_pos",
                y = "methy_prop",
                group = group_col,
                col = group_col),
            method = "loess",
            formula = "y ~ x",
            span = span,
            na.rm = TRUE,
            se = FALSE,
            data = methy_data) +
        palette

    if (flank == 0) {
        labels <- c("start", "end")
        breaks = c(0, 1)
        limits = c(0, 1)
    } else {
        breaks = c(-.33, 0, 1, 1.33)
        limits = c(-0.33, 1.33)
        kb_marker <- round(flank / 1000, 1)
        labels <- c(glue::glue("-{kb_marker}kb"), "start", "end", glue::glue("+{kb_marker}kb"))
        p <- p +
            ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey80") +
            ggplot2::geom_vline(xintercept = 1, linetype = "dashed", color = "grey80")
    }

    p + ggplot2::coord_cartesian(clip = "off") +
        ggplot2::scale_x_continuous(
            name = "Relative Position",
            breaks = breaks,
            limits = limits,
            labels = labels) +
        ggplot2::ylab("Average Methylation Proportion")
}

.get_methy_data <- function(x, regions, flank) {
  query_row_methy <- function(i, methy, regions, flank) {
      # query a larger region such that smoothing doesn't
      # misbehave around ends
      flank <- flank * 1.1
      query_methy(
          methy,
          regions$chr[i],
          regions$start[i] - flank,
          regions$end[i] + flank)
  }

  methy_data <- purrr::map(
      1:nrow(regions),
      query_row_methy,
      methy = methy(x),
      regions = regions,
      flank = flank
  )

  output <- regions %>%
      mutate(methy_data = methy_data)
}

.scale_methy_data <- function(methy_data, stranded, flank) {
    for (i in 1:nrow(methy_data)) {
        m_data <- methy_data$methy_data[[i]]$pos
        within <- m_data >= methy_data$start[i] & m_data <= methy_data$end[i]
        upstream <- m_data < methy_data$start[i]
        downstream <- m_data > methy_data$end[i]

        out <- numeric(length(m_data))
        out[within] <- (m_data[within] - methy_data$start[i]) / (methy_data$end[i] - methy_data$start[i])
        out[upstream] <- (m_data[upstream] - methy_data$start[i]) / flank / 3
        out[downstream] <- 1 + ((m_data[downstream] - methy_data$end[i]) / flank / 3)
        if (stranded && length(methy_data$strand[i]) > 0 && methy_data$strand[i] == "-") {
            out <- 1 - out
        }

        methy_data$methy_data[[i]]$rel_pos <- out
    }

    methy_data
}

.pos_avg <- function(x) {
    average_by_pos <- function(x) {
          x %>%
              dplyr::group_by(.data$sample, .data$chr, .data$strand, .data$rel_pos) %>%
              dplyr::summarise(
                  methy_prop = mean(.data$statistic > 0),
                  .groups = "drop")
      }

      # take feature level average
      x <- x %>%
          dplyr::mutate(methy_data = purrr::map(.data$methy_data, average_by_pos))
}

.unnest_anno <- function(x, samples_anno) {
    x %>%
        dplyr::select(!dplyr::any_of(c("chr", "strand", "start", "end"))) %>%
        tidyr::unnest("methy_data") %>%
        dplyr::left_join(samples_anno, by = "sample")
}

.bin_avg <- function(x, group_col, grid_size = 2^10) {
    min <- -1.1/3
    max <- 1 + 1.1/3
    binned_pos <- seq(min, max, length.out = grid_size + 1)
    binned_pos_df <- data.frame(
        binned_pos = binned_pos[-1],
        interval = cut(binned_pos, breaks = binned_pos)[-1]
    )

    x %>%
        dplyr::mutate(interval = cut(.data$rel_pos, breaks = binned_pos)) %>%
        dplyr::left_join(binned_pos_df, by = "interval") %>%
        dplyr::group_by(dplyr::across(dplyr::all_of(group_col)), binned_pos) %>%
        dplyr::summarise(methy_prop = mean(.data$methy_prop), .groups = "drop")
}

.filter_tabix_chr <- function(nmr, regions) {
    seqs <- get_tabix_sequences(paste0(methy(nmr), ".tbi"))
    filter_bad_seqs <- function(x) {
        bad_seq <- !x$chr %in% seqs

        list(
            regions = x[!bad_seq, ],
            bad_seqs = unique(x$chr[bad_seq])
        )
    }
    res <- purrr::map(regions, filter_bad_seqs)

    regions <- purrr::map(res, ~.x$regions)
    bad_seqs <- purrr::map(res, ~.x$bad_seq)
    bad_seqs <- unique(unlist(bad_seqs))
    if (length(bad_seqs) != 0) {
        warning("requested sequences missing from tabix file:",
                paste(bad_seqs, collapse = ", "))
    }

    regions
}
