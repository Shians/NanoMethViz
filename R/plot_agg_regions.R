#' Plot aggregate regions
#'
#' @param x the NanoMethResult object.
#' @param regions a table of regions containing at least columns chr, strand,
#'   start and end. Any additional columns can be used for grouping.
#' @param binary_threshold the modification probability such that calls with
#'   modification probability above the threshold are considered methylated, and
#'   those with probability equal or below are considered unmethylated.
#' @param group_col the column to group aggregated trends by. This column can
#'   be in from the regions table or samples(x).
#' @param flank the number of flanking bases to add to each side of each region.
#' @param stranded TRUE if negative strand features should have coordinates
#'   flipped to reflect features like transcription start sites.
#' @param span the span for loess smoothing.
#' @param palette the ggplot colour palette used for groups.
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
    binary_threshold = 0.5,
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

    # grouped regions crashes downstream operations
    regions <- ungroup(regions)

    # query methylation data
    methy_data <- plot_agg_regions.get_methy_data(x, regions, flank)

    # process each methy_data list element to obtain methylation proportions at relative positions
    methy_data <- plot_agg_regions.process_methy_data(methy_data, stranded, flank, binary_threshold)

    # unnest methy_data and add sample annotation
    methy_data <- plot_agg_regions.unnest_with_anno(methy_data, NanoMethViz::samples(x))

    # take the average methylation proportion over equal sized positional bins
    methy_data <- plot_agg_regions.avg_over_bins(methy_data, group_col = group_col)

    # change aes spec depending on if group_col is available
    # hacky fix to the deprecation of aes_string which handled a NULL group_col
    # value
    if (!is.null(group_col)) {
        aes_spec <- ggplot2::aes(
                x = .data$binned_pos,
                y = .data$methy_prop,
                group = .data[[group_col]],
                col = .data[[group_col]])
    } else {
        aes_spec <- ggplot2::aes(
                x = .data$binned_pos,
                y = .data$methy_prop)
    }

    # set up plot
    p <- ggplot2::ggplot() +
        ggplot2::ylim(c(0, 1)) +
        ggplot2::theme_bw() +
        ggplot2::stat_smooth(
            aes_spec,
            method = "loess",
            formula = "y ~ x",
            span = span,
            na.rm = TRUE,
            se = FALSE,
            data = methy_data) +
        palette

    # if flank is 0, then then only annotated start and end, else annotate flanks as well
    if (flank == 0) {
        labels <- c("start", "end")
        breaks <- c(0, 1)
        limits <- c(0, 1)
    } else {
        breaks <- c(-.33, 0, 1, 1.33)
        limits <- c(-0.33, 1.33)
        kb_marker <- round(flank / 1000, 1)
        labels <- c(glue::glue("-{kb_marker}kb"), "start", "end", glue::glue("+{kb_marker}kb"))
        p <- p +
            ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey80") +
            ggplot2::geom_vline(xintercept = 1, linetype = "dashed", color = "grey80")
    }

    region_widths <- regions$end - regions$start
    width_avg <- round(mean(region_widths))
    width_std_dev <- round(sd(region_widths))

    p + ggplot2::coord_cartesian(clip = "off") +
        ggplot2::scale_x_continuous(
            name = glue::glue("Relative Position (Avg.Width = {width_avg}, Std.Dev = {width_std_dev})"),
            breaks = breaks,
            limits = limits,
            labels = labels) +
        ggplot2::ylab("Average Methylation Proportion")
}

plot_agg_regions.get_methy_data <- function(x, regions, flank) {
    query_row_methy <- function(i, methy, regions, flank) {
        # query a larger region such that smoothing doesn't
        # misbehave around ends
        flank <- flank * 1.1
        query_methy(
            methy,
            regions$chr[i],
            regions$start[i] - flank,
            regions$end[i] + flank,
            force = TRUE)
    }

    methy_data <- purrr::map(
        seq_len(nrow(regions)),
        query_row_methy,
        methy = x,
        regions = regions,
        flank = flank
    )

    regions %>%
        dplyr::mutate(methy_data = methy_data)
}

plot_agg_regions.process_methy_data <- function(methy_data, stranded, flank, binary_threshold) {

    remove_empty_methy_data <- function(x) {
        x %>%
            dplyr::filter(purrr::map_lgl(.data$methy_data, ~nrow(.x) > 0))
    }

    rescale_positions <- function(x, stranded, flank) {
        # for each region, scale the position values to be between 0 and 1
        for (i in seq_len(nrow(x))) {
            # determine if the positions fall inside the region
            m_data <- x$methy_data[[i]]$pos
            within <- m_data >= x$start[i] & m_data <= x$end[i]
            upstream <- m_data < x$start[i]
            downstream <- m_data > x$end[i]

            # rescale the positions depending on the region
            rel_pos <- numeric(length(m_data))
            rel_pos[within] <- (m_data[within] - x$start[i]) / (x$end[i] - x$start[i])

            # flanks each take up 1/3 of the plot space
            rel_pos[upstream] <- (m_data[upstream] - x$start[i]) / flank / 3
            rel_pos[downstream] <- 1 + ((m_data[downstream] - x$end[i]) / flank / 3)

            # flip positions for negative strand if stranded is TRUE
            if (stranded && length(x$strand[i]) > 0 && x$strand[i] == "-") {
                rel_pos <- 1 - rel_pos
            }

            # store rescaled positions
            x$methy_data[[i]]$rel_pos <- rel_pos
        }

        x
    }

    binarise_statistics <- function(x, binary_threshold) {
        binarise_statistc_in_df <- function(x, binary_threshold) {
            x %>%
                dplyr::mutate(is_methylated = e1071::sigmoid(.data$statistic) > binary_threshold)
        }

        x %>%
            dplyr::mutate(
                methy_data = purrr::map(
                    .data$methy_data,
                    binarise_statistc_in_df, binary_threshold = binary_threshold
                )
            )
    }

    average_binarised_methy_per_pos <- function(x) {
        avg_methy_by_rel_pos <- function(x) {
            x %>%
                dplyr::group_by(.data$sample, .data$chr, .data$strand, .data$rel_pos) %>%
                dplyr::summarise(
                    methy_prop = mean(.data$is_methylated),
                    .groups = "drop")
        }

        x %>%
            dplyr::mutate(
                methy_data = purrr::map(.data$methy_data, avg_methy_by_rel_pos)
            )
    }

    methy_data %>%
        remove_empty_methy_data() %>%
        rescale_positions(stranded, flank) %>%
        binarise_statistics(binary_threshold) %>%
        average_binarised_methy_per_pos()
}


plot_agg_regions.unnest_with_anno <- function(x, samples_anno) {
    x %>%
        dplyr::select(!dplyr::any_of(c("chr", "strand", "start", "end"))) %>%
        tidyr::unnest("methy_data") %>%
        dplyr::inner_join(samples_anno, by = "sample", multiple = "all")
}

plot_agg_regions.avg_over_bins <- function(x, group_col, grid_size = 2^10) {
    min <- -1.1 / 3
    max <- 1 + 1.1 / 3
    binned_pos <- seq(min, max, length.out = grid_size + 1)
    binned_pos_df <- data.frame(
        binned_pos = binned_pos[-1],
        interval = cut(binned_pos, breaks = binned_pos)[-1]
    )

    x %>%
        dplyr::mutate(interval = cut(.data$rel_pos, breaks = binned_pos)) %>%
        dplyr::inner_join(binned_pos_df, by = "interval", multiple = "all") %>%
        dplyr::group_by(dplyr::across(dplyr::all_of(group_col)), binned_pos) %>%
        dplyr::summarise(methy_prop = mean(.data$methy_prop), .groups = "drop")
}
