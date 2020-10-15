#' Plot aggregate regions
#'
#' @param x the NanoMethResult object.
#' @param regions a table of regions or GRanges, or a list of such objects. The
#'   table of regions must contain chr, start and end columns.
#' @param groups if 'features' is a list, a vector of characters of the same
#'   length as the list containing names for each member.
#' @param flank the number of flanking bases to add to each side of each region.
#' @param stranded TRUE if negative strand features should have coordinates
#'   flipped to reflect features like transcription start sites.
#' @param span the span for loess smoothing.
#'
#' @return a ggplot object.
#'
#' @export
plot_agg_regions <- function(
    x,
    regions,
    groups = NULL,
    flank = 2000,
    stranded = TRUE,
    span = 0.05
) {
    is_df_or_granges <- function(x) {
        is.data.frame(x) || is(x, "GRanges")
    }

    assert_that(is_df_or_granges(regions) || is.list(regions))

    # convert GRanges to tibble
    if (is(regions, "GRanges")) {
        regions <- tibble::as_tibble(regions) %>%
            dplyr::rename(chr = "seqnames")
    }

    # convert data.frame regions to single element list
    if (!is.null(dim(regions))) {
        regions <- list(regions)
    }

    # query methylation data
    methy_data <- purrr::map(
        regions,
        function(features) {
            .get_features_with_rel_pos(
                    features,
                    methy = methy(x),
                    flank = flank,
                    stranded = stranded) %>%
                dplyr::bind_rows()
        }
    )

    if (!is.null(groups)) {
        names(methy_data) <- groups
    }

    methy_data <- methy_data %>%
        dplyr::bind_rows(.id = "group")

    if (!is.null(groups)) {
        methy_data <- methy_data %>%
        dplyr::group_by(.data$group, .data$rel_pos) %>%
        dplyr::summarise(methy_prob = mean(.data$methy_prob)) %>%
        dplyr::ungroup()
    } else {
        methy_data <- methy_data %>%
            dplyr::group_by(.data$rel_pos) %>%
            dplyr::summarise(methy_prob = mean(.data$methy_prob)) %>%
            dplyr::ungroup()
    }

    # set up plot
    p <- ggplot2::ggplot() +
        ggplot2::ylim(c(0, 1)) +
        ggplot2::theme_minimal()

    # take binned means
    grid_size <- 2^12
    binned_pos <- seq(-1.1/3, 1 + 1.1/3, length.out = grid_size+ 2)[2:(1 + grid_size)]
    methy_data <- methy_data %>%
        mutate(interval = cut(.data$rel_pos, breaks = grid_size)) %>%
        group_by(.data$interval, .drop = TRUE) %>%
        summarise(methy_prob = mean(.data$methy_prob)) %>%
        ungroup() %>%
        mutate(rel_pos = .data$binned_pos)

    if (!is.null(groups)) {
        p <- p +
            .agg_geom_smooth(methy_data, span = span, group = TRUE) +
            ggplot2::scale_colour_brewer(palette = "Set1") +
            ggplot2::coord_cartesian(clip = "off")
    } else {
        p <- p +
            .agg_geom_smooth(methy_data, span = span, group = FALSE) +
            ggplot2::coord_cartesian(clip = "off")
    }

    kb_marker <- round(flank / 1000, 1)

    labels <- c(glue::glue("-{kb_marker}kb"), "start", "end", glue::glue("+{kb_marker}kb"))

    p +
        geom_vline(xintercept = 0, linetype = "dashed", color = "grey80") +
        geom_vline(xintercept = 1, linetype = "dashed", color = "grey80") +
        ggplot2::scale_x_continuous(
            name = "Relative Position",
            breaks = c(-.33, 0, 1, 1.33),
            limits = c(-0.33, 1.33),
            labels = labels) +
        ggplot2::ylab("Average Methylation Probability")
}

#' Plot aggregate regions with grouped samples
#'
#' @inheritParams plot_agg_regions
#'
#' @return a ggplot plot object.
#'
#' @export
plot_agg_regions_sample_grouped <- function(
    x,
    regions,
    flank = 2000,
    stranded = TRUE,
    span = 0.05
) {
    is_df_or_granges <- function(x) {
        is.data.frame(x) || is(x, "GRanges")
    }

    assert_that(is_df_or_granges(regions) || is.list(regions))

    # convert GRanges to tibble
    if (is(regions, "GRanges")) {
        regions <- tibble::as_tibble(regions) %>%
            dplyr::rename(chr = "seqnames")
    }

    # convert data.frame regions to single element list
    if (!is.null(dim(regions))) {
        regions <- list(regions)
    }

    # query methylation data
    methy_data <- purrr::map(
        regions,
        function(features) {
            .get_features_with_rel_pos(
                    features,
                    methy = methy(x),
                    flank = flank,
                    stranded = stranded,
                    by_sample = TRUE) %>%
                dplyr::bind_rows()
        }
    )

    methy_data <- methy_data %>%
        dplyr::bind_rows()

    methy_data <- methy_data %>%
        dplyr::group_by(sample, .data$rel_pos) %>%
        dplyr::summarise(methy_prob = mean(.data$methy_prob)) %>%
        dplyr::ungroup()

    # set up plot
    p <- ggplot2::ggplot() +
        ggplot2::ylim(c(0, 1)) +
        ggplot2::theme_minimal()

    # take binned means
    grid_size <- 2^12
    binned_pos_df <- tibble::tibble(
        interval = levels(cut(methy_data$rel_pos, breaks = grid_size)),
        binned_pos = seq(-1.1/3, 1 + 1.1/3, length.out = grid_size+ 2)[2:(1 + grid_size)]
    )

    methy_data <- methy_data %>%
        dplyr::mutate(interval = cut(.data$rel_pos, breaks = grid_size)) %>%
        dplyr::group_by(.data$sample, .data$interval) %>%
        dplyr::summarise(methy_prob = mean(.data$methy_prob)) %>%
        dplyr::ungroup() %>%
        dplyr::left_join(binned_pos_df, by = "interval") %>%
        dplyr::mutate(rel_pos = .data$binned_pos) %>%
        dplyr::left_join(samples(x), by = "sample")

    p <- p +
        ggplot2::stat_smooth(
            aes(x = .data$rel_pos, y = .data$methy_prob, group = .data$haplotype, col = .data$haplotype),
            method = "loess",
            span = span,
            na.rm = TRUE,
            se = FALSE,
            data = methy_data) +
        ggplot2::scale_colour_brewer(palette = "Set1") +
        ggplot2::coord_cartesian(clip = "off")

    kb_marker <- round(flank / 1000, 1)

    labels <- c(glue::glue("-{kb_marker}kb"), "start", "end", glue::glue("+{kb_marker}kb"))

    p +
        geom_vline(xintercept = 0, linetype = "dashed", color = "grey80") +
        geom_vline(xintercept = 1, linetype = "dashed", color = "grey80") +
        ggplot2::scale_x_continuous(
            name = "Relative Position",
            breaks = c(-.33, 0, 1, 1.33),
            limits = c(-0.33, 1.33),
            labels = labels) +
        ggplot2::ylab("Average Methylation Probability")
}

.split_rows <- function(x) {
    split(x, seq_len(nrow(x)))
}

.get_features_with_rel_pos <- function(features, methy, flank, stranded, by_sample = FALSE, feature_ids = NULL) {
    .scale_to_feature <- function(x, feature, stranded = stranded) {
        within <- x >= feature$start & x <= feature$end
        upstream <- x < feature$start
        downstream <- x > feature$end

        out <- numeric(length(x))
        out[within] <- (x[within] - feature$start) / (feature$end - feature$start)
        out[upstream] <- (x[upstream] - feature$start) / flank / 3
        out[downstream] <- 1 + ((x[downstream] - feature$end) / flank / 3)
        if (stranded && length(feature$strand) > 0 && feature$strand == "-") {
            out <- 1 - out
        }
        out
    }

    out <- query_methy(
        methy,
        features$chr,
        features$start - flank * 1.10,
        features$end + flank * 1.10,
        simplify = FALSE)

    empty <- purrr::map_lgl(out, function(x) nrow(x) == 0)
    out <- out[!empty]
    features <- features[!empty, ]

    out <- purrr::map2(out, .split_rows(features), function(x, y) {
        if (nrow(x) == 0) {
            return(tibble::add_column(x, rel_pos = numeric()))
        }

        x$rel_pos <- .scale_to_feature(x$pos, y, stranded = stranded)

        if (by_sample) {
            x %>%
                dplyr::group_by(.data$sample, .data$rel_pos) %>%
                dplyr::summarise(methy_prob = mean(e1071::sigmoid(.data$statistic))) %>%
                dplyr::ungroup()
        } else {
            x %>%
                dplyr::group_by(.data$rel_pos) %>%
                dplyr::summarise(methy_prob = mean(e1071::sigmoid(.data$statistic))) %>%
                dplyr::ungroup()
        }
    })

    if (!is.null(feature_ids)) {
        assertthat::assert_that(length(feature_ids) == nrow(out))
        names(out) <- feature_ids
    }

    out
}

.agg_geom_smooth <- function(data, span, group) {
    if (group) {
        ggplot2::stat_smooth(
            aes(x = .data$rel_pos,
                y = .data$methy_prob,
                group = .data$group,
                col = .data$group),
            method = "loess",
            span = span,
            na.rm = TRUE,
            se = FALSE,
            data = data)
    } else {
        ggplot2::stat_smooth(
            aes(x = .data$rel_pos, y = .data$methy_prob),
            method = "loess",
            span = span,
            na.rm = TRUE,
            se = FALSE,
            data = data)
    }
}
