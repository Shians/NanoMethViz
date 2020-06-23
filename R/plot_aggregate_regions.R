#' Plot aggregate regions
#'
#' @param x the NanoMethResults object.
#' @param regions a table of regions or GRanges, or a list of such objects. The table of regions must contain chr, start and end columns.
#' @param groups if 'features' is a list, a vector of characters of the same length as the list containing names for each member.
#' @param flank the number of flanking bases to add to each side of each region.
#'
#' @return
#' @export
plot_aggregate_regions <- function(x, regions, groups = NULL, flank = 2000, stranded = TRUE) {
    is_df_or_granges <- function(x) {
        is.data.frame(x) || is(x, "GRanges")
    }

    assert_that(is_df_or_granges(regions) || is.list(regions))

    if (is(regions, "GRanges")) {
        regions <- tibble::as_tibble(regions) %>%
            dplyr::rename(chr = "seqnames")
    }

    if (!is.null(dim(regions))) {
        regions <- list(regions)
    }

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

    methy_data <- methy_data %>%
        dplyr::group_by(.data$group, .data$rel_pos) %>%
        dplyr::summarise(statistic = mean(.data$statistic)) %>%
        dplyr::mutate(methy_prob = e1071::sigmoid(statistic)) %>%
        dplyr::ungroup()

    methy_data_before <- filter(methy_data, rel_pos <= 0)
    if (nrow(methy_data_before) > 10000) {
        methy_data_before <- dplyr::sample_n(methy_data_before, 10000)
    }

    methy_data_after <- filter(methy_data, rel_pos >= 1)
    if (nrow(methy_data_after) > 10000) {
        methy_data_after <- dplyr::sample_n(methy_data_after, 10000)
    }

    methy_data <- filter(methy_data, between(rel_pos, 0, 1))
    if (nrow(methy_data) > 10000) {
        methy_data <- dplyr::sample_n(methy_data, 10000)
    }

    p <- ggplot2::ggplot() +
        ggplot2::ylim(c(0, 1)) +
        ggplot2::theme_minimal()

    if (!is.null(groups)) {
        p <- p +
            .geom_smooth(methy_data, span = 0.3, group = TRUE) +
            .geom_smooth(methy_data_before, span = 0.35, group = TRUE) +
            .geom_smooth(methy_data_after, span = 0.35, group = TRUE) +
            ggplot2::scale_colour_brewer(palette = "Set1") +
            ggplot2::coord_cartesian(clip = "off")
    } else {
        p <- p +
            .geom_smooth(methy_data, span = 0.3, group = FALSE) +
            .geom_smooth(methy_data_before, span = 0.35, group = FALSE) +
            .geom_smooth(methy_data_after, span = 0.35, group = FALSE) +
            ggplot2::coord_cartesian(clip = "off")
    }

    kb_marker <- round(flank / 1000, 1)

    labels <- c(glue::glue("-{kb_marker}kb"), "start", "end", glue::glue("+{kb_marker}kb"))

    p +
        geom_vline(xintercept = 0, linetype = "dashed", color = "grey80") +
        geom_vline(xintercept = 1, linetype = "dashed", color = "grey80") +
        ggplot2::scale_x_continuous(
            name = "Relative Position",
            breaks = c(-.25, 0, 1, 1.25),
            limits = c(-0.25, 1.25),
            labels = labels) +
        ggplot2::ylab("Average Methylation Probability")
}

.split_rows <- function(x) {
    split(x, 1:nrow(x))
}

#' @importFrom purrr map2
.get_features_with_rel_pos <- function(features, methy, flank, stranded, feature_ids = NULL) {
    .scale_to_feature <- function(x, feature, stranded = stranded) {
        within <- x >= feature$start & x <= feature$end
        upstream <- x < feature$start
        downstream <- x > feature$end

        out <- numeric(length(x))
        out[within] <- (x[within] - feature$start) / (feature$end - feature$start)
        out[upstream] <- (x[upstream] - feature$start) / flank / 4
        out[downstream] <- 1 + ((x[downstream] - feature$end) / flank / 4)
        if (stranded && length(feature$strand) > 0 && feature$strand == "-") {
            out <- 1 - out
        }
        out
    }

    out <- query_methy(
        methy,
        features$chr,
        features$start - flank * 1.1,
        features$end + flank * 1.1,
        simplify = FALSE)

    out <- purrr::map2(out, .split_rows(features), function(x, y) {
        if (nrow(x) == 0) {
            return(tibble::add_column(x, rel_pos = numeric()))
        }

        x$rel_pos <- .scale_to_feature(x$pos, y, stranded = stranded)

        x %>%
            group_by(chr, rel_pos) %>%
            summarise(statistic = mean(statistic)) %>%
            ungroup()
    })

    if (!is.null(feature_ids)) {
        assert_that(length(feature_ids) == nrow(out))
        names(out) <- feature_ids
    }

    out
}

.geom_smooth <- function(data, span, group) {
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
