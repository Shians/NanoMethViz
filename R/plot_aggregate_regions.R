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
        dplyr::summarise(statistic = mean(.data$statistic))

    p <- ggplot2::ggplot(
            ggplot2::aes(
                x = .data$rel_pos,
                y = e1071::sigmoid(.data$statistic)
            ),
            data = methy_data
        ) +
        ggplot2::ylim(c(0, 1)) +
        ggthemes::theme_tufte()

    if (!is.null(groups)) {
        p <- p + ggplot2::stat_smooth(
                aes(group = .data$group, colour = .data$group),
                na.rm = TRUE,
            se = FALSE) +
            ggplot2::scale_colour_brewer(palette = "Set1") +
            ggplot2::coord_cartesian(clip = "off")
    } else {
        p <- p + ggplot2::stat_smooth(se = FALSE, na.rm = TRUE) +
            ggplot2::coord_cartesian(clip = "off")
    }

    kb_marker <- round(flank / 1000, 1)

    labels <- c(glue::glue("-{kb_marker}kb"), "start", "end", glue::glue("+{kb_marker}kb"))

    p +
        ggplot2::scale_x_continuous(
            name = "Relative Position",
            breaks = c(-.25, 0, 1, 1.25),
            labels = labels) +
        ggplot2::ylab("Average Methylation Probability")
}

.split_rows <- function(x) {
    split(x, 1:nrow(x))
}

.row_map <- function(.x, .f, ...) {
    parallel::mclapply(.split_rows(.x), .f, mc.cores = 4L, ...)
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

    out <- query_methy(methy, features$chr, features$start - flank, features$end + flank)
    out <- purrr::map2(out, .split_rows(features), function(x, y) {
        if (nrow(x) == 0) {
            return(tibble::add_column(x, rel_pos = numeric()))
        }

        x$rel_pos <- .scale_to_feature(x$pos, y, stranded = stranded)
        x
    })

    if (!is.null(feature_ids)) {
        assert_that(length(feature_ids) == nrow(out))
        names(out) <- feature_ids
    }

    out
}
