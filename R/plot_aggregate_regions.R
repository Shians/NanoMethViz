#' Plot aggregate regions
#'
#' @param x the NanoMethResults object.
#' @param features a table of regions or GRanges, or a list of such objects. The table of regions must contain chr, start and end columns.
#' @param groups if 'features' is a list, a vector of characters of the same length as the list containing names for each member.
#'
#' @return
#' @export
plot_aggregate_regions <- function(x, features, groups = NULL) {
    is_df_or_granges <- function(x) {
        is.data.frame(x) || is(x, "GRanges")
    }

    assert_that(is_df_or_granges(features) || is.list(features))

    if (is(features, "GRanges")) {
        features <- tibble::as_tibble(features) %>%
            dplyr::rename(chr = "seqnames")
    }

    if (!is.null(dim(features))) {
        features <- list(features)
    }

    methy_data <- purrr::map(
        features,
        function(features) {
            .get_features_with_rel_pos(features, methy = methy(x)) %>%
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

    if (!is.null(group)) {
        p <- p + ggplot2::stat_smooth(
                aes(group = .data$group, colour = .data$group),
                method = "loess",
                na.rm = TRUE) +
            ggplot2::scale_colour_brewer(palette = "Set1")
    } else {
        p <- p + ggplot2::stat_smooth(
                method = "loess",
                na.rm = TRUE)
    }

    p +
        ggplot2::xlab("Relative Position") +
        ggplot2::ylab("Average Methylation Probability")
}

plot_aggregate_regions_by_group <- function(x, features, group_by = NULL) {
    assert_that(is.null(group_by) || is.string(group_by))
    assert_that(features %has_name% group_by)

    methy_data <- .get_features_with_rel_pos(features, methy = methy(x))

    if (!is.null(group_by)) {
        names(methy_data) <- features[[group_by]]
    }

    methy_data <- methy_data %>%
        dplyr::bind_rows(.id = "group")

    p <- ggplot2::ggplot(
            ggplot2::aes(
                x = .data$rel_pos,
                y = e1071::sigmoid(.data$statistic)
            ),
            data = methy_data
        ) +
        ggplot2::ylim(c(0, 1)) +
        ggthemes::theme_tufte()

    if (!is.null(group_by)) {
        p + ggplot2::geom_smooth(aes(group = .data$group, colour = .data$group), na.rm = TRUE)
    } else {
        p + ggplot2::geom_smooth(na.rm = TRUE)
    }
}

plot_aggregate_regions_by_sample <- function(x, features) {
    methy_data <- .get_features_with_rel_pos(features, methy = methy(x))

    methy_data <- methy_data %>%
        dplyr::bind_rows() %>%
        dplyr::inner_join(samples(x), by = "sample")

    p <- ggplot2::ggplot(
            ggplot2::aes(
                x = .data$rel_pos,
                y = e1071::sigmoid(.data$statistic)
            ),
            data = methy_data) +
        ggplot2::ylim(c(0, 1)) +
        ggthemes::theme_tufte()

    if (!is.null(group_by)) {
        p + ggplot2::geom_smooth(
            ggplot2::aes(
                group = .data$group, colour = .data$group
            ),
            na.rm = TRUE)
    } else {
        p + ggplot2::geom_smooth(na.rm = TRUE)
    }
}

.split_rows <- function(x) {
    split(x, 1:nrow(x))
}

.row_map <- function(.x, .f, ...) {
    parallel::mclapply(.split_rows(.x), .f, mc.cores = 4L, ...)
}

#' @importFrom purrr map2
.get_features_with_rel_pos <- function(features, methy, feature_ids = NULL) {
    .scale_to_feature <- function(x, feature) {
        (x - feature$start) / (feature$end - feature$start)
    }

    out <- query_methy(methy, features$chr, features$start, features$end)
    out <- purrr::map2(out, .split_rows(features), function(x, y) {
        if (nrow(x) == 0) {
            return(tibble::add_column(x, rel_pos = numeric()))
        }

        x$rel_pos <- .scale_to_feature(x$pos, y)
        x
    })

    if (!is.null(feature_ids)) {
        assert_that(length(feature_ids) == nrow(out))
        names(out) <- feature_ids
    }

    out
}
