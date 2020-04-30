#' @export
plot_aggregate_regions_by_group <- function(x, features, group_by = NULL) {
    assert_that(is.null(group_by) || is.string(group_by))
    assert_that(features %has_name% group_by)

    methy_data <- .get_features(features, methy = methy(x))

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

#' @export
plot_aggregate_regions_by_sample <- function(x, features) {
    methy_data <- .get_features(features, methy = methy(x))

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
.get_features <- function(features, methy, feature_ids = NULL) {
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
