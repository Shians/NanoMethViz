plot_aggregate_features <- function(x, features, group_by = NULL) {
    assert_that(is.string(group_by))
    assert_that(features %has_name% group_by)

    .split_rows <- function(x) {
        split(x, 1:nrow(x))
    }

    .row_map <- function(.x, .f, ...) {
        purrr::map(.split_rows(.x), .f, ...)
    }

    .get_feature <- function(feature, methy) {
        .scale_to_feature <- function(x) {
            (x - feature$start) / (feature$end - feature$start)
        }
        query_methy(methy, feature$chr, feature$start, feature$end) %>%
            dplyr::mutate(rel_pos = .scale_to_feature(.data$pos))
    }

    methy_data <- .row_map(features, .get_feature, methy = methy(nmr))

    if (!is.null(group_by)) {
        names(methy_data) <- features[[group_by]]
    }

    methy_data <- methy_data %>%
        dplyr::bind_rows(.id = "group")

    p <-
        ggplot(aes(
            x = rel_pos,
            y = e1071::sigmoid(statistic)
        ), data = methy_data) +
        ylim(c(0, 1)) +
        ggthemes::theme_tufte()

    if (!is.null(group_by)) {
        p + geom_smooth(aes(group = group, colour = group), na.rm = TRUE)
    } else {
        p + geom_smooth(na.rm = TRUE)
    }
}
