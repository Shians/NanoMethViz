#' @importFrom stats median
#' @importFrom forcats fct_reorder
lollipop <- function(x, chr, start, end, binary = FALSE) {
    methy_data <- query_methy(x, chr, start, end)

    numeric_median <- function(x) stats::median(as.numeric(x))

    df <- methy_data %>%
        dplyr::group_by(.data$sample, .data$pos) %>%
        dplyr::summarise(val = mean(e1071::sigmoid(.data$statistic))) %>%
        dplyr::mutate(
            pos = as.factor(scales::comma(.data$pos, accuracy = 1)),
            group = str_extract(.data$sample, "bl6|cast")
        ) %>%
        dplyr::ungroup() %>%
        dplyr::mutate(
            group = factor(.data$group),
            sample = forcats::fct_reorder(.data$sample, .data$group, numeric_median),
            state = ifelse(.data$val > 0.5, "methylated", "unmethylated")
        )

    p <- ggplot2::ggplot(df, aes(x = .data$pos, y = .data$sample)) +
        ggplot2::geom_hline(yintercept = df$sample, size = 2) +
        ggplot2::theme(
            axis.text.x = ggplot2::element_blank(),
            axis.ticks.x = ggplot2::element_blank(),
            axis.title.x = ggplot2::element_blank())

    if (binary) {
        p <- p +
            ggplot2::geom_point(aes(fill = .data$state), pch = 21, size = 7) +
            ggplot2::scale_fill_manual(values = c("black", "white"))

    } else {
        p <- p +
            ggplot2::geom_point(aes(fill = .data$val), pch = 21, size = 7) +
            ggplot2::scale_fill_steps(low = "white", high = "blueviolet", n.breaks = 6)
    }

    p
}
