lollipop <- function(x, chr, start, end, binary = FALSE) {
    methy_data <- query_methy(x, chr, start, end)

    df <- methy_data %>%
        group_by(sample, pos) %>%
        summarise(val = mean(e1071::sigmoid(statistic))) %>%
        mutate(
            pos = as.factor(scales::comma(pos, accuracy = 1)),
            group = str_extract(sample, "bl6|cast")
        ) %>%
        ungroup() %>%
        mutate(
            group = factor(group),
            sample = fct_reorder(sample, group, function(x) median(as.numeric(x))),
            state = ifelse(val > 0.5, "methylated", "unmethylated")
        )

    p <- ggplot(df, aes(x = pos, y = sample)) +
        geom_hline(yintercept = df$sample, size = 2) +
        theme(
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.title.x = element_blank())

    if (binary) {
        p <- p +
            geom_point(aes(fill = state), pch = 21, size = 7) +
            scale_fill_manual(values = c("black", "white"))

    } else {
        p <- p +
            geom_point(aes(fill = val), pch = 21, size = 7) +
            scale_fill_steps(low = "white", high = "blueviolet", n.breaks = 6)
    }

    p
}
