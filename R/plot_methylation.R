plot_methy_data <- function(methy_data, chr, start, end, title, spaghetti = TRUE, prop = TRUE) {
    if (prop) {
        plot_data_prop <- methy_data %>%
            mutate(
                haplotype = stringr::str_extract(sample, "cast|bl6"),
                mod_prob = e1071::sigmoid(statistic)
            ) %>%
            group_by(
                haplotype, pos
            ) %>%
            summarise(
                mod_prop = sum(mod_prob > 0.5) / n()
            )
    }

    plot_data <- methy_data %>%
        mutate(
            haplotype = stringr::str_extract(sample, "cast|bl6"),
            mod_prob = e1071::sigmoid(statistic)
        )

    cast_count <- sum(plot_data$haplotype == "cast")
    bl6_count <- sum(plot_data$haplotype == "bl6")
    total_count <- cast_count + bl6_count
    hap_ratio <- cast_count / total_count

    if (prop) {
        p <- ggplot(plot_data, aes(x = pos, col = haplotype)) +
            stat_smooth(data = plot_data_prop, aes(y = mod_prop), geom = "line", method = "loess", size = 4, span = 0.4, formula = y~x)
    } else {
        p <- ggplot(plot_data, aes(x = pos, col = haplotype)) +
            stat_smooth(aes(y = mod_prob), geom = "line", method = "loess", size = 4, span = 0.5, formula = y~x)
    }
    x_min <- max(min(plot_data$pos), start)
    x_max <- min(max(plot_data$pos), end)

    if (spaghetti) {
        p <- p + stat_smooth(aes(y = mod_prob, group = read_name), alpha = 0.25, geom = "line", method = "loess", se = FALSE, span = 1)
    }

    p +
        geom_rug(aes(col = NULL), sides = "b") +
        ggtitle(title) +
        ylim(0, 1) +
        annotate("rect", xmin = start, xmax = end, ymin = 0, ymax = 1, alpha = .1) +
        # xlab() +
        scale_x_continuous(breaks = c(x_min, x_max), label = scales::comma(c(x_min, x_max))) +
        scale_color_brewer(palette = "Set1") +
        scale_fill_brewer(palette = "Set1") +
        ggthemes::theme_tufte()
}

plot_agg <- function(chr, start, end, title, methy, window_prop = 0.3, prop = TRUE) {
    if (is(methy, "character")) {
        con <- DBI::dbConnect(RSQLite::SQLite(), dbname = methy)
        methy <- tbl(con, "methylation")
        need_close <- TRUE
    } else {
        stopifnot(is(methy, "tbl_dbi"))
        need_close <- FALSE
    }

    window <- (end - start) * window_prop

    methy_data <- methy %>%
        dplyr::filter(chr == !!chr) %>%
        dplyr::filter(between(pos, !!start - window, !!end + window)) %>%
        dplyr::select(-strand, -modified) %>%
        as_tibble()

    if (need_close) {
        DBI::dbDisconnect(con)
    }

    plot_methy_data(
        methy_data = methy_data,
        start = start,
        end = end,
        chr = chr,
        title = title,
        spaghetti = FALSE,
        prop = prop
    )
}

plot_spaghetti_and_agg <- function(chr, start, end, title, methy, window_prop = 0.3) {
    if (is(methy, "character")) {
        con <- DBI::dbConnect(RSQLite::SQLite(), dbname = methy)
        methy <- tbl(con, "methylation")
        need_close <- TRUE
    } else {
        stopifnot(is(methy, "tbl_dbi"))
        need_close <- FALSE
    }

    window <- (end - start) * window_prop

    methy_data <- methy %>%
        dplyr::filter(chr == !!chr) %>%
        dplyr::filter(between(pos, !!start - window, !!end + window)) %>%
        dplyr::select(-strand, -modified) %>%
        as_tibble()

    if (need_close) {
        DBI::dbDisconnect(con)
    }

    plot_methy_data(
        methy_data = methy_data,
        start = start,
        end = end,
        chr = chr,
        title = title
    )
}

plot_feature <- function(feature, title = "", methy, window_prop= 0.3, spaghetti = TRUE) {
    if (spaghetti) {
        plot_spaghetti_and_agg(
            feature$chr,
            feature$start,
            feature$end,
            title = title,
            methy = methy,
            window_prop = window_prop
        )
    } else {
        plot_agg(
            feature$chr,
            feature$start,
            feature$end,
            title = title,
            methy = methy,
            window_prop = window_prop
        )
    }
}
