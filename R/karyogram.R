karyogram <- function(nmr, chr_info) {
    assert_has_columns(chr_info, c("chr", "length"))

    n_chr <- nrow(chr_info)
    plots <- list()

    widest <- max(chr_info$length)

    for (i in 1:n_chr) {
        m <- query_methy(nmr, chr_info$chr[i], 1, chr_info$length[i])

        bins <- round(seq(1, chr_info$length[i], length.out = 2^10))

        bin_means <- binMeans(m$statistic > 0, x = m$pos, bx = bins)

        df <- tibble(
            pos = bins[-length(bins)],
            methy_prop = as.numeric(bin_means),
            count = attr(bin_means, "count")
        )

        plots[[i]] <- ggplot(df, aes(x = pos, y = 1, fill = methy_prop)) +
            geom_tile() +
            scico::scale_colour_scico(palette = 'imola') +
            theme_void() +
            theme(plot.margin = margin(0, 2, 0, 0)) +
            ggtitle(chr_info$chr[i]) +
            xlim(1, widest)
    }

    patchwork::wrap_plots(plots, ncol = 1, guides = "collect")
}

