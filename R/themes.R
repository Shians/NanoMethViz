 # heatmap
theme_methy_heatmap <- ggplot2::theme_bw() +
        ggplot2::theme(
            axis.ticks.y = ggplot2::element_blank(),
            axis.title.y = ggplot2::element_blank(),
            axis.text.y = ggplot2::element_blank(),
            panel.grid.major.y = ggplot2::element_blank(),
            panel.grid.minor.y = ggplot2::element_blank()
        )

heatmap_fill_scale <- scico::scale_fill_scico(palette = 'imola', direction = -1)
heatmap_col_scale <- scico::scale_colour_scico(palette = 'imola', direction = -1)

# Future releast candidate
# heatmap_fill_scale <- scico::scale_fill_scico(palette = 'roma', direction = -1)
# heatmap_col_scale <- scico::scale_colour_scico(palette = 'roma', direction = -1)
