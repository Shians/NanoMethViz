plot_gene <- function(gene, methy, exons) {
    exons_anno <- query_exons_symbol(exons, symbol = gene)

    p1 <- with(exons_anno,
         plot_agg(
            chr = unique(chr),
            start = min(start),
            end = max(end),
            title = gene,
            methy = methy
        )
    )

    get_ggplot_range_x <- function(x) {
        ggplot_build(x)$layout$panel_scales_x[[1]]$range$range
    }

    xlim <- get_ggplot_range_x(p1)

    p2 <- plot_gene_annotation(exons_anno, xlim[1], xlim[2])

    n_unique <- function(x) { length(unique(x)) }

    p1/p2 + patchwork::plot_layout(height = c(1, 0.075*n_unique(exons_anno$transcript_id)))
}