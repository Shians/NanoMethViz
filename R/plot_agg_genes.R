#' Plot gene aggregate plot
#'
#' @param genes a character vector of genes to include in aggregate plot, if NULL then all genes are used.
#' @inheritParams plot_agg_regions
#'
#' @return a ggplot object containing the aggregate methylation trend of genes.
#'
#' @export
#'
#' @examples
#' nmr <- load_example_nanomethresult()
#' plot_agg_genes(nmr)
#'
plot_agg_genes <- function(
    x,
    genes = NULL,
    binary_threshold = 0.5,
    group_col = NULL,
    flank = 2000,
    stranded = TRUE,
    span = 0.05,
    palette = ggplot2::scale_colour_brewer(palette = "Set1")
) {
    gene_regions <- exons_to_genes(exons(x))
    if (!is.null(genes)) {
        gene_regions <- gene_regions %>%
            filter(.data$symbol %in% genes)
    }

    plot_agg_regions(
        x,
        regions = gene_regions,
        binary_threshold = binary_threshold,
        group_col = group_col,
        flank = flank,
        stranded = stranded,
        span = span,
        palette = palette
    )
}

plot_agg_tss <- function(
    x,
    genes = NULL,
    binary_threshold = 0.5,
    group_col = NULL,
    flank = 2000,
    stranded = TRUE,
    span = 0.05,
    palette = ggplot2::scale_colour_brewer(palette = "Set1")
) {
    gene_regions <- exons_to_genes(exons(x))
    if (!is.null(genes)) {
        gene_regions <- gene_regions %>%
            filter(.data$symbol %in% genes)
    }
    tss_regions <- gene_regions %>%
        mutate(
            start = case_when(
                strand == "+" ~ .data$start,
                strand == "-" ~ .data$end,
                TRUE ~ .data$start)
        ) %>%
        mutate(
            start = .data$start - flank,
            end = .data$start + 2 * flank
        )
    kb_marker <- round(flank / 1000, 1)
    labels <- c(
        glue::glue("-{kb_marker}kb"),
        "TSS",
        glue::glue("+{kb_marker}kb")
    )

    p <- plot_agg_regions(
        x,
        regions = tss_regions,
        binary_threshold = binary_threshold,
        group_col = group_col,
        flank = 0,
        stranded = stranded,
        span = span,
        palette = palette
    )
    # hack to delete existing to avoid warning
    p$scales$scales[[which(p$scales$find("x"))]] <- NULL

    p + ggplot2::scale_x_continuous(
        name = glue::glue(""),
        breaks = c(0, 0.5, 1),
        limits = c(0, 1),
        labels = labels
    )
}

plot_agg_tes <- function(
    x,
    genes = NULL,
    binary_threshold = 0.5,
    group_col = NULL,
    flank = 2000,
    stranded = TRUE,
    span = 0.05,
    palette = ggplot2::scale_colour_brewer(palette = "Set1")
) {
    gene_regions <- exons_to_genes(exons(x))
    if (!is.null(genes)) {
        gene_regions <- gene_regions %>%
            filter(.data$symbol %in% genes)
    }
    tes_regions <- gene_regions %>%
        mutate(
            start = case_when(
                strand == "+" ~ .data$end,
                strand == "-" ~ .data$start,
                TRUE ~ .data$end)
        ) %>%
        mutate(
            start = .data$start - flank,
            end = .data$start + 2 * flank
        )
    kb_marker <- round(flank / 1000, 1)
    labels <- c(
        glue::glue("-{kb_marker}kb"),
        "TES",
        glue::glue("+{kb_marker}kb")
    )

    p <- plot_agg_regions(
        x,
        regions = tes_regions,
        binary_threshold = binary_threshold,
        group_col = group_col,
        flank = 0,
        stranded = stranded,
        span = span,
        palette = palette
    )
    # hack to delete existing to avoid warning
    p$scales$scales[[which(p$scales$find("x"))]] <- NULL

    p + ggplot2::scale_x_continuous(
        name = glue::glue(""),
        breaks = c(0, 0.5, 1),
        limits = c(0, 1),
        labels = labels
    )
}
