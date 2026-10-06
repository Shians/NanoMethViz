setup_violin_data <- function() {
    nmr <- load_example_nanomethresult()
    regions <- exons_to_genes(NanoMethViz::exons(nmr))
    list(nmr = nmr, regions = regions)
}

test_that("plot_violin errors on unknown group_col", {
    d <- setup_violin_data()

    expect_error(
        plot_violin(d$nmr, d$regions, group_col = "foo"),
        "'foo' could not be found"
    )
})

test_that("plot_violin default call builds a ggplot", {
    d <- setup_violin_data()

    p <- expect_no_warning(plot_violin(d$nmr, d$regions))
    expect_s3_class(p, "ggplot")
    expect_no_warning(ggplot2::ggplot_build(p))

    expect_true(all(p$data$methy_prop >= 0 & p$data$methy_prop <= 1))
    # one value per region per sample
    expect_equal(
        nrow(p$data),
        nrow(dplyr::distinct(p$data, .data$gene_id, .data$sample))
    )
    expect_setequal(p$data$group, samples(d$nmr)$group)
})

test_that("plot_violin groups by sample", {
    d <- setup_violin_data()

    p <- plot_violin(d$nmr, d$regions, group_col = "sample")
    built <- ggplot2::ggplot_build(p)

    expect_equal(
        length(unique(built$data[[1]]$x)),
        length(unique(samples(d$nmr)$sample))
    )
})

test_that("plot_violin draws a single violin when group_col is NULL", {
    d <- setup_violin_data()

    p <- plot_violin(d$nmr, d$regions, group_col = NULL)
    built <- expect_no_error(ggplot2::ggplot_build(p))

    expect_equal(length(unique(built$data[[1]]$x)), 1)
})

test_that("plot_violin groups by a column of the regions table", {
    d <- setup_violin_data()
    regions <- d$regions
    regions$type <- rep(c("a", "b"), length.out = nrow(regions))

    p <- plot_violin(d$nmr, regions, group_col = "type")
    built <- expect_no_error(ggplot2::ggplot_build(p))

    expect_setequal(p$data$type, c("a", "b"))
    expect_equal(length(unique(built$data[[1]]$x)), 2)
})

test_that("plot_violin works with only coordinate columns", {
    d <- setup_violin_data()
    regions <- d$regions[, c("chr", "strand", "start", "end")]

    p <- plot_violin(d$nmr, regions)
    expect_no_error(ggplot2::ggplot_build(p))
})

test_that("plot_violin accepts grouped regions", {
    d <- setup_violin_data()
    grouped <- dplyr::group_by(d$regions, .data$chr)

    p <- plot_violin(d$nmr, grouped)
    expect_no_error(ggplot2::ggplot_build(p))
    expect_equal(p$data, plot_violin(d$nmr, d$regions)$data)
})

test_that("plot_violin drops regions without data", {
    d <- setup_violin_data()
    empty_region <- tibble::tibble(
        gene_id = "empty",
        chr = "chr11",
        strand = "+",
        symbol = "Empty",
        start = 1L,
        end = 100L
    )
    regions <- dplyr::bind_rows(d$regions, empty_region)

    p <- plot_violin(d$nmr, regions)
    expect_no_error(ggplot2::ggplot_build(p))

    expect_false("empty" %in% p$data$gene_id)
    expect_setequal(unique(p$data$gene_id), d$regions$gene_id)
})

test_that("plot_violin binary_threshold changes proportions", {
    d <- setup_violin_data()

    low <- plot_violin(d$nmr, d$regions, binary_threshold = 0.1)$data
    high <- plot_violin(d$nmr, d$regions, binary_threshold = 0.9)$data

    expect_true(all(low$methy_prop >= high$methy_prop))
    expect_gt(mean(low$methy_prop), mean(high$methy_prop))
})
