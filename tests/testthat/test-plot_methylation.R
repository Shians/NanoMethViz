test_that("Plotting gene works", {
    # setup
    nmr <- load_example_nanomethresult()

    # test
    p_gene <- expect_no_warning(plot_gene(nmr, "Peg3"))
    p_gene2 <- expect_no_warning(plot_gene(nmr, "Peg3", spaghetti = TRUE))
    expect_s3_class(p_gene, "patchwork")
    expect_s3_class(p_gene, "ggplot")

    expect_s3_class(p_gene2, "patchwork")
    expect_s3_class(p_gene2, "ggplot")
})

test_that("deprecated 'span' argument warns", {
    # setup
    nmr <- load_example_nanomethresult()
    span_warning <- "the 'span' argument has been deprecated, please use 'smoothing_window' instead"
    methy_data <- query_methy(nmr, "chr7", 6703892, 6730431) %>%
        dplyr::select(-"strand")

    # test
    expect_warning(
        plot_grange(nmr, GenomicRanges::GRanges("chr7:6703892-6730431"), heatmap = FALSE, span = 0.5),
        span_warning,
        fixed = TRUE
    )

    # plot_region() accepts 'span' but does not forward it, so call the
    # internal plotting function directly
    expect_warning(
        plot_methylation_data(
            methy_data = methy_data,
            sample_anno = samples(nmr),
            chr = "chr7",
            start = 6703892,
            end = 6730431,
            title = "Peg3",
            span = 0.1
        ),
        span_warning,
        fixed = TRUE
    )
})
