test_that("Aggregate plotting works", {
    # setup
    nmr <- load_example_nanomethresult()

    # test
    expect_silent(plot_agg_regions(nmr, NanoMethViz::exons(nmr)))
    expect_silent(plot_agg_regions_sample_grouped(nmr, NanoMethViz::exons(nmr)))
})

test_that("Aggregate plotting error checking works", {
    # setup
    nmr <- load_example_nanomethresult()
    missing_chr <- NanoMethViz::exons(nmr)[, -2]
    missing_chr_list <- list(
        NanoMethViz::exons(nmr)[, -2],
        NanoMethViz::exons(nmr),
        NanoMethViz::exons(nmr)[, -2]
    )

    # test
    expect_error(plot_agg_regions(nmr, missing_chr))
    expect_error(plot_agg_regions(nmr, missing_chr_list))
    expect_error(plot_agg_regions_sample_grouped(nmr, missing_chr))
    expect_error(plot_agg_regions_sample_grouped(nmr, missing_chr_list))
})
