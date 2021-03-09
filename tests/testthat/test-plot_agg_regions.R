test_that("Aggregate plotting works", {
    # setup
    nmr <- load_example_nanomethresult()

    # test
    expect_silent(plot_agg_regions(nmr, NanoMethViz::exons(nmr)))
    expect_silent(plot_agg_regions_sample_grouped(nmr, NanoMethViz::exons(nmr)))
})
