test_that("Plotting region heatmap works", {
    # setup
    nmr <- load_example_nanomethresult()

    # test
    expect_silent(p <- plot_region_heatmap(nmr, "chr7", 6703892, 6730431))
    expect_true(is(p, "ggplot"))
})
