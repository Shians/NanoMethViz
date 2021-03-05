test_that("Plotting region works", {
    # setup
    nmr <- load_example_nanomethresult()

    # test
    expect_silent(p <- plot_region(nmr, "chr7", 6703892, 6730431))
    expect_true(is(p, "ggplot"))
})
