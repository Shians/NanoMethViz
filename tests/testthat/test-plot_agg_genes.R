test_that("Plotting gene aggregates works", {
    nmr <- load_example_nanomethresult()
    expect_silent(plot_agg_genes(nmr))
    expect_silent(plot_agg_tss(nmr))
    expect_silent(plot_agg_tes(nmr))
})
