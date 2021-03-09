test_that("Plotting gene methylation heatmap works", {
    # setup
    nmr <- load_example_nanomethresult()

    # test
    expect_silent(p1 <- plot_gene_heatmap(nmr, "Peg3"))
    expect_true(is(p1, "ggplot"))
})
