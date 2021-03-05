test_that("Plotting gene works", {
    # setup
    nmr <- load_example_nanomethresult()
    p_gene <- plot_gene(nmr, "Peg3")

    # test
    expect_true(is(p_gene, "patchwork"))
    expect_true(is(p_gene, "ggplot"))
})
