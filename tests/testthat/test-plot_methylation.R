test_that("Plotting gene works", {
    # setup
    nmr <- load_example_nanomethresult()

    # test
    expect_silent(p_gene <- plot_gene(nmr, "Peg3"))
    expect_silent(p_gene2 <- plot_gene(nmr, "Peg3", spaghetti = TRUE))
    expect_true(is(p_gene, "patchwork"))
    expect_true(is(p_gene, "ggplot"))

    expect_true(is(p_gene2, "patchwork"))
    expect_true(is(p_gene2, "ggplot"))
})
