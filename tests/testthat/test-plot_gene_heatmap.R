test_that("Plotting gene methylation heatmap works", {
    # setup
    nmr <- load_example_nanomethresult()

    # test
    expect_silent(p1 <- plot_gene_heatmap(nmr, "Peg3"))
    expect_true(is(p1, "ggplot"))

    # test for a bug whereby samples not present in data cause function to hang
    nmr_extra_sample <- load_example_nanomethresult()
    samples(nmr_extra_sample) <- bind_rows(
        samples(nmr_extra_sample),
        c(sample = "foo", group = "bar")
    )
    expect_silent(p <- plot_gene_heatmap(nmr_extra_sample, "Peg3"))
    expect_true(is(p, "ggplot"))
})
