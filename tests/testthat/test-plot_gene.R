test_that("Plotting gene works", {
    # setup
    nmr <- load_example_nanomethresult()

    # test
    expect_silent(p <- plot_gene(nmr, "Peg3"))
    expect_true(is(p, "ggplot"))

    params <- expand.grid(
        heatmap = c(TRUE, FALSE),
        spaghetti = c(TRUE, FALSE),
        gene_anno = c(TRUE, FALSE)
    )

    for (i in 1:nrow(params)) {
        expect_silent(
            plot_gene(
                nmr,
                "Peg3",
                heatmap = params$heatmap[i],
                spaghetti = params$spaghetti[i]
            )
        )
    }
})
