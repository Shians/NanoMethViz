test_that("Plotting gene works", {
    # setup
    nmr <- load_example_nanomethresult()
    mbr <- load_example_modbamresult()
    params <- expand.grid(
        heatmap = c(TRUE, FALSE),
        spaghetti = c(TRUE, FALSE),
        gene_anno = c(TRUE, FALSE)
    )

    # test
    for (x in list(nmr, mbr)) {
        expect_silent(p1 <- plot_gene(x, "Peg3"))
        expect_s3_class(p1, "ggplot")

        expect_silent(p2 <- plot_gene(x, "Peg3", heatmap = TRUE))
        expect_s3_class(p2, "ggplot")

        for (bt in c(0, 0.5, 1)) {
            expect_silent(
                p <- plot_gene(
                    x, "Peg3",
                    binary_threshold = bt
                )
            )
            expect_s3_class(p, "ggplot")
        }
        for (i in 1:nrow(params)) {
            expect_silent(
                p <- plot_gene(
                    x, "Peg3",
                    heatmap = params$heatmap[i],
                    spaghetti = params$spaghetti[i]
                )
            )

            expect_s3_class(p, "ggplot")
        }

        expect_error(plot_gene(x, "missing_gene", heatmap = TRUE))
    }
})
