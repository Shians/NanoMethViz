test_that("Plotting region works", {
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
        expect_silent(p1 <- plot_region(x, "chr7", 6703892, 6730431))
        expect_s3_class(p1, "ggplot")

        expect_silent(p2 <- plot_region(x, "chr7", 6703892, 6730431, heatmap = TRUE))
        expect_s3_class(p2, "ggplot")

        expect_silent(p3 <- plot_region(x, factor("chr7"), 6703892, 6730431))
        expect_s3_class(p3, "ggplot")

        expect_silent(p4 <- plot_region(x, factor("chr7"), 6703892, 6730431, heatmap = TRUE))
        expect_s3_class(p4, "ggplot")

        expect_silent(p5 <- plot_region(nmr, "chr7", 6703892, 6730431, heatmap = TRUE, heatmap_subsample = 5))
        expect_s3_class(p5, "ggplot")

        for (i in 1:nrow(params)) {
            expect_silent(
                p <- plot_region(
                    x, "chr7", 6703892, 6730431,
                    heatmap = params$heatmap[i],
                    spaghetti = params$spaghetti[i]
                )
            )

            expect_s3_class(p, "ggplot")
        }
    }

    expect_warning(expect_error(plot_region(nmr, "unknown_chr", 6703892, 6730431)))
    expect_error(plot_region(mbr, "unknown_chr", 6703892, 6730431))
})
