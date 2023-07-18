test_that("Plotting gene works", {
    # setup
    nmr <- load_example_nanomethresult()
    mbr <- load_example_modbamresult()
    params <- expand.grid(
        heatmap = c(TRUE, FALSE),
        spaghetti = c(TRUE, FALSE),
        gene_anno = c(TRUE, FALSE)
    )

    # test plot_gene() ----
    for (x in list(nmr, mbr)) {
        expect_silent(p <- plot_gene(x, "Peg3"))
        expect_s3_class(p, "ggplot")

        expect_silent(p <- plot_gene(x, "Peg3", heatmap = TRUE))
        expect_s3_class(p, "ggplot")

        expect_silent(p <- plot_gene(x, "Peg3", heatmap = TRUE, heatmap_subsample = 5))
        expect_s3_class(p, "ggplot")

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

    # test plot_region() ----
    for (x in list(nmr, mbr)) {
        expect_silent(p <- plot_region(x, "chr7", 6703892, 6730431))
        expect_s3_class(p, "ggplot")

        expect_silent(p <- plot_region(x, factor("chr7"), 6703892, 6730431))
        expect_s3_class(p, "ggplot")

        expect_silent(p <- plot_region(nmr, "chr7", 6703892, 6730431, heatmap = TRUE, heatmap_subsample = 5))
        expect_s3_class(p, "ggplot")

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
        expect_error(plot_region(x, chr = "chr7", start = 0, end = 0), "Elements 1 of start < end are not true")
        expect_warning(plot_region(x, chr = "chr7", start = 0, end = 10), "no methylation data in region")
    }

    # test plot_region() with unknown chromosome ----
    expect_warning(expect_error(plot_region(nmr, "unknown_chr", 6703892, 6730431)))
    expect_error(plot_region(mbr, "unknown_chr", 6703892, 6730431))

    # test plot_grange() ----
    grange <- GenomicRanges::GRanges("chr7:6703892-6730431")
    for (x in list(nmr, mbr)) {
        for (i in 1:nrow(params)) {
            expect_silent(
                p <- plot_grange(
                    x, grange,
                    heatmap = params$heatmap[i],
                    spaghetti = params$spaghetti[i]
                )
            )

            expect_s3_class(p, "ggplot")
        }
    }
})
