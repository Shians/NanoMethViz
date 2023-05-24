test_that("Plotting region works", {
    # setup
    nmr <- load_example_nanomethresult()
    mbr <- ModBamResult(
        methy = ModBamFiles(
            paths = system.file(package = "NanoMethViz", "peg3.bam"),
            samples = "sample1"
        ),
        samples = tibble::tibble(
            sample = "sample1",
            group = "group1"
        )
    )

    # test
    expect_silent(p1 <- plot_region(nmr, "chr7", 6703892, 6730431))
    expect_true(is(p1, "ggplot"))
    expect_silent(p2 <- plot_region(mbr, "chr7", 6703892, 6730431))
    expect_true(is(p2, "ggplot"))

    expect_silent(plot_region(nmr, factor("chr7"), 6703892, 6730431))
    expect_silent(plot_region(mbr, factor("chr7"), 6703892, 6730431))

    params <- expand.grid(
        heatmap = c(TRUE, FALSE),
        spaghetti = c(TRUE, FALSE),
        gene_anno = c(TRUE, FALSE)
    )

    for (i in 1:nrow(params)) {
        expect_silent(
            plot_region(
                nmr, "chr7", 6703892, 6730431,
                heatmap = params$heatmap[i],
                spaghetti = params$spaghetti[i]
            )
        )
    }

    expect_warning(expect_error(plot_region(nmr, "unknown_chr", 6703892, 6730431)))
    expect_error(plot_region(mbr, "unknown_chr", 6703892, 6730431))
})
