test_that("Plotting region methylation heatmap works", {
    # setup
    nmr <- load_example_nanomethresult()

    # test
    expect_silent(p <- plot_region_heatmap(nmr, "chr7", 6703892, 6730431))
    expect_true(is(p, "ggplot"))
})

test_that("Plotting region methylation works without exons", {
    # setup
    nmr <- load_example_nanomethresult()
    nmr@exons <- tibble(
        gene_id = character(0),
        chr = character(0),
        strand = character(0),
        start = integer(0),
        end = integer(0),
        transcript_id = character(0),
        symbol = character(0))

    expect_silent(p <- plot_region_heatmap(nmr, "chr7", 6703892, 6730431))
})

