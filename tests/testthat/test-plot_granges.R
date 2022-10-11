test_that("Plotting GRanges works", {
    nmr <- load_example_nanomethresult()

    expect_silent(plot_grange(nmr, GenomicRanges::GRanges("chr7:6703892-6730431")))
    expect_silent(plot_grange(nmr, GenomicRanges::GRanges("chr7:6703892-6730431", spaghetti = TRUE)))
    expect_silent(plot_grange(nmr, GenomicRanges::GRanges("chr7:6703892-6730431", heatmap = TRUE)))
    expect_silent(plot_grange(nmr, GenomicRanges::GRanges("chr7:6703892-6730431", spagheti = TRUE, heatmap = TRUE)))
    expect_silent(plot_grange_heatmap(nmr, GenomicRanges::GRanges("chr7:6703892-6730431")))
})
