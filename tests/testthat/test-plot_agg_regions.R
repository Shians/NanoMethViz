test_that("Aggregate plotting works", {
    # setup
    nmr <- load_example_nanomethresult()
    gene_anno <- exons_to_genes(NanoMethViz::exons(nmr))

    # test
    expect_silent(plot_agg_regions(nmr, gene_anno))
    expect_silent(plot_agg_regions(nmr, gene_anno, "sample"))
    expect_silent(plot_agg_regions(nmr, gene_anno, "group"))
})

test_that("Aggregate plotting error checking works", {
    # setup
    nmr <- load_example_nanomethresult()
    missing_chr <- NanoMethViz::exons(nmr)[, -2]
    missing_chr_list <- list(
        NanoMethViz::exons(nmr)[, -2],
        NanoMethViz::exons(nmr),
        NanoMethViz::exons(nmr)[, -2]
    )

    # test
    expect_error(plot_agg_regions(nmr, gene_anno, "foo"))
})
