test_that("Aggregate plotting works", {
    # setup
    nmr <- load_example_nanomethresult()
    gene_anno <- exons_to_genes(NanoMethViz::exons(nmr))

    # test
    expect_silent(plot_agg_regions(nmr, gene_anno))
    expect_silent(plot_agg_regions(nmr, list(gene_anno, gene_anno)))
    expect_silent(plot_agg_regions(nmr, list(set1=gene_anno, set2=gene_anno)))
    expect_silent(plot_agg_regions_sample_grouped(nmr, gene_anno))
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
    expect_error(plot_agg_regions(nmr, missing_chr))
    expect_error(plot_agg_regions(nmr, missing_chr_list))
    expect_error(plot_agg_regions_sample_grouped(nmr, missing_chr))
    expect_error(plot_agg_regions_sample_grouped(nmr, missing_chr_list))
})
