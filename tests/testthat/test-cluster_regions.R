test_that("Plotting feature clustering", {
    nmr <- load_example_nanomethresult()
    gene_anno <- exons_to_genes(NanoMethViz::exons(nmr))
    expect_silent(cluster_regions(nmr, gene_anno, centers = 2, grid_method = "uniform"))
})
