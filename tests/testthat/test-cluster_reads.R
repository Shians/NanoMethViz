test_that("cluster_reads assertions", {
    x <- load_example_modbamresult()
    chr <- "chr7"
    start <- 6713892
    end <- 6720421
    min_pts <- 5

    # Assertion tests
    expect_error(cluster_reads(x, chr, start, end, "min_pts"))

    expect_error(cluster_reads(x, chr, "start", end, min_pts))

    expect_error(cluster_reads(x, chr, start, end, -2))

    # Successful assertion test
    expect_no_error(cluster_reads(x, chr, start, end, min_pts))
    plot_clustered_reads(x, chr, start, end)
})
