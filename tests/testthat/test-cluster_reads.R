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
    expect_no_warning(plot_clustered_reads(x, chr, start, end))
    expect_error(plot_clustered_reads(x, chr, start, end, min_pts = 1))
})

test_that("cluster_reads errors on empty region", {
    x <- load_example_modbamresult()

    expect_error(
        cluster_reads(x, "chr7", 1, 1000),
        "No reads containing methylation data"
    )
})

test_that("cluster_reads errors when too few reads before filtering", {
    x <- load_example_modbamresult()

    expect_error(
        cluster_reads(x, "chr7", 6713892, 6720421, min_pts = 20),
        "Insufficient reads for clustering"
    )
})

test_that("cluster_reads errors when too few reads after filtering", {
    x <- load_example_modbamresult()

    expect_error(
        cluster_reads(x, "chr7", 6713892, 6720421, min_pts = 7),
        "Insufficient reads after filtering"
    )
})

test_that("cluster_reads keeps matrix shape when a single read survives filtering", {
    x <- load_example_modbamresult()

    # only one of the two reads in this region passes the missingness filter
    expect_error(
        cluster_reads(x, "chr7", 6718511, 6718843, min_pts = 2),
        "Insufficient reads after filtering: 1 reads remaining"
    )
})
