test_that("cluster_reads assertions", {
    x <- ModBamResult(
        methy = ModBamFiles(
            paths = system.file(package = "NanoMethViz", "peg3.bam"),
            samples = "sample1"
        ),
        samples = tibble::tibble(
            sample = "sample1",
            group = "group1"
        )
    )
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
})
