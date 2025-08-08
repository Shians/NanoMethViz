test_that("Plotting region methylation heatmap works", {
    # setup
    nmr <- load_example_nanomethresult()

    # test
    expect_no_warning(p <- plot_region_heatmap(nmr, "chr7", 6703892, 6730431))
    expect_s3_class(p, "ggplot")

    expect_no_warning(p <- plot_region_heatmap(nmr, "chr7", 6703892, 6730431, pos_style = "compact"))
    expect_s3_class(p, "ggplot")
})

test_that("plot_region_heatmap works with different parameters", {
    nmr <- load_example_nanomethresult()

    # test with different window_prop values
    expect_no_warning(p <- plot_region_heatmap(nmr, "chr7", 6703892, 6730431, window_prop = 0.1))
    expect_s3_class(p, "ggplot")

    # test with vector window_prop
    expect_no_warning(p <- plot_region_heatmap(nmr, "chr7", 6703892, 6730431, window_prop = c(0.05, 0.15)))
    expect_s3_class(p, "ggplot")

    # test with different subsample value
    expect_no_warning(p <- plot_region_heatmap(nmr, "chr7", 6703892, 6730431, subsample = 30))
    expect_s3_class(p, "ggplot")

    # test with factor chromosome
    chr_factor <- as.factor("chr7")
    expect_no_warning(p <- plot_region_heatmap(nmr, chr_factor, 6703892, 6730431, pos_style = "to_scale"))
    expect_s3_class(p, "ggplot")
})

test_that("Plotting works with GRanges", {
    nmr <- load_example_nanomethresult()
    gr <- GenomicRanges::GRanges(
        seqnames = "chr7",
        ranges = IRanges::IRanges(start = 6703892, end = 6730431),
        strand = "*"
    )

    expect_no_warning(p <- plot_grange_heatmap(nmr, gr))
    expect_s3_class(p, "ggplot")

    expect_no_warning(p <- plot_grange_heatmap(nmr, gr, pos_style = "compact"))
    expect_s3_class(p, "ggplot")
})

test_that("Plotting region methylation works without exons", {
    # setup
    nmr <- load_example_nanomethresult()
    nmr@exons <- tibble::tibble(
        gene_id = character(0),
        chr = character(0),
        strand = character(0),
        start = integer(0),
        end = integer(0),
        transcript_id = character(0),
        symbol = character(0))

    expect_no_warning(p <- plot_region_heatmap(nmr, "chr7", 6703892, 6730431))
})

test_that("plot_region_heatmap window calculation works correctly", {
    nmr <- load_example_nanomethresult()

    # test that window calculations work correctly by testing edge regions
    # small region to test window calculation
    expect_no_warning(p <- plot_region_heatmap(nmr, "chr7", 6720000, 6720100, window_prop = 0.5))
    expect_s3_class(p, "ggplot")

    # test zero window
    expect_no_warning(p <- plot_region_heatmap(nmr, "chr7", 6720000, 6720100, window_prop = 0))
    expect_s3_class(p, "ggplot")
})
