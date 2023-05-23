test_that("make_granges returns correct GRanges object", {
    chr <- c("chr1", "chr2", "chr3")
    start <- c(100, 200, 300)
    end <- c(200, 300, 400)

    # Call the make_granges function
    granges <- make_granges(chr, start, end)

    # Test assertions
    expect_s4_class(granges, "GRanges")
    expect_equal(length(granges), length(chr))
    expect_equal(names(granges), NULL)  # Expect names to be NULL

    # Check each range individually
    for (i in seq_along(granges)) {
        expect_equal(as.character(GenomicRanges::seqnames(granges))[i], chr[i])
        expect_equal(GenomicRanges::start(granges)[i], start[i])
        expect_equal(GenomicRanges::end(granges)[i], end[i])
    }
})

test_that("make_granges throws error for unequal input lengths", {
    chr <- c("chr1", "chr2", "chr3")
    start <- c(100, 200)
    end <- c(200, 300, 400)

    # Expect an error to be thrown for unequal input lengths
    expect_error(make_granges(chr, start, end))
})
