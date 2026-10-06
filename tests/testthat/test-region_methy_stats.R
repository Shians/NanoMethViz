test_that("region_methy_stats calculates methylation statistics for each region", {
    # Create sample data
    nmr <- load_example_nanomethresult()
    regions <- exons_to_genes(exons(nmr))

    # Call the function
    result <- region_methy_stats(nmr, regions, threshold = 0.6)

    # Assert the result has the expected structure and columns
    expect_s3_class(result, "data.frame")
    expect_identical(colnames(result), c("gene_id", "chr", "strand", "symbol", "start", "end", "mean_methy_prob", "prevalence"))

    # Assert the number of rows in the result matches the number of input regions
    expect_equal(nrow(result), nrow(regions))

    # Assert that the calculated statistics are within the expected range
    expect_true(all(result$mean_methy_prob >= 0 & result$mean_methy_prob <= 1))
    expect_true(all(result$prevalence >= 0 & result$prevalence <= 1))
})

test_that("region_methy_stats returns NA statistics for regions without data", {
    # setup
    nmr <- load_example_nanomethresult()
    regions <- tibble::tibble(
        chr = c("chr7", "chr7"),
        start = c(6703892, 1),
        end = c(6730431, 10)
    )

    # test
    result <- region_methy_stats(nmr, regions)
    expect_equal(nrow(result), 2)
    expect_false(is.na(result$mean_methy_prob[1]))
    expect_false(is.na(result$prevalence[1]))
    expect_true(is.na(result$mean_methy_prob[2]))
    expect_true(is.na(result$prevalence[2]))
})
