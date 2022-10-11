test_that("Filtering methylation files works", {
    nmr <- load_example_nanomethresult()
    output_file <- paste0(tempfile(), ".tsv.bgz")

    expect_no_error(filter_methy(nmr, output_file = output_file, chr == "chrX"))
    expect_no_error(filter_methy(nmr, output_file = output_file, chr == "chrX", pos > 1))
    expect_no_error(filter_methy(methy(nmr), output_file = output_file, chr == "chrX"))
})
