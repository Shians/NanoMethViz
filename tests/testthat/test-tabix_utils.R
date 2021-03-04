test_that("multiplication works", {
    # setup
    methy_calls <- system.file(package = "NanoMethViz",
        c("sample1_nanopolish.tsv.gz", "sample2_nanopolish.tsv.gz"))
    temp_file <- paste0(tempfile(), ".tsv.bgz")
    withr::defer(file.remove(temp_file))

    # test
    expect_message(create_tabix_file(methy_calls, temp_file))
    expect_true(is(methy_to_bsseq(temp_file), "BSseq"))
})
