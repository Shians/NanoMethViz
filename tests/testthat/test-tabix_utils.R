test_that("Importing nanopolish works", {
    # setup
    methy_calls <- system.file(package = "NanoMethViz",
        c("sample1_nanopolish.tsv.gz", "sample2_nanopolish.tsv.gz"), mustWork = FALSE)
    temp_file <- paste0(tempfile(), ".tsv.bgz")
    withr::defer(file.remove(temp_file))

    # test
    expect_message(create_tabix_file(methy_calls, temp_file))
    expect_s4_class(methy_to_bsseq(temp_file), "BSseq")
})

test_that("Importing megalodon works", {
    # setup
    methy_calls <- system.file(package = "NanoMethViz",
       "megalodon_calls.txt.gz", mustWork = FALSE)
    temp_file <- paste0(tempfile(), ".tsv.bgz")
    withr::defer(file.remove(temp_file))

    # test
    expect_message(create_tabix_file(methy_calls, temp_file))
})

test_that("sort_methy_file stops when the sort command fails", {
    skip_on_os("windows")

    x <- tempfile()
    writeLines(c("sample1\tchr1\t1\t+\t0.5\tread1"), x)
    withr::defer(unlink(x))

    Sys.chmod(x, "0000")
    expect_error(sort_methy_file(x), "Failed to sort methylation file")
})
