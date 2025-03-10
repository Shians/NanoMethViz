test_that("Modbam to tabix conversion works", {
    out_file <- paste0(tempfile(), ".tsv.bgz")
    mbr <- ModBamResult(
        methy = ModBamFiles(
            samples = "sample1",
            paths = system.file("peg3.bam", package = "NanoMethViz", mustWork = FALSE)
        ),
        samples = data.frame(
            sample = "sample1",
            group = "group1"
        )
    )

    expect_no_error(modbam_to_tabix(mbr, out_file))
    expect_true(file_exists(out_file))

    tabix_data <- expect_no_error(read_tsv(out_file, col_names = methy_col_names()))
    expect_equal(nrow(tabix_data), 10371)
    expect_equal(ncol(tabix_data), 6)
    expect_equal(unique(tabix_data$sample), "sample1")

    fs::file_delete(out_file)
})


test_that("Modbam to tabix error checking works", {
    out_file <- paste0(tempfile())
    mbr <- ModBamResult(
        methy = ModBamFiles(
            samples = "sample1",
            paths = system.file("peg3.bam", package = "NanoMethViz", mustWork = FALSE)
        ),
        samples = data.frame(
            sample = "sample1",
            group = "group1"
        )
    )

    expect_error(modbam_to_tabix(mbr, out_file), "output_file must end with .bgz extension.")
    expect_false(file_exists(out_file))

    out_folder <- paste0(tempfile(), ".tsv.bgz")
    fs::dir_create(out_folder)

    expect_error(modbam_to_tabix(mbr, out_folder), "output_file exists and is not a file")
})
