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
    expect_true(file_exists(paste0(out_file, ".tbi")))

    tabix_data <- expect_no_error(read_tsv(out_file, col_names = methy_col_names()))
    expect_equal(nrow(tabix_data), 10371)
    expect_equal(ncol(tabix_data), 6)
    expect_equal(unique(tabix_data$sample), "sample1")

    expect_no_error(query_methy(out_file, "chr7", 6713552, 6730431))

    fs::file_delete(out_file)
    fs::file_delete(paste0(out_file, ".tbi"))
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

test_that("Modbam to tabix creates missing output directories", {
    out_file <- file.path(withr::local_tempdir(), "a", "b", "out.tsv.bgz")
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

    expect_false(fs::dir_exists(fs::path_dir(out_file)))
    expect_no_error(modbam_to_tabix(mbr, out_file))
    expect_true(fs::file_exists(out_file))
    expect_true(fs::file_exists(paste0(out_file, ".tbi")))
})

test_that("Modbam to tabix overwrites existing output file", {
    out_file <- file.path(withr::local_tempdir(), "out.tsv.bgz")
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

    writeLines("stale content", out_file)

    expect_message(modbam_to_tabix(mbr, out_file), "will overwrite")
    expect_true(fs::file_exists(paste0(out_file, ".tbi")))

    tabix_data <- read_tsv(out_file, col_names = methy_col_names(), show_col_types = FALSE)
    expect_equal(nrow(tabix_data), 10371)
    expect_equal(unique(tabix_data$sample), "sample1")
})
