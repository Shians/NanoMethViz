test_that("Modbam to tabix conversion works", {
    out_file <- paste0(tempfile(), ".tsv.bgz")
    mbr <- ModBamResult(
        methy = ModBamFiles(
            samples = "sample1",
            paths = system.file("peg3.bam", package = "NanoMethViz")
        ),
        samples = data.frame(
            sample = "sample1",
            group = "group1"
        )
    )
    
    expect_no_error(modbam_to_tabix(mbr, out_file))
    expect_true(file_exists(out_file))
    
    expect_no_error(tabix_data <- read_tsv(out_file, col_names = methy_col_names()))
    expect_equal(nrow(tabix_data), 10371)
    expect_equal(ncol(tabix_data), 6)
    expect_equal(unique(tabix_data$sample), "sample1")
    
    fs::file_delete(out_file)
})
