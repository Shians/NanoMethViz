test_that("assert_has_index accepts .bai and .csi indexes", {
    bam <- tempfile()
    file.create(bam)
    withr::defer(unlink(c(bam, paste0(bam, ".bai"), paste0(bam, ".csi"))))

    expect_error(assert_has_index(bam), "[.]bai or [.]csi")

    file.create(paste0(bam, ".csi"))
    expect_error(assert_has_index(bam), NA)

    unlink(paste0(bam, ".csi"))
    file.create(paste0(bam, ".bai"))
    expect_error(assert_has_index(bam), NA)
})

test_that("assert_has_index reports missing index for multiple files", {
    bams <- tempfile(pattern = c("bam1", "bam2"))
    file.create(bams)
    withr::defer(unlink(c(bams, paste0(bams, ".bai"), paste0(bams, ".csi"))))

    expect_error(assert_has_index(bams), "are missing their index files")
})
