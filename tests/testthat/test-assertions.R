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

test_that("assert_readable accepts existing files", {
    paths <- tempfile(pattern = c("file1", "file2"))
    file.create(paths)
    withr::defer(unlink(paths))

    expect_no_error(assert_readable(paths))
})

test_that("assert_readable reports a single missing file", {
    existing <- tempfile()
    file.create(existing)
    withr::defer(unlink(existing))
    missing <- tempfile()

    expect_error(
        assert_readable(c(existing, missing)),
        paste0("File path '", missing, "' does not exist"),
        fixed = TRUE
    )
})

test_that("assert_readable reports multiple missing files", {
    missing <- tempfile(pattern = c("missing1", "missing2"))

    expect_error(
        assert_readable(missing),
        "File paths '.+', '.+' do not exist"
    )
})

test_that("assert_valid_methy_samples errors on unreadable file", {
    samples <- data.frame(sample = "A", group = "1")
    missing <- tempfile(fileext = ".tsv.bgz")

    expect_error(
        suppressWarnings(assert_valid_methy_samples(missing, samples)),
        "Failed to read methylation data file"
    )
})

test_that("assert_valid_methy_samples errors when no samples match", {
    methy <- system.file("methy_subset.tsv.bgz", package = "NanoMethViz")
    samples <- data.frame(sample = c("X", "Y"), group = c("1", "2"))

    expect_error(
        assert_valid_methy_samples(methy, samples),
        "No sample names from the data file match the sample annotation"
    )
})

methy_subset_samples <- c(
    "B6Cast_Prom_1_bl6", "B6Cast_Prom_1_cast",
    "B6Cast_Prom_2_bl6", "B6Cast_Prom_2_cast",
    "B6Cast_Prom_3_bl6", "B6Cast_Prom_3_cast"
)

test_that("assert_valid_methy_samples succeeds when all samples match", {
    methy <- system.file("methy_subset.tsv.bgz", package = "NanoMethViz")
    samples <- data.frame(sample = methy_subset_samples, group = "1")

    expect_message(
        assert_valid_methy_samples(methy, samples),
        "Successfully matched 6 samples"
    )
})

test_that("assert_valid_methy_samples warns about data samples missing from annotation", {
    methy <- system.file("methy_subset.tsv.bgz", package = "NanoMethViz")
    samples <- data.frame(sample = methy_subset_samples[1:2], group = c("1", "2"))

    expect_message(
        expect_warning(
            assert_valid_methy_samples(methy, samples),
            "Found 4 samples in data file not present in annotation"
        ),
        "Successfully matched 2 samples"
    )
})

test_that("assert_valid_methy_samples warns about annotation samples missing from data", {
    methy <- system.file("methy_subset.tsv.bgz", package = "NanoMethViz")
    samples <- data.frame(sample = c(methy_subset_samples, "extra"), group = "1")

    expect_message(
        expect_warning(
            assert_valid_methy_samples(methy, samples),
            "Found 1 samples in annotation not present in data file: extra"
        ),
        "Successfully matched 6 samples"
    )
})
