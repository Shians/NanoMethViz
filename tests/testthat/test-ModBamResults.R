test_that("ModBamResults getters and setters work", {
    # setup
    mbr <- load_example_modbamresult()

    # test
    expect_s4_class(NanoMethViz::methy(mbr), "ModBamFiles")
    expect_s3_class(NanoMethViz::exons(mbr), "data.frame")
    expect_s3_class(NanoMethViz::samples(mbr), "data.frame")

    expect_no_warning(
        ModBamResult(
            NanoMethViz::methy(mbr),
            NanoMethViz::samples(mbr)
        )
    )

    expect_no_warning(methy(mbr) <- methy(mbr))
    expect_no_warning(samples(mbr) <- samples(mbr))
    expect_no_warning(exons(mbr) <- exons(mbr))
    expect_error(methy(mbr) <- "invalid_path")
    expect_error(
        exons(mbr) <- dplyr::select(exons(mbr), -"strand"),
        regexp = "columns missing from .+: strand"
    )

    expect_error(
        ModBamFiles(
            paths = system.file(package = "NanoMethViz", "missing.bam", mustWork = FALSE),
            samples = "sample1"
        ),
        regexp = "File path .+ does not exist"
    )

    expect_error(
        ModBamFiles(
            paths = system.file(package = "NanoMethViz", "no_index.bam", mustWork = FALSE),
            samples = "sample1"
        ),
        regexp = ".+ is missing its index file"
    )
})

test_that("ModBamFiles validates its inputs", {
    # setup
    bam_path <- system.file(package = "NanoMethViz", "peg3.bam", mustWork = TRUE)

    # test
    expect_error(
        ModBamFiles(samples = c("sample1", "sample2"), paths = bam_path),
        regexp = "Length of samples \\(2\\) must equal length of paths \\(1\\)"
    )
    expect_error(
        ModBamFiles(samples = character(), paths = character()),
        regexp = "At least one sample and path must be provided"
    )
    expect_error(
        ModBamFiles(samples = NA_character_, paths = bam_path),
        regexp = "Sample names cannot be empty or NA"
    )
    expect_error(
        ModBamFiles(samples = "", paths = bam_path),
        regexp = "Sample names cannot be empty or NA"
    )
    expect_error(
        ModBamFiles(samples = c("sample1", "sample1"), paths = c(bam_path, bam_path)),
        regexp = "Found duplicate sample names: sample1"
    )
})

test_that("ModBamFiles show method works", {
    # setup
    mbr <- load_example_modbamresult()

    # test
    expect_output(show(methy(mbr)), "A ModBamFiles object containing 1 samples")
})

test_that("ModBamResult mod_code setter works", {
    # setup
    mbr <- load_example_modbamresult()

    # test
    expect_identical(mod_code(mbr), "m")

    mod_code(mbr) <- "h"
    expect_identical(mod_code(mbr), "h")

    expect_error(
        mod_code(mbr) <- "",
        regexp = "Modification code cannot be empty or NA"
    )
    expect_error(
        mod_code(mbr) <- NA_character_,
        regexp = "Modification code cannot be empty or NA"
    )
    expect_error(
        mod_code(mbr) <- "mh",
        regexp = "Modification code must be a single character. Got: 'mh'"
    )
})

test_that("ModBamResult constructor validates its inputs", {
    # setup
    mbr <- load_example_modbamresult()
    bam_files <- methy(mbr)
    sample_anno <- samples(mbr)

    # test
    expect_error(
        ModBamResult(methy = as.data.frame(bam_files), samples = sample_anno),
        regexp = "The 'methy' argument must be a ModBamFiles object"
    )
    expect_error(
        ModBamResult(methy = bam_files, samples = sample_anno, mod_code = ""),
        regexp = "Modification code cannot be empty or NA"
    )
    expect_error(
        ModBamResult(methy = bam_files, samples = sample_anno, mod_code = NA_character_),
        regexp = "Modification code cannot be empty or NA"
    )
    expect_error(
        ModBamResult(methy = bam_files, samples = sample_anno, mod_code = "mh"),
        regexp = "Modification code must be a single character. Got: 'mh'"
    )
})

test_that("ModBamResult constructor checks sample matching", {
    # setup
    bam_path <- system.file(package = "NanoMethViz", "peg3.bam", mustWork = TRUE)
    one_bam <- ModBamFiles(samples = "sample1", paths = bam_path)
    two_bams <- ModBamFiles(
        samples = c("sample1", "sample2"),
        paths = c(bam_path, bam_path)
    )
    one_anno <- tibble::tibble(sample = "sample1", group = "group1")
    two_anno <- tibble::tibble(
        sample = c("sample1", "sample2"),
        group = c("group1", "group2")
    )

    # test
    expect_error(
        ModBamResult(
            methy = one_bam,
            samples = tibble::tibble(sample = "other", group = "group1")
        ),
        regexp = "No sample names match between ModBamFiles and sample annotation"
    )

    expect_warning(
        suppressMessages(ModBamResult(methy = two_bams, samples = one_anno)),
        regexp = "Found 1 samples in ModBamFiles not in annotation: sample2"
    )

    expect_warning(
        suppressMessages(ModBamResult(methy = one_bam, samples = two_anno)),
        regexp = "Found 1 samples in annotation not in ModBamFiles: sample2"
    )

    expect_message(
        mbr <- ModBamResult(methy = two_bams, samples = two_anno),
        regexp = "Successfully created ModBamResult with 2 matched samples"
    )
    expect_s4_class(mbr, "ModBamResult")
    expect_s3_class(samples(mbr)$group, "factor")
})
