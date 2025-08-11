test_that("genomic coordinate validation works", {
    # Test valid coordinates
    expect_silent(assert_valid_genomic_coords("chr1", 1, 100))
    expect_silent(assert_valid_genomic_coords(c("chr1", "chr2"), c(1, 50), c(100, 150)))

    # Test invalid: start > end
    expect_error(
        assert_valid_genomic_coords("chr1", 100, 50),
        "Start positions must be < end positions"
    )

    # Test invalid: negative coordinates
    expect_error(
        assert_valid_genomic_coords("chr1", -5, 100),
        "Start positions must be positive"
    )

    expect_error(
        assert_valid_genomic_coords("chr1", 1, -5),
        "End positions must be positive"
    )

    # Test invalid: NA chromosome
    expect_error(
        assert_valid_genomic_coords(NA, 1, 100),
        "Chromosome cannot be NA"
    )

    # Test invalid: NA coordinates
    expect_error(
        assert_valid_genomic_coords("chr1", NA, 100),
        "Start and end positions cannot be NA"
    )

    # Test string to numeric conversion
    expect_silent(assert_valid_genomic_coords("chr1", "1", "100"))

    # Test invalid string conversion
    expect_error(
        assert_valid_genomic_coords("chr1", "abc", 100),
        "Start position must be numeric"
    )

    # Test allow_equal = TRUE
    expect_silent(assert_valid_genomic_coords("chr1", 100, 100, allow_equal = TRUE))
    expect_error(
        assert_valid_genomic_coords("chr1", 100, 100, allow_equal = FALSE),
        "Start positions must be < end positions"
    )
})

test_that("sample validation works", {
    # Test valid samples
    valid_samples <- data.frame(sample = c("A", "B", "C"), group = c("1", "2", "1"))
    expect_silent(assert_valid_samples(valid_samples))

    # Test duplicate sample names
    dup_samples <- data.frame(sample = c("A", "B", "A"), group = c("1", "2", "1"))
    expect_error(
        assert_valid_samples(dup_samples),
        "Found duplicate sample names"
    )

    # Test empty sample names
    empty_samples <- data.frame(sample = c("A", "", "C"), group = c("1", "2", "1"))
    expect_error(
        assert_valid_samples(empty_samples),
        "Sample names cannot be empty or NA"
    )

    # Test missing required columns
    missing_cols <- data.frame(sample = c("A", "B", "C"))
    expect_error(
        assert_valid_samples(missing_cols),
        "columns missing from"
    )

    # Test non-data.frame input
    expect_error(
        assert_valid_samples("not_a_df"),
        "must be a data.frame"
    )

    # Test empty data.frame
    empty_df <- data.frame(sample = character(0), group = character(0))
    expect_error(
        assert_valid_samples(empty_df),
        "data.frame is empty"
    )
})

test_that("exon validation works", {
    # Test valid exons
    valid_exons <- data.frame(
        gene_id = "GENE1",
        chr = "chr1",
        strand = "+",
        start = 1000,
        end = 2000,
        transcript_id = "T1",
        symbol = "GENE1"
    )
    expect_silent(assert_valid_exons(valid_exons))

    # Test invalid strand
    invalid_strand <- valid_exons
    invalid_strand$strand <- "X"
    expect_error(
        assert_valid_exons(invalid_strand),
        "Invalid strand values"
    )

    # Test invalid coordinates in exons
    invalid_coords <- valid_exons
    invalid_coords$start <- 2000
    invalid_coords$end <- 1000
    expect_error(
        assert_valid_exons(invalid_coords),
        "Invalid genomic coordinates found in exon annotation"
    )

    # Test missing required columns
    missing_cols <- valid_exons[, -which(names(valid_exons) == "symbol")]
    expect_error(
        assert_valid_exons(missing_cols),
        "columns missing from"
    )

    # Test non-data.frame input
    expect_error(
        assert_valid_exons("not_a_df"),
        "must be a data.frame"
    )

    # Test empty exons (should not error)
    empty_exons <- data.frame(
        gene_id = character(0),
        chr = character(0),
        strand = character(0),
        start = integer(0),
        end = integer(0),
        transcript_id = character(0),
        symbol = character(0)
    )
    expect_silent(assert_valid_exons(empty_exons))
})

test_that("error messages provide helpful guidance", {
    # Test that error messages contain helpful suggestions
    expect_error(
        assert_valid_genomic_coords("chr1", 0, 100),
        "Genomic coordinates are 1-based in this package"
    )

    expect_error(
        assert_valid_samples(data.frame(sample = c("A", "A"), group = c("1", "2"))),
        "Each sample must have a unique identifier"
    )
})
