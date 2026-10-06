test_that("reformat_modkit converts 0-based ref_position to 1-based pos", {
    modkit_data <- tibble::tibble(
        read_id = c("read1", "read2", "read3"),
        forward_read_position = c(10, 20, 30),
        ref_position = c(0, 5, -1),
        chrom = c("chr1", "chr1", "chr1"),
        mod_strand = c("+", "+", "+"),
        ref_strand = c("+", "+", "+"),
        ref_mod_strand = c("+", "+", "+"),
        fw_soft_clipped_start = c(0, 0, 0),
        fw_soft_clipped_end = c(100, 100, 100),
        read_length = c(100, 100, 100),
        mod_qual = c(0.9, 0.5, 0.8),
        mod_code = c("m", "m", "m"),
        base_qual = c(30, 30, 30),
        ref_kmer = c("ACG", "CGT", "GCA"),
        query_kmer = c("ACG", "CGT", "GCA"),
        canonical_base = c("C", "C", "C"),
        modified_primary_base = c("m", "m", "m"),
        inferred = c(FALSE, FALSE, FALSE),
        flag = c(0, 0, 0)
    )

    result <- reformat_modkit(modkit_data, sample = "sample1")

    expect_equal(nrow(result), 2)
    expect_equal(result$pos[result$read_name == "read1"], 1L)
    expect_equal(result$pos[result$read_name == "read2"], 6L)
})

test_that("reformat_modkit filters rows by mod_code", {
    make_modkit_data <- function() {
        tibble::tibble(
            read_id = c("read1", "read2", "read3", "read1", "read2", "read3"),
            forward_read_position = c(10, 20, 30, 40, 50, 60),
            ref_position = c(0, 5, 10, 15, 20, -1),
            chrom = "chr1",
            mod_strand = "+",
            ref_strand = "+",
            ref_mod_strand = "+",
            fw_soft_clipped_start = 0,
            fw_soft_clipped_end = 100,
            read_length = 100,
            mod_qual = c(0.9, 0.5, 0.8, 0.7, 0.6, 0.4),
            mod_code = c("m", "m", "h", "a", "m", "h"),
            base_qual = 30,
            ref_kmer = "ACG",
            query_kmer = "ACG",
            canonical_base = "C",
            modified_primary_base = "m",
            inferred = FALSE,
            flag = 0
        )
    }

    result <- reformat_modkit(make_modkit_data(), sample = "sample1")
    expect_equal(nrow(result), 3)
    expect_setequal(result$read_name, c("read1", "read2"))
    expect_setequal(result$pos, c(1L, 6L, 21L))

    result_h <- reformat_modkit(make_modkit_data(), sample = "sample1", mod_code = "h")
    expect_equal(nrow(result_h), 1)
    expect_equal(result_h$read_name, "read3")
    expect_equal(result_h$pos, 11L)

    result_a <- reformat_modkit(make_modkit_data(), sample = "sample1", mod_code = "a")
    expect_equal(nrow(result_a), 1)
    expect_equal(result_a$read_name, "read1")
    expect_equal(result_a$pos, 16L)

    result_no_match <- reformat_modkit(
        make_modkit_data(),
        sample = "sample1",
        mod_code = "nonexistent"
    )
    expect_equal(nrow(result_no_match), 0)
})

test_that("resolve_mod_code defaults to 5mC for modkit and NULL otherwise", {
    expect_equal(resolve_mod_code(NULL, "modkit"), "m")
    expect_equal(resolve_mod_code(NULL, c("nanopolish", "modkit")), "m")
    expect_null(resolve_mod_code(NULL, c("nanopolish", "f5c")))
    expect_equal(resolve_mod_code("h", "modkit"), "h")
    expect_equal(resolve_mod_code("21839", "modkit"), "21839")
})

test_that("resolve_mod_code rejects mod_code without modkit input", {
    expect_error(
        resolve_mod_code("m", c("nanopolish", "megalodon")),
        "detected: nanopolish, megalodon"
    )
})

test_that("resolve_mod_code rejects malformed mod_code", {
    expect_error(resolve_mod_code(c("m", "h"), "modkit"), "single non-empty string")
    expect_error(resolve_mod_code(NA_character_, "modkit"), "single non-empty string")
    expect_error(resolve_mod_code("", "modkit"), "single non-empty string")
})
