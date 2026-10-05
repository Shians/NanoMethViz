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
