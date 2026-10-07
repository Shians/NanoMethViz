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

test_that("reformat_modkit takes strand from the reference and shifts reverse CG calls", {
    # simplex reads have mod_strand "+" whatever their alignment, so strand
    # must come from ref_mod_strand; query_kmer is read-oriented
    modkit_data <- tibble::tibble(
        read_id = c("fwd_cg", "rev_cg", "rev_ca", "rev_6ma"),
        forward_read_position = c(10, 20, 30, 40),
        ref_position = c(99, 200, 300, 400),
        chrom = "chr1",
        mod_strand = "+",
        ref_strand = c("+", "-", "-", "-"),
        ref_mod_strand = c("+", "-", "-", "-"),
        fw_soft_clipped_start = 0,
        fw_soft_clipped_end = 0,
        read_length = 100,
        mod_qual = 0.9,
        mod_code = c("m", "m", "m", "a"),
        base_qual = 30,
        ref_kmer = ".",
        query_kmer = c("GACGG", "TACGT", "TACAT", "GTAGC"),
        canonical_base = c("C", "C", "C", "A"),
        modified_primary_base = c("C", "C", "C", "A"),
        inferred = FALSE,
        flag = c(0, 16, 16, 16)
    )

    result_m <- reformat_modkit(modkit_data, sample = "sample1")
    expect_equal(result_m$read_name, c("fwd_cg", "rev_cg", "rev_ca"))
    expect_equal(as.character(result_m$strand), c("+", "-", "-"))
    # forward CG keeps 1-based pos, reverse CG moves from G (201) to C (200),
    # reverse non-CG stays on its aligned base
    expect_equal(result_m$pos, c(100L, 200L, 301L))

    # non-C modifications are never shifted
    result_a <- reformat_modkit(modkit_data, sample = "sample1", mod_code = "a")
    expect_equal(result_a$pos, 401L)
    expect_equal(as.character(result_a$strand), "-")
})

test_that("modkit import matches modbam_to_tabix on the same BAM", {
    skip_if(Sys.which("modkit") == "", "modkit not installed")

    bam <- system.file("peg3.bam", package = "NanoMethViz")
    extract_file <- withr::local_tempfile(fileext = ".tsv")
    status <- system2(
        "modkit",
        c("extract", "full", "--mapped-only", bam, extract_file),
        stdout = FALSE, stderr = FALSE
    )
    skip_if(status != 0, "modkit extract failed")

    modkit_out <- withr::local_tempfile(fileext = ".tsv")
    convert_methy_format(extract_file, modkit_out, samples = "s1", verbose = FALSE)

    mbr <- ModBamResult(
        ModBamFiles("s1", bam),
        data.frame(sample = "s1", group = "g")
    )
    modbam_out <- withr::local_tempfile(fileext = ".tsv.bgz")
    suppressMessages(modbam_to_tabix(mbr, modbam_out))

    read_calls <- function(path) {
        readr::read_tsv(
            path,
            col_names = methy_col_names(),
            col_types = methy_col_types()
        ) %>%
            dplyr::mutate(strand = as.character(.data$strand)) %>%
            dplyr::select("read_name", "pos", "strand") %>%
            dplyr::arrange(.data$read_name, .data$pos)
    }

    expect_equal(read_calls(modkit_out), read_calls(modbam_out))
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

# Write tab-separated lines to a temporary file cleaned up with the caller
write_methy_lines <- function(header, rows = character(), env = parent.frame()) {
    path <- withr::local_tempfile(fileext = ".tsv", .local_envir = env)
    lines <- c(paste(header, collapse = "\t"), rows)
    writeLines(lines, path)
    path
}

f5c_header <- c(
    "chromosome", "strand", "start", "end", "read_name", "log_lik_ratio",
    "log_lik_methylated", "log_lik_unmethylated", "num_calling_strands",
    "num_cpgs", "sequence"
)

modkit_header <- c(
    "read_id", "forward_read_position", "ref_position", "chrom",
    "mod_strand", "ref_strand", "ref_mod_strand", "fw_soft_clipped_start",
    "fw_soft_clipped_end", "read_length", "mod_qual", "mod_code", "base_qual",
    "ref_kmer", "query_kmer", "canonical_base", "modified_primary_base",
    "inferred", "flag"
)

test_that("guess_methy_source detects f5c", {
    f5c_file <- write_methy_lines(
        f5c_header,
        "chr1\t+\t100\t100\tread1\t2.5\t-10.0\t-12.5\t1\t1\tAAAAACGTTTT"
    )

    expect_equal(guess_methy_source(f5c_file), "f5c")
})

test_that("guess_methy_source detects nanopolish", {
    nanopolish_file <- write_methy_lines(
        c(
            "chromosome", "strand", "start", "end", "read_name",
            "log_lik_ratio", "log_lik_methylated", "log_lik_unmethylated",
            "num_calling_strands", "num_motifs", "sequence"
        ),
        "chr1\t-\t100\t100\tread1\t-5.91\t-100.38\t-94.47\t1\t1\tCATTACGTTTC"
    )

    expect_equal(guess_methy_source(nanopolish_file), "nanopolish")
})

test_that("guess_methy_source detects megalodon", {
    megalodon_file <- write_methy_lines(
        c("read_id", "chrm", "strand", "pos", "mod_log_prob", "can_log_prob", "mod_base"),
        "read1\tIII\t+\t7641294\t-0.665\t-0.722\tY"
    )

    expect_equal(guess_methy_source(megalodon_file), "megalodon")
})

test_that("guess_methy_source detects modkit", {
    modkit_file <- write_methy_lines(
        modkit_header,
        "read1\t10\t0\tchr1\t+\t+\t+\t0\t100\t100\t0.9\tm\t30\tACG\tACG\tC\tm\tfalse\t0"
    )

    expect_equal(guess_methy_source(modkit_file), "modkit")
})

test_that("guess_methy_source errors on unrecognised header", {
    unknown_file <- write_methy_lines(c("foo", "bar", "baz"), "1\t2\t3")

    expect_error(guess_methy_source(unknown_file), "Format not recognised")
})

test_that("convert_methy_format expands multi-CpG f5c calls", {
    # f5c starts are 0-based; the sequence carries 5 bases of flanking context,
    # so the CpGs in the second row sit at 0-based 200 and 204
    f5c_file <- write_methy_lines(
        f5c_header,
        c(
            "chr1\t+\t100\t100\tread1\t2.5\t-10.0\t-12.5\t1\t1\tAAAAACGTTTT",
            "chr1\t+\t200\t204\tread2\t-1.5\t-20.0\t-18.5\t1\t2\tAAAAACGTTCGTTTTT"
        )
    )
    output_file <- withr::local_tempfile(fileext = ".tsv")

    convert_methy_format(f5c_file, output_file, samples = "s1", verbose = FALSE)
    result <- readr::read_tsv(
        output_file,
        col_names = methy_col_names(),
        col_types = methy_col_types()
    )

    expect_equal(nrow(result), 3)
    expect_equal(as.character(result$sample), rep("s1", 3))
    expect_equal(as.character(result$chr), rep("chr1", 3))
    expect_equal(result$pos, c(101L, 201L, 205L))
    expect_equal(as.character(result$strand), rep("*", 3))
    expect_equal(result$statistic, c(2.5, -1.5, -1.5))
    expect_equal(result$read_name, c("read1", "read2", "read2"))
})

test_that("convert_methy_format reports no mod codes for header-only modkit file", {
    modkit_file <- write_methy_lines(modkit_header)
    output_file <- withr::local_tempfile(fileext = ".tsv")

    expect_error(
        convert_methy_format(modkit_file, output_file, samples = "s1", verbose = FALSE),
        "Mod codes found in file: none"
    )
})
