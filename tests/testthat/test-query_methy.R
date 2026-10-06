test_that("Querying methylation works", {
    # setup
    nmr <- load_example_nanomethresult()
    mbr <- load_example_modbamresult()

    # test
    # each query is half of Peg3
    queries <- tibble(
        chr = c("chr7", "chr7"),
        start = c(6703892, 6717162),
        end = c(6717161, 6730431)
    )

    # test basic operation
    methy_data_reg <- expect_no_warning(query_methy(nmr, "chr7", 6703892, 6730431))
    expect_equal(
        colnames(methy_data_reg),
        c("sample", "chr", "pos", "strand", "statistic", "read_name", "mod_prob")
    )
    expect_no_warning(query_methy(nmr, queries$chr, queries$start, queries$end))
    expect_no_warning(query_methy(nmr, queries$chr, queries$start, queries$end, simplify = FALSE))

    # test working on direct methy
    query_methy(methy(nmr), queries$chr[1], queries$start[1], queries$end[2])

    # test working on methods
    methy_data_fct <- expect_no_warning(query_methy(nmr, factor(queries$chr[1]), queries$start[1], queries$end[2]))
    expect_equal(methy_data_reg, methy_data_fct)

    # test working on gene
    expect_no_warning(methy_data_gene <- query_methy_gene(nmr, "Peg3"))
    expect_equal(methy_data_fct, methy_data_gene)

    # test warnings and errors
    expect_warning(expect_error(query_methy(nmr, "Missing", 1, 1000)))
    expect_error(query_methy_gene(nmr, "Missing"))

    # test working on table queries
    expect_no_warning(query_methy_df(methy(nmr), queries))

    # test working on GRanges
    regions_gr <- GenomicRanges::GRanges(queries)
    methy_data_gr <- expect_no_warning(query_methy_gr(nmr, regions_gr))
    expect_equal(methy_data_reg, methy_data_gr)

    # test working on modbam
    methy_data_modbam <- expect_no_warning(query_methy(mbr, queries$chr[1], queries$start[1], queries$end[2]))
    expect_equal(colnames(methy_data_reg), colnames(methy_data_modbam))

    queries_warn <- tibble(
        chr = c("chr7", "chr7", "chr7", "foo"),
        start = c(6703892, 6717162, 10000, 1),
        end = c(6717161, 6730431, 20000, 2)
    )

    # test ordering of unsimplified output
    expect_warning(
        methy_data_list <- query_methy(nmr, queries_warn$chr, queries_warn$start, queries_warn$end, simplify = FALSE),
        "requested sequences missing from tabix file and will be excluded from query:foo"
    )
    expect_true(length(methy_data_list) == 4)

    # test to make sure that reading from text connection works for methylation data
    methy_data_with_inf <- readLines(system.file("methy_data_with_inf.tsv", package = "NanoMethViz", mustWork = FALSE))
    methy_data <- expect_no_warning(read_methy_lines(methy_data_with_inf))
    expect_equal(nrow(methy_data), length(methy_data_with_inf))

    # test when query regions are empty
    queries_with_empty <- tibble(
        chr = c("chr7", "chr7", "chr7"),
        start = c(6703892, 6717162, 1),
        end = c(6717161, 6730431, 10)
    )
    expect_s3_class(query_methy_df(mbr, queries_with_empty), "data.frame")
})

get_query_test_modbamresult <- function() {
    ModBamResult(
        methy = ModBamFiles(
            paths = system.file("peg3.bam", package = "NanoMethViz", mustWork = FALSE),
            samples = "sample1"
        ),
        samples = tibble::tibble(sample = "sample1", group = "group1")
    )
}

test_that("duplicated regions return one output per region", {
    # setup
    methy_path <- system.file("methy_subset.tsv.bgz", package = "NanoMethViz", mustWork = FALSE)
    mbr <- get_query_test_modbamresult()
    chr <- c("chr7", "chr7", "chr7")
    start <- c(6703892, 6717162, 6703892)
    end <- c(6717161, 6730431, 6717161)

    # test
    out_tabix <- query_methy(methy_path, chr, start, end, simplify = FALSE)
    expect_length(out_tabix, length(chr))
    expect_identical(out_tabix[[1]], out_tabix[[3]])

    out_modbam <- query_methy(mbr, chr, start, end, simplify = FALSE)
    expect_length(out_modbam, length(chr))
    expect_identical(out_modbam[[1]], out_modbam[[3]])

    # duplicated regions should not double count data
    single <- query_methy(mbr, chr[1], start[1], end[1], truncate = FALSE)
    dup <- query_methy(mbr, chr[1:2], start[1:2], end[1:2], truncate = FALSE, simplify = FALSE)
    expect_length(dup, 2)
    expect_identical(dup[[1]], single)
})

test_that("modbam regions with no reads return typed empty output", {
    # setup
    mbr <- get_query_test_modbamresult()
    chr <- c("chr7", "chr7", "chr7")
    start <- c(6703892, 1, 6717162)
    end <- c(6717161, 10, 6730431)

    # test
    out <- query_methy(mbr, chr, start, end, simplify = FALSE)
    expect_length(out, length(chr))
    expect_equal(nrow(out[[2]]), 0)
    expect_true(all(methy_col_names() %in% colnames(out[[2]])))

    out_internal <- query_methy_modbam(mbr, chr[2], start[2], end[2], mod_code(mbr))
    expect_s3_class(out_internal[[1]], "data.frame")
    expect_equal(colnames(out_internal[[1]]), methy_col_names())
    expect_equal(nrow(out_internal[[1]]), 0)
})

test_that("force returns empty output for sequences missing from modbam", {
    # setup
    mbr <- get_query_test_modbamresult()

    # test hard error without force
    expect_warning(
        expect_error(query_methy(mbr, "chrZZ", 1, 10), "no chromosome matches between query and modbam file"),
        "requested sequences missing from modbam file and will be excluded from query:chrZZ"
    )

    # force returns empty output with warning
    expect_warning(
        out <- query_methy(mbr, "chrZZ", 1, 10, force = TRUE),
        "requested sequences missing from modbam file and will be excluded from query:chrZZ"
    )
    expect_equal(nrow(out), 0)
    expect_equal(colnames(out), c(methy_col_names(), "mod_prob"))

    # mixed valid and invalid sequences
    expect_warning(
        out_mixed <- query_methy(mbr, c("chr7", "chrZZ"), c(6703892, 1), c(6717161, 10), force = TRUE, simplify = FALSE),
        "chrZZ"
    )
    expect_length(out_mixed, 2)
    expect_gt(nrow(out_mixed[[1]]), 0)
    expect_equal(nrow(out_mixed[[2]]), 0)
})

test_that("modbam regions on interleaved chromosomes keep their input order", {
    # setup
    mbr <- get_query_test_modbamresult()
    chr <- c("chr7", "chr1", "chr7")
    start <- c(6703892, 1e6, 6717162)
    end <- c(6717161, 2e6, 6730431)

    # test
    out <- query_methy(mbr, chr, start, end, simplify = FALSE, truncate = FALSE)
    expect_length(out, 3)
    expect_equal(nrow(out[[2]]), 0)
    expect_identical(
        out[[1]],
        query_methy(mbr, chr[1], start[1], end[1], truncate = FALSE)
    )
    expect_identical(
        out[[3]],
        query_methy(mbr, chr[3], start[3], end[3], truncate = FALSE)
    )
})

test_that("read_methy_lines parses numeric chromosomes consistently", {
    # setup
    lines_num <- "sample1\t1\t100\t+\t0.5\tread1"
    lines_chr <- "sample1\tchr1\t200\t-\t-0.5\tread2"

    # test
    methy_num <- read_methy_lines(lines_num)
    methy_chr <- read_methy_lines(lines_chr)
    expect_equal(nrow(methy_num), 1)
    expect_false(is.numeric(methy_num$chr))
    expect_true("1" %in% as.character(methy_num$chr))

    # chunks with numeric and character chromosome names combine cleanly
    combined <- dplyr::bind_rows(methy_num, methy_chr)
    expect_equal(nrow(combined), 2)
    expect_equal(as.character(combined$chr), c("1", "chr1"))
})

test_that("site_filter validation rejects invalid values", {
    # setup
    methy_path <- system.file("methy_subset.tsv.bgz", package = "NanoMethViz", mustWork = FALSE)

    # test
    expect_error(
        query_methy(methy_path, "chr7", 6703892, 6717161, site_filter = 0),
        "site_filter must be a single number greater than or equal to 1"
    )
    expect_error(
        query_methy(methy_path, "chr7", 6703892, 6717161, site_filter = "a"),
        "site_filter"
    )
    expect_no_error(query_methy(methy_path, "chr7", 6703892, 6717161, site_filter = 1))
})
