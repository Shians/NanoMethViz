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
