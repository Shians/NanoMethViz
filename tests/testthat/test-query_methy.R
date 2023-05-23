test_that("Querying methylation works", {
    # setup
    nmr <- load_example_nanomethresult()
    mbr <- ModBamResult(
        methy = ModBamFiles(
            paths = system.file(package = "NanoMethViz", "peg3.bam"),
            samples = "sample1"
        ),
        samples = tibble::tibble(
            sample = "sample1",
            group = "group1"
        )
    )

    # test
    expect_silent(methy_data1 <- query_methy(methy(nmr), "chr7", 6703892, 6730431))
    expect_equal(
        colnames(methy_data1),
        c("sample", "chr", "pos", "strand", "statistic", "read_name", "mod_prob")
    )

    expect_silent(methy_data2 <- query_methy(nmr, "chr7", 6703892, 6730431))
    expect_equal(methy_data1, methy_data2)


    expect_silent(gene_methy_data <- query_methy_gene(nmr, "Peg3"))
    expect_equal(methy_data1, gene_methy_data)

    expect_warning(expect_error(query_methy(nmr, "Missing", 1, 1000)))
    expect_error(query_methy_gene(nmr, "Missing"))

    regions <- data.frame(
        chr = c("chr7", "chr7"),
        start = c(6703892, 6717161),
        end = c(6717162, 6730431)
    )

    expect_silent(query_methy_df(methy(nmr), regions))

    regions_gr <- GenomicRanges::GRanges(regions)

    query_methy_gr(methy(nmr), regions_gr)

    expect_silent(query_methy(mbr, "chr7", 6703892, 6730431))
})
