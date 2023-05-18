test_that("Querying methylation works", {
    # setup
    nmr <- load_example_nanomethresult()

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
})
