test_that("Stacking works", {
    # setup
    nmr <- load_example_nanomethresult()
    methy <- query_methy_gene(nmr, gene = "Peg3")
    reads <- methy %>%
        dplyr::group_by(read_name) %>%
        dplyr::summarise(start = min(pos), end = max(pos))

    stacked_intervals <- stacked_intervals(reads)
    stacked_interval_inds <- stacked_interval_inds(reads)

    # test
    expect_equal(nrow(stacked_intervals), nrow(reads))
    expect_lte(dplyr::n_distinct(stacked_intervals$group), nrow(reads))
    expect_s3_class(stacked_interval_inds, "factor")
})
