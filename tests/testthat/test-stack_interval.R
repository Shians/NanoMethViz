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

    # assertion 1: the number of unique groups should be less than or equal to the number of rows
    expect_true(length(unique(stacked_intervals$group)) <= nrow(reads),
        info = "the number of unique groups should be less than or equal to the number of rows")

    # assertion 2: reads within group should not overlap
    for (i in unique(stacked_intervals$group)) {
        group_reads <- stacked_intervals %>%
            dplyr::filter(group == i) %>%
            dplyr::arrange(start)

        if (nrow(group_reads) == 1) {
            next
        }

        for (j in 1:(nrow(group_reads) - 1)) {
            expect_true(group_reads$end[j] < group_reads$start[j+1],
                info = glue::glue("no reads overlap within each group: {group_reads$end[j]} > {group_reads$start[j+1]}"))
        }
    }
})
