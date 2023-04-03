test_that("plot_pca works correctly", {
  # create test data
    nmr <- load_example_nanomethresult()
    bss <- methy_to_bsseq(nmr)
    x <- bsseq_to_log_methy_ratio(bss)
    groups <- samples(nmr)$group
    labels <- colnames(x)

    expect_no_error(plot_pca(x))
    expect_no_error(plot_pca(x, labels = labels, groups = groups))
    expect_error(plot_pca(x, groups = c(groups, "foo")))
    expect_error(plot_pca(x, labels = c(labels, "foo")))
})
