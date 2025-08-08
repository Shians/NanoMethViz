test_that("Plotting gene methylation heatmap works", {
    # setup
    nmr <- load_example_nanomethresult()

    # test
    p <- expect_no_warning(plot_gene_heatmap(nmr, "Peg3"))
    expect_s3_class(p, "ggplot")

    # test for a bug whereby samples not present in data cause function to hang
    nmr_extra_sample <- load_example_nanomethresult()
    samples(nmr_extra_sample) <- bind_rows(
        samples(nmr_extra_sample),
        c(sample = "foo", group = "bar")
    )
    p <- expect_no_warning(plot_gene_heatmap(nmr_extra_sample, "Peg3"))
    expect_s3_class(p, "ggplot")
})

test_that("plot_gene_heatmap works with different parameters", {
    nmr <- load_example_nanomethresult()

    # test with different window_prop values
    p <- expect_no_warning(plot_gene_heatmap(nmr, "Peg3", window_prop = 0.5))
    expect_s3_class(p, "ggplot")

    # test with vector window_prop
    p <- expect_no_warning(plot_gene_heatmap(nmr, "Peg3", window_prop = c(0.2, 0.4)))
    expect_s3_class(p, "ggplot")

    # test with compact pos_style
    p <- expect_no_warning(plot_gene_heatmap(nmr, "Peg3", pos_style = "compact"))
    expect_s3_class(p, "ggplot")

    # test with different subsample value
    p <- expect_no_warning(plot_gene_heatmap(nmr, "Peg3", subsample = 25))
    expect_s3_class(p, "ggplot")
})

test_that("plot_gene_heatmap error handling", {
    nmr <- load_example_nanomethresult()

    # test with empty exons
    nmr_no_exons <- nmr
    nmr_no_exons@exons <- tibble::tibble(
        gene_id = character(0),
        chr = character(0),
        strand = character(0),
        start = integer(0),
        end = integer(0),
        transcript_id = character(0),
        symbol = character(0)
    )

    expect_error(
        plot_gene_heatmap(nmr_no_exons, "Peg3"),
        "exons\\(x\\) is empty, gene cannot be queried"
    )

    # test with gene not in annotation
    expect_error(
        plot_gene_heatmap(nmr, "NonExistentGene"),
        "gene NonExistentGene not found in exon annotation"
    )
})
