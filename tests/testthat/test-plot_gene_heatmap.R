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
        "No exon annotation found in the data object"
    )

    # test with gene not in annotation
    expect_error(
        plot_gene_heatmap(nmr, "NonExistentGene"),
        "Gene 'NonExistentGene' not found in exon annotation"
    )
})

test_that("plot_gene_heatmap errors on empty or NA gene symbol", {
    nmr <- load_example_nanomethresult()

    expect_error(plot_gene_heatmap(nmr, ""), "Gene symbol cannot be empty or NA")
    expect_error(plot_gene_heatmap(nmr, NA_character_), "Gene symbol cannot be empty or NA")
})

test_that("plot_gene_heatmap gives helpful hints for unknown genes", {
    nmr <- load_example_nanomethresult()

    # similar genes hint, with lines separated by real newlines
    err <- expect_error(plot_gene_heatmap(nmr, "Pegx"), "Similar genes found: Peg3", fixed = TRUE)
    expect_match(conditionMessage(err), "annotation.\nPlease check", fixed = TRUE)
    expect_no_match(conditionMessage(err), "\\n", fixed = TRUE)

    # available genes hint when no prefix match
    err <- expect_error(plot_gene_heatmap(nmr, "Zzz1"), "Available genes: ")
    expect_no_match(conditionMessage(err), "Similar genes found")
})

test_that("plot_gene_heatmap points to full symbol list when many genes are annotated", {
    nmr <- load_example_nanomethresult()
    extra_exons <- exons(nmr)[rep(1, 30), ]
    extra_exons$symbol <- sprintf("Gene%02d", 1:30)
    extra_exons$gene_id <- sprintf("id%02d", 1:30)
    nmr@exons <- dplyr::bind_rows(exons(nmr), extra_exons)

    err <- expect_error(
        plot_gene_heatmap(nmr, "Zzz1"),
        "Use unique(exons(your_object)$symbol) to see all 36 available genes.",
        fixed = TRUE
    )
    expect_no_match(conditionMessage(err), "Available genes:")

    # more than 10 similar genes are truncated and the count hint is still given
    err <- expect_error(plot_gene_heatmap(nmr, "Genx"), "Similar genes found: Gene01")
    expect_match(conditionMessage(err), "(20 more)", fixed = TRUE)
    expect_match(conditionMessage(err), "see all 36 available genes", fixed = TRUE)
})
