test_that("Plotting gene aggregates works", {
    nmr <- load_example_nanomethresult()
    expect_no_error(plot_agg_genes(nmr))
    expect_no_error(plot_agg_tss(nmr))
    expect_no_error(plot_agg_tes(nmr))

    # empty exons should throw error
    exons(nmr) <- tibble::tibble(
        gene_id = character(),
        chr = character(),
        strand = character(),
        start = integer(),
        end = integer(),
        transcript_id = character(),
        symbol = character()
    )
    expect_error(plot_agg_genes(nmr), "no exon annotations found in object")
    expect_error(plot_agg_tss(nmr), "no exon annotations found in object")
    expect_error(plot_agg_tes(nmr), "no exon annotations found in object")
})

test_that("Gene aggregate plots can be subset to selected genes", {
    nmr <- load_example_nanomethresult()
    subset_genes <- c("Peg3", "Impact")
    expect_true(all(subset_genes %in% NanoMethViz::exons(nmr)$symbol))

    # a wider span avoids loess fit failures on the sparse TSS/TES windows
    for (plot_fn in list(plot_agg_genes, plot_agg_tss, plot_agg_tes)) {
        p_all <- plot_fn(nmr, span = 0.2)
        p_subset <- plot_fn(nmr, genes = subset_genes, span = 0.2)
        expect_no_warning(ggplot2::ggplot_build(p_subset))

        expect_false(isTRUE(all.equal(
            p_all$layers[[1]]$data,
            p_subset$layers[[1]]$data
        )))
    }
})

test_that("TSS and TES aggregate plots are labelled around the site", {
    nmr <- load_example_nanomethresult()

    built_tss <- ggplot2::ggplot_build(plot_agg_tss(nmr, span = 0.2))
    built_tes <- ggplot2::ggplot_build(plot_agg_tes(nmr, span = 0.2))

    expect_equal(
        built_tss$layout$panel_params[[1]]$x$get_labels(),
        c("-2kb", "TSS", "+2kb")
    )
    expect_equal(
        built_tes$layout$panel_params[[1]]$x$get_labels(),
        c("-2kb", "TES", "+2kb")
    )
})

test_that("Gene aggregate plots error when no genes match", {
    nmr <- load_example_nanomethresult()

    # currently fails inside dplyr::inner_join() rather than with an
    # informative message about the unmatched gene names
    expect_error(plot_agg_genes(nmr, genes = "not_a_gene"))
    expect_error(plot_agg_tss(nmr, genes = "not_a_gene"))
    expect_error(plot_agg_tes(nmr, genes = "not_a_gene"))
})
