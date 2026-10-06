test_that("Aggregate plotting works", {
    # setup
    nmr <- load_example_nanomethresult()
    mbr <- NanoMethViz:::load_example_modbamresult()
    gene_anno <- exons_to_genes(NanoMethViz::exons(nmr))

    # test
    for (x in list(nmr, mbr)) {
        expect_no_warning(plot_agg_regions(x, gene_anno))
        expect_no_warning(plot_agg_regions(x, gene_anno, group_col = "sample"))
        expect_no_warning(plot_agg_regions(x, gene_anno, group_col = "group"))
    }
})

test_that("Aggregate plotting error checking works", {
    # setup
    nmr <- load_example_nanomethresult()
    gene_anno <- exons_to_genes(NanoMethViz::exons(nmr))

    # test
    expect_error(
        plot_agg_regions(nmr, gene_anno, group_col = "foo"),
        "'foo' could not be found in columns of 'regions' or samples\\(x\\)"
    )
})

test_that("Grouped aggregate plot builds with one line per group", {
    nmr <- load_example_nanomethresult()
    gene_anno <- exons_to_genes(NanoMethViz::exons(nmr))

    p <- plot_agg_regions(nmr, gene_anno, group_col = "group")
    built <- ggplot2::ggplot_build(p)

    n_groups <- length(unique(NanoMethViz::samples(nmr)$group))
    expect_gt(n_groups, 1)
    expect_length(unique(built$data[[1]]$group), n_groups)
})

test_that("Aggregate plot with zero flank labels only start and end", {
    nmr <- load_example_nanomethresult()
    gene_anno <- exons_to_genes(NanoMethViz::exons(nmr))

    p <- plot_agg_regions(nmr, gene_anno, flank = 0)
    built <- ggplot2::ggplot_build(p)

    x_scale <- built$layout$panel_params[[1]]$x
    expect_equal(x_scale$get_labels(), c("start", "end"))
    # region boundary lines are only drawn when there are flanks
    expect_length(p$layers, 1)
})

test_that("Aggregate plot with flank labels flanks in kb", {
    nmr <- load_example_nanomethresult()
    gene_anno <- exons_to_genes(NanoMethViz::exons(nmr))

    p <- plot_agg_regions(nmr, gene_anno, flank = 1500)
    built <- ggplot2::ggplot_build(p)

    x_scale <- built$layout$panel_params[[1]]$x
    expect_equal(x_scale$get_labels(), c("-1.5kb", "start", "end", "+1.5kb"))
    expect_length(p$layers, 3)
})

test_that("Unstranded aggregate plot does not flip negative strand regions", {
    nmr <- load_example_nanomethresult()
    gene_anno <- exons_to_genes(NanoMethViz::exons(nmr))
    # flipping only matters if there are negative strand regions
    expect_true(any(gene_anno$strand == "-"))

    p_stranded <- plot_agg_regions(nmr, gene_anno, stranded = TRUE)
    p_unstranded <- plot_agg_regions(nmr, gene_anno, stranded = FALSE)
    expect_no_error(ggplot2::ggplot_build(p_unstranded))

    expect_false(isTRUE(all.equal(
        p_stranded$layers[[1]]$data,
        p_unstranded$layers[[1]]$data
    )))
})

test_that("Binary threshold controls methylation calls", {
    nmr <- load_example_nanomethresult()
    gene_anno <- exons_to_genes(NanoMethViz::exons(nmr))

    p_low <- plot_agg_regions(nmr, gene_anno, binary_threshold = 0.1)
    p_high <- plot_agg_regions(nmr, gene_anno, binary_threshold = 0.9)
    expect_no_error(ggplot2::ggplot_build(p_high))

    prop_low <- p_low$layers[[1]]$data$methy_prop
    prop_high <- p_high$layers[[1]]$data$methy_prop
    expect_true(all(prop_low >= prop_high))
    expect_gt(mean(prop_low), mean(prop_high))
})
