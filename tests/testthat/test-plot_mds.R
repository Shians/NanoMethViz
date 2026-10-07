test_that("Plot MDS works", {
    nmr <- load_example_nanomethresult()
    bss <- methy_to_bsseq(nmr)
    lmr <- bsseq_to_log_methy_ratio(bss)

    expect_no_warning(plot_mds(lmr))

    expect_no_warning(plot_mds(lmr, plot_dims = c(2, 3)))
    expect_no_warning(plot_mds(lmr, labels = paste0("samples", 1:6)))
    expect_no_warning(plot_mds(lmr, labels = paste0("samples", 1:6), groups = rep(c("A", "B"), 3)))

    expect_no_warning(plot_mds(lmr, groups = rep(c("A", "B"), 3)))
    expect_no_warning(plot_mds(lmr, groups = rep(c("A", "B"), 3), legend_name = "group_name"))

    expect_message(plot_mds(lmr, groups = 1:6), "Ignoring labels as groups is numeric. Set `labels=NULL` to suppress this message.")
    expect_no_warning(plot_mds(lmr, groups = 1:6, labels = NULL))
})

test_that("Plot MDS coordinates match limma plotMDS output", {
    nmr <- load_example_nanomethresult()
    bss <- methy_to_bsseq(nmr)
    lmr <- bsseq_to_log_methy_ratio(bss)

    for (plot_dims in list(c(1, 2), c(2, 3))) {
        mds <- limma::plotMDS(lmr, top = 500, dim.plot = plot_dims, plot = FALSE)
        p <- plot_mds(lmr, plot_dims = plot_dims)
        layer_data <- ggplot2::ggplot_build(p)$data[[1]]

        expect_equal(p$data$dim1, mds$x)
        expect_equal(p$data$dim2, mds$y)
        expect_equal(layer_data$x, mds$x)
        expect_equal(layer_data$y, mds$y)
    }
})

# Text drawn in the rendered legend, independent of how the title was set.
legend_text <- function(p) {
    gtable <- ggplot2::ggplotGrob(p)
    grob_fields <- unlist(gtable$grobs[grepl("guide-box", gtable$layout$name)])
    unname(grob_fields[grepl("label$", names(grob_fields))])
}

test_that("Plot MDS handles label and group combinations", {
    nmr <- load_example_nanomethresult()
    bss <- methy_to_bsseq(nmr)
    lmr <- bsseq_to_log_methy_ratio(bss)

    # no labels, no groups
    p <- plot_mds(lmr, labels = NULL)
    expect_s3_class(p, "ggplot")
    expect_s3_class(p$layers[[1]]$geom, "GeomPoint")
    expect_no_error(ggplot2::ggplot_build(p))

    # no labels, discrete groups with custom legend name
    p <- plot_mds(lmr, labels = NULL, groups = rep(c("A", "B"), 3), legend_name = "condition")
    expect_s3_class(p$layers[[1]]$geom, "GeomPoint")
    built <- ggplot2::ggplot_build(p)
    expect_true(built$plot$scales$get_scales("colour")$is_discrete())
    expect_contains(legend_text(p), "condition")
})

test_that("Plot MDS rejects labels and groups of the wrong length", {
    nmr <- load_example_nanomethresult()
    bss <- methy_to_bsseq(nmr)
    lmr <- bsseq_to_log_methy_ratio(bss)

    expect_error(plot_mds(lmr, labels = paste0("samples", 1:5)))
    expect_error(plot_mds(lmr, labels = paste0("samples", 1:7)))
    expect_error(plot_mds(lmr, groups = rep(c("A", "B"), 2)))
    expect_error(plot_mds(lmr, groups = 1:7))
})
