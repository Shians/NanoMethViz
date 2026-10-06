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

# Labels as resolved by ggplot2; get_labs() exists from ggplot2 3.5.2 onwards.
plot_labs <- function(p) {
    if (exists("get_labs", envir = asNamespace("ggplot2"))) {
        ggplot2::get_labs(p)
    } else {
        p$labels
    }
}

# Text drawn in the rendered legend, independent of how the title was set.
legend_text <- function(p) {
    gtable <- ggplot2::ggplotGrob(p)
    grob_fields <- unlist(gtable$grobs[grepl("guide-box", gtable$layout$name)])
    unname(grob_fields[grepl("label$", names(grob_fields))])
}

test_that("plot_pca handles label and group combinations", {
    nmr <- load_example_nanomethresult()
    bss <- methy_to_bsseq(nmr)
    x <- bsseq_to_log_methy_ratio(bss)
    groups <- samples(nmr)$group

    # no labels, no groups
    p <- plot_pca(x, labels = NULL)
    expect_s3_class(p, "ggplot")
    expect_s3_class(p$layers[[1]]$geom, "GeomPoint")
    expect_no_error(ggplot2::ggplot_build(p))

    # no labels, character groups with custom legend name
    p <- plot_pca(x, labels = NULL, groups = groups, legend_name = "condition")
    expect_s3_class(p$layers[[1]]$geom, "GeomPoint")
    built <- ggplot2::ggplot_build(p)
    expect_true(built$plot$scales$get_scales("colour")$is_discrete())
    expect_contains(legend_text(p), "condition")

    # numeric groups with labels ignores labels
    expect_message(
        p <- plot_pca(x, groups = 1:6, legend_name = "score"),
        "Ignoring labels"
    )
    expect_s3_class(p$layers[[1]]$geom, "GeomPoint")
    built <- ggplot2::ggplot_build(p)
    colour_scale <- built$plot$scales$get_scales("colour")
    expect_false(colour_scale$is_discrete())
    expect_equal(colour_scale$name, "score")

    # numeric groups without labels is silent
    expect_no_message(p <- plot_pca(x, labels = NULL, groups = 1:6))
    expect_no_error(ggplot2::ggplot_build(p))
})

test_that("plot_pca axis labels follow plot_dims", {
    nmr <- load_example_nanomethresult()
    bss <- methy_to_bsseq(nmr)
    x <- bsseq_to_log_methy_ratio(bss)

    labs <- plot_labs(plot_pca(x))
    expect_equal(as.character(labs$x), "PCA Dim 1")
    expect_equal(as.character(labs$y), "PCA Dim 2")

    p <- plot_pca(x, plot_dims = c(2, 3))
    expect_no_error(ggplot2::ggplot_build(p))
    labs <- plot_labs(p)
    expect_equal(as.character(labs$x), "PCA Dim 2")
    expect_equal(as.character(labs$y), "PCA Dim 3")
})
