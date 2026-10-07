# helper to build a small exon annotation table
make_exons <- function(transcript_id, start, end, strand = "+", gene_id = transcript_id) {
    tibble::tibble(
        gene_id = gene_id,
        transcript_id = transcript_id,
        symbol = paste0("sym_", gene_id),
        strand = strand,
        start = start,
        end = end
    )
}

# layer indices in plot_gene_annotation() output
layer_gap_pos <- 1
layer_gap_neg <- 3
layer_gap_none <- 5
layer_exons <- 6

test_that("plot_gene_annotation handles empty exons", {
    exons <- make_exons(character(), numeric(), numeric(), strand = character())

    p <- plot_gene_annotation(exons, 0, 1000)

    expect_s3_class(p, "ggplot")
    expect_equal(attr(p, "plot_height"), 0)
})

test_that("plot_gene_annotation handles single-exon transcripts", {
    exons <- make_exons("tx1", 100, 200)

    p <- plot_gene_annotation(exons, 0, 1000)

    expect_s3_class(p, "ggplot")
    # y_offset = 2.5 * 1 - 1 = 1.5, plot_height = 2 + 1.5
    expect_equal(attr(p, "plot_height"), 3.5)
    expect_no_error(ggplot2::ggplot_build(p))
    expect_equal(nrow(ggplot2::layer_data(p, layer_exons)), 1)
    expect_equal(nrow(ggplot2::layer_data(p, layer_gap_none)), 0)
})

test_that("plot_gene_annotation builds multi-exon transcripts on all strands", {
    exons <- dplyr::bind_rows(
        make_exons("tx_pos", c(100, 400, 700), c(200, 500, 800), strand = "+"),
        make_exons("tx_neg", c(150, 450), c(250, 550), strand = "-"),
        make_exons("tx_none", c(300, 600), c(350, 650), strand = "*")
    )

    p <- plot_gene_annotation(exons, 0, 1000)

    expect_s3_class(p, "ggplot")
    # three transcripts: max y_offset = 2.5 * 3 - 1 = 6.5
    expect_equal(attr(p, "plot_height"), 8.5)
    expect_no_error(ggplot2::ggplot_build(p))
    expect_equal(nrow(ggplot2::layer_data(p, layer_gap_pos)), 2)
    expect_equal(nrow(ggplot2::layer_data(p, layer_gap_neg)), 1)
    expect_equal(nrow(ggplot2::layer_data(p, layer_gap_none)), 1)
    expect_equal(nrow(ggplot2::layer_data(p, layer_exons)), 7)
})

test_that("plot_gene_annotation handles exons outside the plot window", {
    exons <- make_exons("tx1", c(100, 400), c(200, 500))

    # summarising the emptied exon table currently warns from min()/max()
    p <- suppressWarnings(plot_gene_annotation(exons, 5000, 6000))

    expect_s3_class(p, "ggplot")
    expect_equal(attr(p, "plot_height"), 0)
    expect_no_error(ggplot2::ggplot_build(p))
})

test_that("plot_gene_annotation computes gaps per transcript for interleaved transcripts", {
    # exons supplied out of order and interleaved between transcripts
    exons <- dplyr::bind_rows(
        make_exons("tx1", c(700, 100, 400), c(800, 200, 500), strand = "*"),
        make_exons("tx2", c(550, 250), c(650, 350), strand = "*")
    )
    exons <- exons[c(3, 4, 1, 5, 2), ]

    p <- plot_gene_annotation(exons, 0, 1000)
    gaps <- ggplot2::layer_data(p, layer_gap_none)
    # connector lines are drawn from gap end (x) to gap start (xend)
    gaps <- gaps[order(gaps$y, gaps$xend), c("xend", "x", "y")]

    # tx1 has y_offset 1.5, tx2 has y_offset 4
    expect_equal(gaps$xend, c(200, 500, 350))
    expect_equal(gaps$x, c(400, 700, 550))
    expect_equal(gaps$y, c(1.5, 1.5, 4) + 0.275)
})

test_that("plot_gene_annotation truncates gaps identically on both strands", {
    starts <- c(100, 300, 900)
    ends <- c(200, 400, 1000)
    plot_start <- 500
    plot_end <- 800

    p_pos <- plot_gene_annotation(make_exons("tx", starts, ends, "+"), plot_start, plot_end)
    p_neg <- plot_gene_annotation(make_exons("tx", starts, ends, "-"), plot_start, plot_end)

    gap_pos <- ggplot2::layer_data(p_pos, layer_gap_pos)
    gap_neg <- ggplot2::layer_data(p_neg, layer_gap_neg)

    # gap 200-300 lies outside the window and is dropped; gap 400-900 spans
    # the window and is clipped to it, with arrows pointing in strand direction
    expect_equal(nrow(gap_pos), 1)
    expect_equal(nrow(gap_neg), 1)
    expect_equal(c(gap_pos$x, gap_pos$xend), c(500, 650))
    expect_equal(c(gap_neg$x, gap_neg$xend), c(800, 650))

    # partially overlapping minus-strand gap is kept and clipped
    p_partial <- plot_gene_annotation(make_exons("tx", c(100, 400), c(200, 500), "-"), 300, 1000)
    gap_partial <- ggplot2::layer_data(p_partial, layer_gap_neg)
    expect_equal(c(gap_partial$x, gap_partial$xend), c(400, 350))
})
