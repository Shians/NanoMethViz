test_that("methy_to_bsseq works", {
    # setup
    nmr <- load_example_nanomethresult()
    bss <- methy_to_bsseq(nmr)

    # test
    expect_no_error(methy_to_bsseq(methy(nmr)))
    expect_s4_class(methy_to_bsseq(nmr), "BSseq")
    expect_equal(ncol(bss), 6)
    expect_equal(
        nrow(SummarizedExperiment::colData(bss)),
        nrow(NanoMethViz::samples(nmr))
    )
    expect_equal(
        ncol(SummarizedExperiment::colData(bss)),
        ncol(NanoMethViz::samples(nmr))
    )
    expect_equal(
        colnames(SummarizedExperiment::colData(bss)),
        colnames(NanoMethViz::samples(nmr))
    )
})

test_that("repeated methy_to_bsseq calls don't share state", {
    methy_file <- system.file("methy_subset.tsv.bgz", package = "NanoMethViz", mustWork = FALSE)

    out_folder1 <- fs::path(tempdir(), "dss-1")
    out_folder2 <- fs::path(tempdir(), "dss-2")
    fs::dir_create(out_folder1)
    fs::dir_create(out_folder2)

    bss1 <- suppressMessages(methy_to_bsseq(methy_file, out_folder1, verbose = FALSE))
    bss2 <- suppressMessages(methy_to_bsseq(methy_file, out_folder2, verbose = FALSE))

    expect_equal(ncol(bss2), ncol(bss1))
    expect_equal(
        sort(unique(bsseq::getBSseq(bss2, "M"))),
        sort(unique(bsseq::getBSseq(bss1, "M")))
    )

    # second call's intermediate files must contain exactly one header
    dss_files <- fs::dir_ls(out_folder2, glob = "*.txt")
    for (f in dss_files) {
        lines <- readLines(f)
        expect_equal(lines[1], "chr\tpos\ttotal\tmethylated")
        expect_false(any(lines[-1] == "chr\tpos\ttotal\tmethylated"))
    }
})

test_that("bsseq_to_* works", {
    nmr <- load_example_nanomethresult()
    bss <- methy_to_bsseq(NanoMethViz::methy(nmr))

    # test
    edger_counts <- expect_no_warning(bsseq_to_edger(bss))
    expect_ncol(edger_counts, 12)
    expect_nrow(edger_counts, 4778)
    lmr <- expect_no_warning(bsseq_to_log_methy_ratio(bss))
    expect_ncol(lmr, 6)
    expect_nrow(lmr, 4778)

    expect_equal(ncol(edger_counts), 2 * ncol(bss))
    expect_equal(ncol(lmr), ncol(bss))
})
