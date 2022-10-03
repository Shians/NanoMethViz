test_that("NanoMethResults getters work", {
    # setup
    nmr <- load_example_nanomethresult()

    # test
    expect_true(is(NanoMethViz::methy(nmr), "character"))
    expect_true(fs::file_exists(NanoMethViz::methy(nmr)))
    expect_true(is(NanoMethViz::exons(nmr), "data.frame"))
    expect_true(is(NanoMethViz::samples(nmr), "data.frame"))

    expect_silent(
        NanoMethResult(
            NanoMethViz::methy(nmr),
            NanoMethViz::samples(nmr)
        )
    )

    expect_silent(methy(nmr) <- methy(nmr))
    expect_error(methy(nmr) <- "invalid_path")
    expect_silent(samples(nmr) <- samples(nmr))
    expect_silent(exons(nmr) <- exons(nmr))
    expect_error(exons(nmr) <- exons(nmr)[, -"strand"])
})
