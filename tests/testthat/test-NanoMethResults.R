test_that("NanoMethResults getters work", {
    # setup
    nmr <- load_example_nanomethresult()

    # test
    expect_true(is(NanoMethViz::methy(nmr), "character"))
    expect_true(fs::file_exists(NanoMethViz::methy(nmr)))

    expect_true(is(NanoMethViz::exons(nmr), "data.frame"))

    expect_true(is(NanoMethViz::samples(nmr), "data.frame"))
})
