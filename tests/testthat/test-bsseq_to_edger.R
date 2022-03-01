test_that("Conversion works", {
    nmr <- load_example_nanomethresult()
    bsseq <- methy_to_bsseq(methy(nmr))
    regions <- exons_to_genes(NanoMethViz::exons(nmr))

    edger_site <- expect_silent(bsseq_to_edger(bsseq))
    expect_ncol(edger_site, 2*ncol(bsseq))
    expect_nrow(edger_site, nrow(bsseq))

    edger_region <- expect_silent(bsseq_to_edger(bsseq, regions))
    expect_ncol(edger_region, 2*ncol(bsseq))
    expect_nrow(edger_region, nrow(regions))

    lmr_regions <- expect_silent(bsseq_to_log_methy_ratio(bsseq, regions))
    expect_ncol(lmr_regions, ncol(bsseq))
    expect_nrow(lmr_regions, nrow(regions))
})
