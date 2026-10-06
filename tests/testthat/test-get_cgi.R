test_that("CGI getters work", {
    cgi_getters <- list(
        get_cgi_hg19,
        get_cgi_hg38,
        get_cgi_t2t,
        get_cgi_mm10,
        get_cgi_grcm39
    )

    for (cgi_fn in cgi_getters) {
        cgi_anno <- expect_no_warning(cgi_fn())
        expect_s3_class(cgi_anno, "data.frame")
        expect_gt(nrow(cgi_anno), 0)
        expect_gt(ncol(cgi_anno), 0)
        expect_contains(colnames(cgi_anno), c("gene_id", "chr", "strand", "start", "end", "transcript_id", "symbol"))
        expect_true(all(startsWith(cgi_anno$gene_id, "CpG:")))
    }
})

test_that("get_cgi rejects unknown genomes", {
    expect_error(get_cgi("not_a_genome"), "genome must be one of")
})
