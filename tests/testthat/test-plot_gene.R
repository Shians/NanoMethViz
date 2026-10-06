empty_exons <- function() {
    tibble::tibble(
        gene_id = character(0),
        chr = character(0),
        strand = character(0),
        start = integer(0),
        end = integer(0),
        transcript_id = character(0),
        symbol = character(0)
    )
}

test_that("plot_gene warns when deprecated span argument is used", {
    nmr <- load_example_nanomethresult()

    expect_warning(
        plot_gene(nmr, "Peg3", heatmap = FALSE, span = 0.5),
        "the 'span' argument has been deprecated"
    )
})

test_that("plot_gene errors on empty or NA gene symbol", {
    nmr <- load_example_nanomethresult()
    mbr <- load_example_modbamresult()

    for (x in list(nmr, mbr)) {
        expect_error(plot_gene(x, ""), "Gene symbol cannot be empty or NA")
        expect_error(plot_gene(x, NA_character_), "Gene symbol cannot be empty or NA")
    }
})

test_that("plot_gene errors when exon annotation is empty", {
    nmr <- load_example_nanomethresult()
    mbr <- load_example_modbamresult()
    nmr@exons <- empty_exons()
    mbr@exons <- empty_exons()

    for (x in list(nmr, mbr)) {
        expect_error(plot_gene(x, "Peg3"), "No exon annotation found in the data object")
    }
})

test_that("plot_gene suggests similar genes when gene is not found", {
    nmr <- load_example_nanomethresult()
    mbr <- load_example_modbamresult()

    for (x in list(nmr, mbr)) {
        expect_error(plot_gene(x, "Pegx"), "Gene 'Pegx' not found in exon annotation")
        expect_error(plot_gene(x, "Pegx"), "Similar genes found: Peg3", fixed = TRUE)
        expect_error(plot_gene(x, "pegx"), "Similar genes found: Peg3", fixed = TRUE)
        expect_error(plot_gene(x, "Pegx"), "Available genes:", fixed = TRUE)
    }
})

test_that("plot_gene lists available genes when no similar gene exists", {
    nmr <- load_example_nanomethresult()

    err <- expect_error(plot_gene(nmr, "Zzz1"), "Available genes: ")
    expect_no_match(conditionMessage(err), "Similar genes found")
    expect_match(conditionMessage(err), "Peg3", fixed = TRUE)
})

test_that("plot_gene treats gene prefix literally rather than as a regex", {
    nmr <- load_example_nanomethresult()

    err <- expect_error(plot_gene(nmr, ".*x"), "Gene '.*x' not found", fixed = TRUE)
    expect_no_match(conditionMessage(err), "Similar genes found")
})
