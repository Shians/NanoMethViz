test_that("query_exons_region works with NanoMethResult", {
    nmr <- load_example_nanomethresult()

    # Test basic region query
    result <- query_exons_region(nmr, "chr7", 6703000, 6730000)
    expect_s3_class(result, "tbl_df")
    expect_true(nrow(result) > 0)
    expect_true(all(c("gene_id", "chr", "start", "end", "symbol") %in% colnames(result)))
    expect_true(all(result$chr == "chr7"))

    # Test region with no overlapping genes
    result_empty <- query_exons_region(nmr, "chr7", 1, 100)
    expect_s3_class(result_empty, "tbl_df")
    expect_equal(nrow(result_empty), 0)

    # Test region with partial overlap
    result_partial <- query_exons_region(nmr, "chr7", 6720000, 6725000)
    expect_s3_class(result_partial, "tbl_df")
})

test_that("query_exons_region works with data.frame", {
    nmr <- load_example_nanomethresult()
    exons_df <- exons(nmr)

    # Test with data.frame directly
    result <- query_exons_region(exons_df, "chr7", 6703000, 6730000)
    expect_s3_class(result, "tbl_df")
    expect_true(nrow(result) > 0)
    expect_true(all(result$chr == "chr7"))
})

test_that("query_exons_gene_id works correctly", {
    nmr <- load_example_nanomethresult()

    # Get a valid gene_id from the exons
    valid_gene_id <- exons(nmr)$gene_id[1]

    # Test with valid gene_id
    result <- query_exons_gene_id(nmr, valid_gene_id)
    expect_s3_class(result, "tbl_df")
    expect_true(nrow(result) > 0)
    expect_true(all(result$gene_id == valid_gene_id))

    # Test with multiple gene_ids
    gene_ids <- exons(nmr)$gene_id[1:2]
    result_multiple <- query_exons_gene_id(nmr, gene_ids)
    expect_s3_class(result_multiple, "tbl_df")
    expect_true(nrow(result_multiple) > 0)
    expect_true(all(result_multiple$gene_id %in% gene_ids))

    # Test with data.frame
    exons_df <- exons(nmr)
    result_df <- query_exons_gene_id(exons_df, valid_gene_id)
    expect_s3_class(result_df, "tbl_df")
    expect_true(all(result_df$gene_id == valid_gene_id))
})

test_that("query_exons_gene_id error handling", {
    nmr <- load_example_nanomethresult()

    # Test with non-existent gene_id
    expect_error(
        query_exons_gene_id(nmr, "NonExistentGeneID"),
        "gene NonExistentGeneID not found in exon annotation"
    )

    # Test with vector of non-existent gene_ids
    expect_error(
        query_exons_gene_id(nmr, c("Gene1", "Gene2")),
        "not found in exon annotation"
    )
})

test_that("query_exons_symbol works correctly", {
    nmr <- load_example_nanomethresult()

    # Get a valid symbol from the exons
    valid_symbol <- exons(nmr)$symbol[1]

    # Test with valid symbol
    result <- query_exons_symbol(nmr, valid_symbol)
    expect_s3_class(result, "tbl_df")
    expect_true(nrow(result) > 0)
    expect_true(all(result$symbol == valid_symbol))

    # Test with multiple symbols
    symbols <- unique(exons(nmr)$symbol)[1:2]
    result_multiple <- query_exons_symbol(nmr, symbols)
    expect_s3_class(result_multiple, "tbl_df")
    expect_true(nrow(result_multiple) > 0)
    expect_true(all(result_multiple$symbol %in% symbols))

    # Test with data.frame
    exons_df <- exons(nmr)
    result_df <- query_exons_symbol(exons_df, valid_symbol)
    expect_s3_class(result_df, "tbl_df")
    expect_true(all(result_df$symbol == valid_symbol))
})

test_that("query_exons_symbol error handling", {
    nmr <- load_example_nanomethresult()

    # Test with non-existent symbol
    expect_error(
        query_exons_symbol(nmr, "NonExistentSymbol"),
        "gene NonExistentSymbol not found in exon annotation"
    )

    # Test with vector of non-existent symbols
    expect_error(
        query_exons_symbol(nmr, c("Symbol1", "Symbol2")),
        "not found in exon annotation"
    )
})

test_that("get_exons helper function works correctly", {
    nmr <- load_example_nanomethresult()
    exons_df <- exons(nmr)

    # Test with NanoMethResult
    result_nmr <- NanoMethViz:::get_exons(nmr)
    expect_identical(result_nmr, exons(nmr))

    # Test with data.frame
    result_df <- NanoMethViz:::get_exons(exons_df)
    expect_identical(result_df, exons_df)

    # Test error handling with invalid object
    expect_error(
        NanoMethViz:::get_exons("invalid_object"),
        "x must be a NanoMethResult, ModBamResult or data.frame"
    )

    expect_error(
        NanoMethViz:::get_exons(123),
        "x must be a NanoMethResult, ModBamResult or data.frame"
    )
})

test_that("query_exons functions work with edge cases", {
    nmr <- load_example_nanomethresult()

    # Test query_exons_region with exact boundaries
    first_gene <- exons(nmr) %>%
        dplyr::group_by(gene_id) %>%
        dplyr::summarise(start = min(start), end = max(end), chr = first(chr)) %>%
        dplyr::slice(1)

    result_exact <- query_exons_region(nmr, first_gene$chr, first_gene$start, first_gene$end)
    expect_true(nrow(result_exact) > 0)
    expect_true(first_gene$gene_id %in% result_exact$gene_id)

    # Test query_exons_region with chromosome that doesn't exist
    result_no_chr <- query_exons_region(nmr, "chrNonExistent", 1, 1000)
    expect_equal(nrow(result_no_chr), 0)

    # Test with empty exons data.frame
    empty_exons <- tibble::tibble(
        gene_id = character(0),
        chr = character(0),
        strand = character(0),
        start = integer(0),
        end = integer(0),
        transcript_id = character(0),
        symbol = character(0)
    )

    expect_warning(
        result_empty_df <- query_exons_region(empty_exons, "chr1", 1, 1000),
        "no non-missing arguments"
    )
    expect_equal(nrow(result_empty_df), 0)

    expect_error(
        query_exons_gene_id(empty_exons, "AnyGene"),
        "gene AnyGene not found in exon annotation"
    )

    expect_error(
        query_exons_symbol(empty_exons, "AnySymbol"),
        "gene AnySymbol not found in exon annotation"
    )
})
