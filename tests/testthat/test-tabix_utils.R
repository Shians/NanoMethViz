test_that("Importing nanopolish works", {
    # setup
    methy_calls <- system.file(package = "NanoMethViz",
        c("sample1_nanopolish.tsv.gz", "sample2_nanopolish.tsv.gz"), mustWork = FALSE)
    temp_file <- paste0(tempfile(), ".tsv.bgz")
    withr::defer(file.remove(temp_file))

    # test
    expect_message(create_tabix_file(methy_calls, temp_file))
    expect_s4_class(methy_to_bsseq(temp_file), "BSseq")
})

test_that("Importing megalodon works", {
    # setup
    methy_calls <- system.file(package = "NanoMethViz",
       "megalodon_calls.txt.gz", mustWork = FALSE)
    temp_file <- paste0(tempfile(), ".tsv.bgz")
    withr::defer(file.remove(temp_file))

    # test
    expect_message(create_tabix_file(methy_calls, temp_file))
})

test_that("create_tabix_file rejects mod_code for non-modkit input", {
    methy_calls <- system.file(package = "NanoMethViz",
        "sample1_nanopolish.tsv.gz", mustWork = FALSE)
    temp_file <- paste0(tempfile(), ".tsv.bgz")

    expect_error(
        create_tabix_file(methy_calls, temp_file, mod_code = "m", verbose = FALSE),
        "mod_code only applies to modkit input"
    )
    expect_false(file.exists(temp_file))
})

test_that("sort_methy_file stops when the sort command fails", {
    skip_on_os("windows")

    x <- tempfile()
    writeLines(c("sample1\tchr1\t1\t+\t0.5\tread1"), x)
    withr::defer(unlink(x))

    Sys.chmod(x, "0000")
    expect_error(sort_methy_file(x), "Failed to sort methylation file")
})

test_that("create_tabix_file defaults to 5mC and threads mod_code to modkit input", {
    modkit_calls <- tempfile(fileext = ".tsv")
    write_modkit_call <- function(read_id, ref_position, mod_code, mod_qual) {
        tibble::tibble(
            read_id = read_id,
            forward_read_position = ref_position,
            ref_position = ref_position,
            chrom = "chr1",
            mod_strand = "+",
            ref_strand = "+",
            ref_mod_strand = "+",
            fw_soft_clipped_start = 0,
            fw_soft_clipped_end = 100,
            read_length = 100,
            mod_qual = mod_qual,
            mod_code = mod_code,
            base_qual = 30,
            ref_kmer = "ACG",
            query_kmer = "ACG",
            canonical_base = "C",
            modified_primary_base = "m",
            inferred = FALSE,
            flag = 0
        )
    }
    modkit_df <- dplyr::bind_rows(
        write_modkit_call("read1", 9, "m", 0.9),
        write_modkit_call("read2", 19, "h", 0.5),
        write_modkit_call("read3", 29, "a", 0.8)
    )
    readr::write_tsv(modkit_df, modkit_calls)
    withr::defer(unlink(modkit_calls))

    temp_file <- paste0(tempfile(), ".tsv.bgz")
    withr::defer({
        unlink(temp_file)
        unlink(paste0(temp_file, ".tbi"))
    })

    expect_message(
        create_tabix_file(modkit_calls, temp_file, samples = "sample1")
    )
    methy <- readr::read_tsv(
        gzfile(temp_file),
        col_names = methy_col_names(),
        col_types = methy_col_types()
    )
    expect_equal(nrow(methy), 1)
    expect_equal(methy$pos, 10L)

    temp_file_h <- paste0(tempfile(), ".tsv.bgz")
    withr::defer({
        unlink(temp_file_h)
        unlink(paste0(temp_file_h, ".tbi"))
    })
    create_tabix_file(modkit_calls, temp_file_h, samples = "sample1", mod_code = "h", verbose = FALSE)
    methy_h <- readr::read_tsv(
        gzfile(temp_file_h),
        col_names = methy_col_names(),
        col_types = methy_col_types()
    )
    expect_equal(nrow(methy_h), 1)
    expect_equal(methy_h$pos, 20L)

    temp_file_none <- paste0(tempfile(), ".tsv.bgz")
    expect_error(
        create_tabix_file(modkit_calls, temp_file_none, samples = "sample1", mod_code = "c", verbose = FALSE),
        "No calls with mod_code 'c'.*found in file: 'm', 'h', 'a'"
    )
    expect_false(file.exists(temp_file_none))
})
