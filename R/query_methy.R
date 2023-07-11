#' Query methylation data
#'
#' @param x the NanoMethResults object or a path to the methylation data
#'   (tabix-bgzipped).
#' @param chr the vector of chromosomes
#' @param start the vector of start positions
#' @param end the vector of end positions
#' @param simplify whether returned results should be row-concatenated
#' @param force whether to force empty output when query region 'chr' does not
#'   appear in data. Without 'force', an empty result indicates that the
#'   requested 'chr' appears in the data but no data overlaps with requested
#'   region, and an invalid 'chr' will cause an error.
#' @param truncate when querying from ModBamFiles, whether or not to truncate
#'   returned results to only those within the specified region. Otherwise
#'   methylation data for entire reads overlapping the region will be returned.
#'
#' @return A table containing the data within the queried regions. If simplify
#'   is TRUE (default) then all data is contained within one table, otherwise it
#'   is a list of tables where each element is the data from one region.
#'
#'
#' @examples
#' nmr <- load_example_nanomethresult()
#' query_methy(methy(nmr), "chr7", 6703892, 6730431)
#'
#' @importFrom Rsamtools TabixFile scanTabix
#'
#' @export
query_methy <- function(x, chr, start, end, simplify = TRUE, force = FALSE, truncate = TRUE) {
    if (is(x, "NanoMethResult")) {
        x <- methy(x)
    }

    if (is(x, "ModBamResult")) {
        mod_code <- mod_code(x)
    }

    assert_that(
        same_length(chr, start, end),
        msg = "vectors 'chr', 'start' and 'end' must be the same length"
    )

    if (is(x, "ModBamResult")) {
        out <- query_methy_modbam(x, chr, start, end, mod_code)
    } else if (can_open_tabix(x)) {
        out <- query_methy_tabix(x, chr, start, end, force = force)
    } else {
        stop("'x' is not a recognised file type")
    }

    out <- map(out, function(x) {
        dplyr::mutate(x, mod_prob = sigmoid(.data$statistic))
    })

    if (truncate) {
        pos_list <- purrr::map2(start, end, ~c(start = .x, end = .y))

        truncate_fn <- function(x, pos_range) {
            x %>%
                dplyr::filter(
                    .data$pos >= pos_range["start"],
                    .data$pos <= pos_range["end"]
                )
        }
        out <- purrr::map2(out, pos_list, truncate_fn)
    }

    if (simplify) {
        out <- dplyr::bind_rows(out)
    }

    out
}

query_methy_df <- function(x, regions, simplify = TRUE, force = FALSE) {
    assert_has_columns(regions, c("chr", "start", "end"))
    query_methy(x, regions$chr, regions$start, regions$end, simplify, force)
}

query_methy_gr <- function(x, regions, simplify = TRUE, force = FALSE) {
    assert_that(is(regions, "GRanges"))
    query_methy(
        x,
        as.character(GenomicRanges::seqnames(regions)),
        as.numeric(GenomicRanges::start(regions)),
        as.numeric(GenomicRanges::end(regions)),
        simplify,
        force
    )
}

query_methy_gene <- function(x, gene, window_prop = 0, simplify = TRUE) {
    if (!gene %in% exons(x)$symbol) {
        stop(glue::glue("gene '{gene}' not found in NanoMethViz::exons(x)"))
    }

    if (length(window_prop) == 1) {
        # convert to two sided window
        window_prop <- c(window_prop, window_prop)
    }

    pos_range <- gene_pos_range(x, gene)

    gene_width <- pos_range[2] - pos_range[1]
    window_left <- gene_width * window_prop[1]
    window_right <- gene_width * window_prop[2]

    chr <- exons(x) %>%
        dplyr::filter(.data$symbol == gene) %>%
        dplyr::slice(1) %>%
        dplyr::pull(chr)

    query_methy(
        x,
        chr = chr,
        start = pos_range[1] - window_left,
        end = pos_range[2] + window_right,
        simplify = simplify
    )
}

can_open_tabix <- function(x) {
    assert_readable(x)
    out <- TRUE

    tryCatch(
        Rsamtools::TabixFile(x),
        warning = function(x) { out <<- FALSE },
        error = function(x) { out <<- FALSE }
    )

    return(out)
}

empty_methy_query_output <- function() {
    tibble::tibble(
        "sample" = character(),
        "chr" = character(),
        "pos" = integer(),
        "strand" = character(),
        "statistic" = numeric(),
        "read_name" = character()
    )
}

#' @importFrom utils read.table
query_methy_tabix <- function(x, chr, start, end, force) {
    tabix_file <- Rsamtools::TabixFile(x)

    tabix_seqs <- get_tabix_sequences(paste0(x, ".tbi"))
    miss <- which(!chr %in% tabix_seqs)
    miss_seqs <- unique(chr[miss])
    if (length(miss_seqs) != 0) {
        warning(
            "requested sequences missing from tabix file:",
            paste(miss_seqs, collapse = ", ")
        )
        chr <- chr[-miss]
        start <- start[-miss]
        end <- end[-miss]
    }

    if (length(chr) == 0) {
        if (!force) {
            stop("no chromosome matches between query and tabix file, please check chromosome format matches between query and methylation file.")
        } else {
            return(empty_methy_query_output())
        }
    }

    query <- make_granges(chr, start, end)

    col_names <- methy_col_names()
    col_types <- methy_col_types()

    query_result <- Rsamtools::scanTabix(tabix_file, param = query)

    parse_tabix <- function(x) {
        if (length(x) == 0) {
            return(empty_methy_query_output())
        }
        if (length(x) == 1) {
            x <- paste0(x, "\n")
        }

        # using readr::read_tsv on character vectors seems to leak memory
        as_tibble(
            utils::read.table(
                textConnection(x),
                col.names = col_names,
                sep = "\t",
                colClasses = c(
                    "character",
                    "character",
                    "integer",
                    "character",
                    "numeric",
                    "character"
                ),
                header = FALSE
            )
        )
    }

    lapply(
        query_result,
        parse_tabix
    )
}

query_methy_modbam <- function(x, chr, start, end, mod_code) {
    assertthat::assert_that(
        is(x, "ModBamResult") ||
        is(x, "ModBamFiles")
    )
    if (is(x, "ModBamResult")) {
        x <- methy(x)
    }

    assert_readable(x$path)

    # query each file
    x <- data.frame(
        sample = x$sample,
        path = x$path
    )
    out <- x %>%
        dplyr::mutate(
            mod_table = map_rows(x, function(x) {
                read_modbam_table(
                    x$path,
                    chr = chr,
                    start = start,
                    end = end,
                    sample = x$sample,
                    mod_code = mod_code)
            })
        )

    # reduce list nesting by one level
    tables <- do.call(c, out$mod_table)

    # assign bind together tables from the same regions
    nms <- names(tables)

    if (is.null(nms)) {
        warning(glue::glue("no data found in {chr}:{start}-{end}"))
    }
    out <- list()
    for (nm in unique(nms)) {
        out[[nm]] <- do.call(rbind, tables[names(tables) == nm])
    }

    out
}

