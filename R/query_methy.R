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
#' @param site_filter the minimum amount of coverage to report a site. This
#'   filters the queried data such that any site with less than the filter is
#'   not returned. The default is 1, which means that all sites are returned.
#'   This option can be set globally using the `options(site_filter = ...)`
#'   which will affect all plotting functions in NanoMethviz.
#'
#' @return A table containing the data within the queried regions. If simplify
#'   is TRUE (default) then returns all data in a single table, otherwise returns
#'   a list of tables where each table is the data from one region.
#'
#' @details
#' The argument `site_filter` can be set globally using the `options(site_filter
#' = ...)` command. The same data entry may appear multiple times in the output
#' if it overlaps multiple regions.
#'
#' @examples
#' nmr <- load_example_nanomethresult()
#' query_methy(methy(nmr), "chr7", 6703892, 6730431)
#'
#' @importFrom Rsamtools TabixFile scanTabix
#'
#' @export
query_methy <- function(
    x,
    chr,
    start,
    end,
    simplify = TRUE,
    force = FALSE,
    truncate = TRUE,
    site_filter = getOption("NanoMethViz.site_filter", 3L)
) {
    # validate input
    chr <- as.character(chr)
    start <- as.numeric(start)
    end <- as.numeric(end)

    # validate genomic coordinates
    assert_valid_genomic_coords(chr, start, end)

    # validate site_filter
    if (!is.numeric(site_filter) || site_filter < 0) {
        stop(glue::glue(
            "site_filter must be a non-negative number. Got: {site_filter}\n",
            "This parameter filters sites with coverage below the threshold. Set using `options(site_filter = ...)`."
        ))
    }

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

    assert_that(site_filter > 0)

    input_type <- guess_input_type(x)
    if (input_type == "tabix") {
        out <- query_methy_tabix(x, chr, start, end, force = force)

    } else if (input_type == "modbam") {
        out <- query_methy_modbam(x, chr, start, end, mod_code)
        if (truncate) {
            truncate_fn <- function(x, start, end) {
                x %>%
                    dplyr::filter(
                        .data$pos >= start,
                        .data$pos <= end
                    )
            }

            # truncate each query by its start and end positions
            out <- purrr::pmap(list(out, start, end), truncate_fn)
        }
    } else {
        stop("'x' is not a recognised file type")
    }

    # compute modification probability column
    out <- map(out, function(x) {
        x %>%
            dplyr::mutate(mod_prob = sigmoid(.data$statistic))
    })

    # filter output for coverage
    out <- purrr::map(out, function(x) {
        pos_count_filter <- x %>%
            dplyr::count(.data$chr, .data$pos) %>%
            dplyr::rename("coverage" = "n") %>%
            dplyr::filter(.data$coverage >= site_filter) %>%
            dplyr::select(-"coverage")

        x %>%
            dplyr::inner_join(
                pos_count_filter,
                by = c("chr", "pos")
            )
    })

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

#' @importFrom utils read.table
query_methy_tabix <- function(x, chr, start, end, force) {
    get_missing_seqs <- function(x) {
        # query tabix index for missing sequences
        tabix_seqs <- get_tabix_sequences(paste0(x, ".tbi"))

        which(!chr %in% tabix_seqs)
    }

    # set up tabix file
    tabix_file <- Rsamtools::TabixFile(x)

    # set up output
    out <- list()
    for (i in seq_along(chr)) {
        nm <- paste0(chr[i], ":", start[i], "-", end[i])
        out[[nm]] <- empty_methy_query_output()
    }

    miss <- get_missing_seqs(x)
    # if missing sequences, warn and remove from query
    if (length(miss) != 0) {
        miss_seqs <- unique(chr[miss])
        warning(
            "requested sequences missing from tabix file and will be excluded from query:",
            paste(miss_seqs, collapse = ", ")
        )

        # remove queries with missing sequences
        chr <- chr[-miss]
        start <- start[-miss]
        end <- end[-miss]
    }

    # if no sequences left, return empty output
    if (length(chr) == 0) {
        if (!force) {
            stop("no chromosome matches between query and tabix file, please check chromosome format matches between query and methylation file.")
        }
    }

    # make query into a granges object
    query <- make_granges(chr, start, end)
    if (length(query) > 0) {
        query_result <- Rsamtools::scanTabix(tabix_file, param = query)
        # helper function to parse tabix output
        parse_tabix <- function(x) {
            if (length(x) == 0) {
                return(empty_methy_query_output())
            }

            # help readr function recognise input as data text rather than path
            if (length(x) == 1) {
                x <- paste0(x, "\n")
            }

            read_methy_lines(x)
        }

        methy_data <- lapply(
            query_result,
            parse_tabix
        )

        for (i in seq_along(methy_data)) {
            out[[names(query_result)[i]]] <- methy_data[[i]]
        }
    }

    out
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

    # get data for each region from each file
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

    # extract and bind data from each region from each file
    extract_and_bind_data <- function(nm, tables) {
        out <- do.call(rbind, tables[names(tables) == nm])
        if (is.null(out)) {
            out <- empty_methy_query_output()
        }
        out
    }

    out <- purrr::map(
        unique(nms),
        extract_and_bind_data,
        tables = tables
    )

    out
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

guess_input_type <- function(x) {
    if (is(x, "ModBamResult")) {
        return("modbam")
    }

    if (is(x, "character")) {
        if (all(map_lgl(x, can_open_tabix))) {
            return("tabix")
        }
    }

    return("uknown")
}

read_methy_lines <- function(x) {
    read.delim(
        textConnection(x),
        header = FALSE,
        col.names = methy_col_names()
    ) %>%
        tibble::tibble()
}
