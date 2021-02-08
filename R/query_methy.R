#' Query methylation data
#'
#' @param x the path to the methylation data (tabix-bgzipped)
#' @param chr the vector of chromosomes
#' @param start the vector of start positions
#' @param end the vector of end positions
#' @param simplify whether returned results should be row-concatenated
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
#' @importFrom RSQLite dbConnect SQLite SQLITE_RO dbDisconnect dbGetQuery
#' @importFrom Rsamtools TabixFile scanTabix
#'
#' @export
query_methy <- function(x, chr, start, end, simplify = TRUE) {
    if (is(x, "NanoMethResult")) {
        x <- methy(x)
    }

    if (can_open_sql(x)) {
        out <- query_methy_sqlite(x, chr, start, end)
    } else if (can_open_tabix(x)) {
        out <- query_methy_tabix(x, chr, start, end)
    } else {
        stop("'x' is not a recognised file of type sqlite3 or tabix")
    }

    if (simplify) {
        out <- dplyr::bind_rows(out)
    }

    out
}

query_methy_gene <- function(x, gene, simplify = TRUE) {
    if (!gene %in% exons(x)$symbol) {
        stop(glue::glue("gene {gene} not found in exon annotation"))
    }

    pos_range <- gene_pos_range(x, gene)
    chr <- exons(x) %>%
        dplyr::filter(symbol == gene) %>%
        dplyr::slice(1) %>%
        dplyr::pull(chr)

    query_methy(
        x,
        chr = chr,
        start = pos_range[1],
        end = pos_range[2],
        simplify = simplify
    )
}

can_open_sql <- function(x) {
    assertthat::is.readable(x)
    out <- TRUE

    tryCatch(
        RSQLite::dbConnect(
            x,
            drv = RSQLite::SQLite(),
            flags = RSQLite::SQLITE_RO
        ),
        warning = function(x) { out <<- FALSE },
        error = function(x) { out <<- FALSE }
    )

    return(out)
}

can_open_tabix <- function(x) {
    assertthat::is.readable(x)
    out <- TRUE

    tryCatch(
        Rsamtools::TabixFile(x),
        warning = function(x) { out <<- FALSE },
        error = function(x) { out <<- FALSE }
    )

    return(out)
}

query_methy_sqlite <- function(x, chr, start, end) {
    db <- RSQLite::dbConnect(
        x,
        drv = RSQLite::SQLite(),
        flags = RSQLite::SQLITE_RO
    )

    query <- glue::glue("SELECT * FROM methylation
                        WHERE chr = '{chr}'
                        AND pos BETWEEN {start} AND {end}")

    out <- RSQLite::dbGetQuery(db, query)

    RSQLite::dbDisconnect(db)

    tibble::as_tibble(out)
}

#' @importFrom utils read.table
query_methy_tabix <- function(x, chr, start, end) {
    tabix_file <- Rsamtools::TabixFile(x)

    query <- GenomicRanges::GRanges(glue::glue("{chr}:{start}-{end}"))

    col_names <- methy_col_names()
    col_types <- methy_col_types()

    query_result <- Rsamtools::scanTabix(tabix_file, param = query)

    parse_tabix <- function(x) {
        if (length(x) == 0) {
            return(
                tibble::tibble(
                    "sample" = character(),
                    "chr" = character(),
                    "pos" = integer(),
                    "strand" = character(),
                    "statistic" = numeric(),
                    "read_name" = character()
                )
            )
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
