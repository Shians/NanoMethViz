#' @importFrom RSQLite dbConnect SQLite SQLITE_RO dbDisconnect dbGetQuery
#' @importFrom Rsamtools TabixFile scanTabix

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

query_methy <- function(x, chr, start, end) {
    if (can_open_sql(x)) {
        query_methy_sqlite(x, chr, start, end)
    } else if (can_open_tabix(x)) {
        query_methy_tabix(x, chr, start, end)
    } else {
        stop("'x' is not a recognised file of type sqlite3 or tabix")
    }
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


#' Column names for methylation data
#'
#' @return column names for methylation data
#' @export
#'
#' @examples
#' methy_data_cols()
methy_data_cols <- function() {
    c(
        "sample",
        "chr",
        "pos",
        "strand",
        "modified",
        "statistic",
        "read_name"
    )
}

query_methy_tabix <- function(x, chr, start, end) {
    tabix_file <- Rsamtools::TabixFile(x)

    query <- GenomicRanges::GRanges(glue::glue("{chr}:{start}-{end}"))

    col_names <- methy_data_cols()

    col_types <- readr::cols(
        "sample" = readr::col_character(),
        "chr" = readr::col_character(),
        "pos" = readr::col_integer(),
        "strand" = readr::col_character(),
        "modified" = readr::col_logical(),
        "statistic" = readr::col_double(),
        "read_name" = readr::col_character()
    )

    query_result <- Rsamtools::scanTabix(tabix_file, param = query)

    parse_tabix <- function(x) {
        if (length(x) == 0) {
            return(
                tibble::tibble(
                    "sample" = character(),
                    "chr" = character(),
                    "pos" = integer(),
                    "strand" = character(),
                    "modified" = logical(),
                    "statistic" = numeric(),
                    "read_name" = character()
                )
            )
        }
        if (length(x) == 1) {
            x <- paste0(x, "\n")
        }
        readr::read_tsv(x, col_names = col_names, col_types = col_types)
    }

    purrr::map(
        query_result,
        parse_tabix
    )
}
