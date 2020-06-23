#' Column names for methylation data
#'
#' @return column names for methylation data
#' @export
#'
#' @examples
#' methy_col_names()
methy_col_names <- function() {
    c(
        "sample",
        "chr",
        "pos",
        "strand",
        "statistic",
        "read_name"
    )
}

methy_col_types <- function() {
    readr::cols_only(
        sample = readr::col_factor(),
        chr = readr::col_factor(),
        pos = readr::col_integer(),
        strand = readr::col_factor(levels = c("+", "-", "*")),
        statistic = readr::col_double(),
        read_name = readr::col_character()
    )
}

f5c_col_types <- function() {
    readr::cols_only(
        chromosome = col_character(),
        start = col_integer(),
        read_name = col_character(),
        log_lik_ratio = col_double(),
        num_cpgs = col_double(),
        sequence = col_character()
    )
}

nanopolish_col_types <- function() {
    readr::cols_only(
        chromosome = readr::col_character(),
        start = readr::col_integer(),
        strand = readr::col_factor(levels = c("+", "-", "*")),
        read_name = readr::col_character(),
        log_lik_ratio = readr::col_double(),
        num_motifs = readr::col_double(),
        sequence = readr::col_character()
    )
}

megalodon_col_types <- function() {
    readr::cols_only(
        read_id = readr::col_character(),
        chrm = readr::col_character(),
        strand = readr::col_integer(),
        pos = readr::col_integer(),
        mod_log_prob = readr::col_double()
    )
}
