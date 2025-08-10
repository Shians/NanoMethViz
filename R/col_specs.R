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
        chromosome = readr::col_character(),
        start = readr::col_integer(),
        read_name = readr::col_character(),
        log_lik_ratio = readr::col_double(),
        num_cpgs = readr::col_double(),
        sequence = readr::col_character()
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

megalodon_col_types_old <- function() {
    readr::cols_only(
        read_id = readr::col_character(),
        chrm = readr::col_character(),
        strand = readr::col_integer(),
        pos = readr::col_integer(),
        mod_log_prob = readr::col_double()
    )
}

megalodon_col_types <- function() {
    readr::cols_only(
        read_id = readr::col_character(),
        chrm = readr::col_character(),
        strand = readr::col_character(),
        pos = readr::col_integer(),
        mod_log_prob = readr::col_double()
    )
}

bedmethy_col_names <- function() {
    c(
        "chrom",
        "start",
        "end",
        "modified_base_code",
        "score",
        "strand",
        "start_incl",
        "end_incl",
        "color",
        "n_valid_cov",
        "fraction_modified",
        "n_mod",
        "n_canonical",
        "n_other_mod",
        "n_delete",
        "n_fail",
        "n_diff",
        "n_nocall"
    )
}

modkit_col_types <- function() {
    readr::cols_only(
        read_id = readr::col_character(),
        forward_read_position = readr::col_double(),
        ref_position = readr::col_double(),
        chrom = readr::col_character(),
        mod_strand = readr::col_character(),
        ref_strand = readr::col_character(),
        ref_mod_strand = readr::col_character(),
        fw_soft_clipped_start = readr::col_double(),
        fw_soft_clipped_end = readr::col_double(),
        read_length = readr::col_double(),
        mod_qual = readr::col_double(),
        mod_code = readr::col_character(),
        base_qual = readr::col_double(),
        ref_kmer = readr::col_character(),
        query_kmer = readr::col_character(),
        canonical_base = readr::col_character(),
        modified_primary_base = readr::col_character(),
        inferred = readr::col_logical(),
        flag = readr::col_double()
    )
}
