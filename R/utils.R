#' Load an example NanoMethResult object
#'
#' Load an example NanoMethResult object for demonstration of plotting
#' functions. Run `load_example_nanomethresults` without the function call to
#' see how the object is constructed.
#'
#' @return a NanoMethResults object
#'
#' @export
#'
#' @examples
#' nmr <- load_example_nanomethresult()
load_example_nanomethresult <- function() {
    methy <- system.file(package = "NanoMethViz", "methy_subset.tsv.bgz")

    sample <- c(
        "B6Cast_Prom_1_bl6",
        "B6Cast_Prom_1_cast",
        "B6Cast_Prom_2_bl6",
        "B6Cast_Prom_2_cast",
        "B6Cast_Prom_3_bl6",
        "B6Cast_Prom_3_cast"
    )
    group <- c(
        "bl6",
        "cast",
        "bl6",
        "cast",
        "bl6",
        "cast"
    )
    sample_anno <- data.frame(sample, group, stringsAsFactors = FALSE)

    exon_tibble <- get_example_exons_mus_musculus()

    NanoMethResult(methy, sample_anno, exon_tibble)
}

#' Load an example ModBamResult object
#'
#' Load an example ModBamResult object for demonstration of plotting
#' functions. Run `load_example_modbamresult` without the function call to
#' see how the object is constructed.
#'
#' @return a ModBamResult object
#'
#' @export
#'
#' @examples
#' mbr <- load_example_modbamresult()
load_example_modbamresult <- function() {
    ModBamResult(
        methy = ModBamFiles(
            paths = system.file(package = "NanoMethViz", "peg3.bam"),
            samples = "sample1"
        ),
        samples = tibble::tibble(
            sample = "sample1",
            group = "group1"
        ),
        exons = get_example_exons_mus_musculus()
    )
}

.get_ggplot_range_x <- function(x) {
    # get x-axis range from a ggplot object
    # returns c(x_min, x_max)
    ggplot2::ggplot_build(x)$layout$panel_scales_x[[1]]$range$range
}

# create a list where the nth element contains the nth values of the original
# vectors
#' @importFrom stats setNames
vec_zip <- function(..., .names = NULL) {
    x <- do.call(data.frame, list(..., stringsAsFactors = FALSE))
    stats::setNames(split(x, seq_len(nrow(x))), .names)
}

extract_file_names <- function(x) {
    fs::path_ext_remove(fs::path_file(x))
}

logit <- function(p) {
    log(p / (1-p))
}

assert_has_columns <- function(x, cols) {
    if (!all(cols %in% colnames(x))) {
        stop(glue::glue(
            "columns missing from {input}: {missing_cols}",
            input = deparse(substitute(x)),
            missing_cols = paste(setdiff(cols, colnames(x)), collapse = ", ")
        ))
    }
}

timed_log <- function(...) {
    time_stamp <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ")
    message(time_stamp, ...)
}

gene_pos_range <- function(nmr, gene) {
    NanoMethViz::exons(nmr) %>%
        dplyr::filter(.data$symbol == !!gene) %>%
        dplyr::select("start", "end") %>%
        unlist() %>%
        range()
}

get_tabix_sequences <- function(file) {
    f <- gzfile(file, "rb")

    # read magic
    readChar(f, 4)

    # get n_ref and skip other fields
    n_ref <- readBin(f, "integer")
    for (i in seq_len(7)) readBin(f, "integer")

    # read sequences
    seqs <- character(n_ref)
    for (i in seq_len(n_ref)) {
        seqs[i] <- readBin(f, "character")
    }

    close(f)
    seqs
}

stack_plots <- function(spaghetti, heat) {
    if (is(spaghetti, "patchwork")) {
        spaghetti[[1]] + heat + spaghetti[[2]] +
            patchwork::plot_layout(
                nrow = 3,
                heights = c(spaghetti$patches$layout$heights[1], 1, spaghetti$patches$layout$heights[2])
            )
    } else {
        spaghetti + heat + patchwork::plot_layout(nrow = 2)
    }
}

map_rows <- function(x, .f, ...) {
    if (!is.data.frame(x) || nrow(x) == 0) {
        stop("'x' must be a non-empty data.frame")
    }

    # Apply the function .f to each row of x using lapply() and store the results in a list.
    out <- lapply(seq_len(nrow(x)), function(i) .f(x[i, ], ...))

    out
}

neg_strand_offset <- function(x, offset) {
    assert_has_columns(x, c("strand", "pos"))

    x %>%
        dplyr::mutate(pos = case_when(
            strand == "-" ~ pos + offset,
            TRUE ~ pos
        ))
}

missingness <- function(x) {
    mean(is.na(x))
}

mat_row_map <- function(x, f) {
    apply(x, 1, f)
}

mat_col_map <- function(x, f) {
    apply(x, 2, f)
}

df_to_matrix <- function(x, rownames = 1) {
    if (!is.null(rownames)) {
        rnames <- x[, rownames, drop = TRUE]
        x <- x[, -rownames]
    }

    x <- as.matrix(x)
    rownames(x) <- rnames
    x
}

is_coordinate <- function(x) {
    is.numeric(x) && x == round(x) && x >= 0
}

assertthat::on_failure(is_coordinate) <- function(call, env) {
  glue::glue("'{deparse(call$x)}' ({eval(call$x)}) is not a valid coordinate")
}

same_length <- function(...) {
    nulls <- purrr::map_lgl(list(...), is.null)
    lengths <- purrr::map_dbl(list(...), length)
    !any(nulls) && all(map_lgl(lengths, ~ . == lengths[1]))
}

#' @importFrom IRanges IRanges
make_granges <- function(chr, start, end) {
    assertthat::assert_that(same_length(chr, start, end))
    assertthat::assert_that(all(start < end))
    GenomicRanges::GRanges(
        seqnames = chr,
        ranges = IRanges::IRanges(start, end)
    )
}

get_bam_total_reads <- function(path) {
    sum(Rsamtools::idxstatsBam(path)[, c("mapped", "unmapped")])
}
