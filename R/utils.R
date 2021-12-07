#' Load an example NanoMethResult object
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
    spaghetti[[1]] + heat + spaghetti[[2]] +
        patchwork::plot_layout(
            nrow = 3,
            heights = c(spaghetti$patches$layout$heights[1], 1, spaghetti$patches$layout$heights[2])
        )
}
