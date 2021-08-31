# expand multiple motifs
expand_motifs <- function(x) {
    x_mult <- x[x$num_cpgs > 1, ]
    x <- x[x$num_cpgs == 1, ]
    x$num_cpgs <- NULL

    mult_expand <- with(x_mult, {
        start_offsets <- stringr::str_locate_all(
            stringr::str_sub(sequence, start = 6),
            "CG"
        )
        start_offsets <- lapply(
            start_offsets,
            function(x) {
                x[, 1] - 1
            }
        )

        start_offsets <- unlist(start_offsets)
        tidyr::uncount(x_mult, .data$num_cpgs) %>%
            dplyr::mutate(start = .data$start + start_offsets) %>%
            dplyr::select(!"num_cpgs")
    })

    rbind(x, mult_expand) %>%
        dplyr::select(-sequence)
}

reformat_f5c <- function(x, sample) {
    x <- x %>%
        expand_motifs() %>%
        add_column(sample = sample, .before = 1)

    x %>%
        dplyr::transmute(
            sample = factor(.data$sample),
            chr = factor(.data$chromosome),
            pos = as.integer(.data$start),
            strand = factor("*", levels = c("+", "-", "*")),
            statistic = .data$log_lik_ratio,
            read_name = .data$read_name
        )
}

reformat_nanopolish <- function(x, sample) {
    x <- x %>%
        dplyr::rename(num_cpgs = "num_motifs") %>%
        expand_motifs() %>%
        add_column(sample = sample, .before = 1)

    x %>%
        dplyr::transmute(
            sample = factor(.data$sample),
            chr = factor(.data$chromosome),
            pos = as.integer(.data$start),
            strand = .data$strand,
            statistic = .data$log_lik_ratio,
            read_name = .data$read_name
        )
}

reformat_megalodon_old <- function(x, sample) {
    x %>%
        rename(
            chr = .data$chrm,
            statistic = .data$mod_log_prob,
            read_name = .data$read_id) %>%
        add_column(sample = sample, .before = 1) %>%
        mutate(
            sample = as.factor(.data$sample),
            chr = factor(.data$chr),
            statistic = logit(exp(.data$statistic)),
            strand = case_when(
                strand == 1 ~ "+",
                strand == -1 ~ "-",
                TRUE ~ "*"),
            pos = as.integer(.data$pos) + 1,
            strand = factor(.data$strand, levels = c("+", "-", "*"))) %>%
        select(methy_col_names())
}

reformat_megalodon <- function(x, sample) {
    x %>%
        rename(
            chr = .data$chrm,
            statistic = .data$mod_log_prob,
            read_name = .data$read_id) %>%
        add_column(sample = sample, .before = 1) %>%
        mutate(
            sample = as.factor(.data$sample),
            chr = factor(.data$chr),
            statistic = logit(exp(.data$statistic)),
            pos = as.integer(.data$pos) + 1,
            strand = factor(.data$strand, levels = c("+", "-", "*"))) %>%
        select(methy_col_names())
}

guess_methy_source <- function(methy_file) {
    assert_that(is.readable(methy_file))

    first_line <- readr::read_lines(methy_file, n_max = 1)

    switch (
        first_line,
        "chromosome\tstart\tend\tread_name\tlog_lik_ratio\tlog_lik_methylated\tlog_lik_unmethylated\tnum_calling_strands\tnum_cpgs\tsequence" = "f5c",
        "chromosome\tstrand\tstart\tend\tread_name\tlog_lik_ratio\tlog_lik_methylated\tlog_lik_unmethylated\tnum_calling_strands\tnum_motifs\tsequence" = "nanopolish",
        "read_id\tchrm\tstrand\tpos\tmod_log_prob\tcan_log_prob\tmod_base\tmotif" = "megalodon",
        "read_id\tchrm\tstrand\tpos\tmod_log_prob\tcan_log_prob\tmod_base" = "megalodon",
        stop("Format not recognised.")
    )
}

#' Convert methylation calls to NanoMethViz format
#' @keywords internal
#'
#' @param input_files the files to convert
#' @param output_file the output file to write results to
#' @param samples the names of samples corresponding to each file
#' @param verbose TRUE if progress messages are to be printed
#'
#' @return invisibly returns the output file path, creates a tabix file (.bgz)
#'   and its index (.bgz.tbi)
convert_methy_format <- function(
    input_files,
    output_file,
    samples = fs::path_ext_remove(fs::path_file(input_files)),
    verbose = TRUE
) {
    for (f in input_files) {
        assert_that(is.readable(f))
    }

    assert_that(
        is.character(output_file)
    )

    assert_that(is.dir(fs::path_dir(output_file)))
    file.create(path.expand(output_file))
    assert_that(is.writeable(output_file))

    for (element in vec_zip(file = input_files, sample = samples)) {
        if (verbose) {
            message(glue::glue("processing {element$file}..."))
        }
        methy_source <- guess_methy_source(element$file)
        if (verbose) {
            message(glue::glue("guessing file is produced by {methy_source}..."))
        }

        col_types <- switch (
            methy_source,
            "nanopolish" = nanopolish_col_types(),
            "f5c" = f5c_col_types(),
            "megalodon" = megalodon_col_types()
        )

        reformatter <- switch (
            methy_source,
            "nanopolish" = reformat_nanopolish,
            "f5c" = reformat_f5c,
            "megalodon" = reformat_megalodon
        )

        callback_f <- function(x, i) {
            data.table::fwrite(
                reformatter(x, sample = element$sample),
                file = output_file,
                sep = "\t",
                append = TRUE,
                scipen = 999L
            )
        }

        readr::read_tsv_chunked(
            element$file,
            col_types = col_types,
            readr::SideEffectChunkCallback$new(callback_f)
        )
    }

    invisible(output_file)
}
