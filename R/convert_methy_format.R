# expand multiple motifs
expand_motifs <- function(x) {
    x_mult <- x[x$num_cpgs > 1, ]
    x <- x[x$num_cpgs == 1, ]
    x$num_cpgs <- NULL

    mult_expand <- with(x_mult, {
        start_offsets <- stringr::str_locate_all(
            stringr::str_sub(sequence, start = 6),
            "CG"
        ) %>%
        purrr::map(
            function(x) {
                x[, 1] - 1
            }
        ) %>%
        unlist()

        tidyr::uncount(x_mult, .data$num_cpgs) %>%
            dplyr::mutate(start = start_offsets + .data$start) %>%
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
            pos = as.integer(.data$start) + 1,
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
            pos = as.integer(.data$start) + 1,
            strand = .data$strand,
            statistic = .data$log_lik_ratio,
            read_name = .data$read_name
        )
}

reformat_megalodon <- function(x, sample) {
    x %>%
        rename(
            chr = "chrm",
            statistic = "mod_log_prob",
            read_name = "read_id") %>%
        add_column(sample = sample, .before = 1) %>%
        mutate(
            sample = as.factor(.data$sample),
            chr = factor(.data$chr),
            statistic = logit(exp(.data$statistic)),
            pos = as.integer(.data$pos) + 1,
            strand = factor(.data$strand, levels = c("+", "-", "*"))) %>%
        select(methy_col_names())
}

reformat_modkit <- function(x, sample) {
    x %>%
        filter(ref_position > 0) %>% # remove unmapped positions
        add_column(sample = sample, .before = 1) %>%
        rename(
            chr = "chrom",
            pos = "ref_position",
            strand = "mod_strand",
            statistic = "mod_qual",
            read_name = "read_id"
        ) %>%
        mutate(
            sample = as.factor(.data$sample),
            chr = factor(.data$chr),
            pos = as.integer(.data$pos),
            strand = factor(.data$strand, levels = c("+", "-", "*")),
            statistic = logit(.data$statistic)
        ) %>%
        select(
            "sample",
            "chr",
            "pos",
            "strand",
            "statistic",
            "read_name"
        )
}

guess_methy_source <- function(methy_file) {
    assert_readable(methy_file)

    readr::local_edition(1) # temporary fix for vroom bad value
    first_line <- readr::read_lines(methy_file, n_max = 1)

    # Be robust to BOM and whitespace
    header <- stringr::str_replace(first_line, "^\ufeff", "")
    header <- stringr::str_trim(header)

    # Detect delimiter (prefer tab, fallback to comma, then whitespace)
    delim <- if (stringr::str_detect(header, "\t")) "\t"
        else if (stringr::str_detect(header, ",")) ","
        else "\\s+"

    cols <- stringr::str_split(header, delim)[[1]]
    cols <- stringr::str_trim(cols)
    cols_lower <- stringr::str_to_lower(cols)

    has_all <- function(required) all(required %in% cols_lower)

    # Identify by minimal, distinctive column sets; allow extras and any order
    if (has_all(c("read_id", "chrm", "strand", "pos", "mod_log_prob"))) {
        return("megalodon")
    }

    if (has_all(c("read_id", "ref_position", "chrom", "mod_strand", "mod_qual"))) {
        return("modkit")
    }

    if (has_all(c("chromosome", "strand", "start", "read_name", "log_lik_ratio", "num_motifs", "sequence"))) {
        return("nanopolish")
    }

    if (has_all(c("chromosome", "start", "read_name", "log_lik_ratio", "num_cpgs", "sequence"))) {
        return("f5c")
    }

    stop("Format not recognised.")
}

#' Convert methylation calls to NanoMethViz format
#' @keywords internal
#'
#' @param input_files the files to convert
#' @param output_file the output file to write results to (must end in .bgz)
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
        assert_readable(f)
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

        col_types <- switch(
            methy_source,
            "nanopolish" = nanopolish_col_types(),
            "f5c" = f5c_col_types(),
            "megalodon" = megalodon_col_types(),
            "modkit" = modkit_col_types()
        )

        reformatter <- switch(
            methy_source,
            "nanopolish" = reformat_nanopolish,
            "f5c" = reformat_f5c,
            "megalodon" = reformat_megalodon,
            "modkit" = reformat_modkit
        )

        writer_fn <- function(x, i) {
            readr::write_tsv(
                reformatter(x, sample = element$sample),
                file = output_file,
                append = TRUE
            )
        }
        readr::local_edition(1) # temporary fix for vroom bad value
        readr::read_tsv_chunked(
            element$file,
            col_types = col_types,
            readr::SideEffectChunkCallback$new(writer_fn)
        )
    }

    invisible(output_file)
}
