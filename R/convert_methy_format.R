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
        dplyr::rename(
            chr = "chrm",
            statistic = "mod_log_prob",
            read_name = "read_id") %>%
        add_column(sample = sample, .before = 1) %>%
        dplyr::mutate(
            sample = as.factor(.data$sample),
            chr = factor(.data$chr),
            statistic = logit(exp(.data$statistic)),
            pos = as.integer(.data$pos) + 1,
            strand = factor(.data$strand, levels = c("+", "-", "*"))) %>%
        dplyr::select(methy_col_names())
}

# TRUE where the read base after the centre of the read-oriented kmer is G,
# i.e. the called base is the C of a CG on the read
is_read_cg <- function(canonical_base, query_kmer) {
    centre <- (nchar(query_kmer) + 1) %/% 2
    next_base <- toupper(substr(query_kmer, centre + 1, centre + 1))
    !is.na(next_base) & canonical_base == "C" & next_base == "G"
}

reformat_modkit <- function(x, sample, mod_code = "m") {
    x %>%
        dplyr::filter(ref_position >= 0) %>% # remove unmapped positions
        dplyr::filter(.data$mod_code == !!mod_code) %>%
        add_column(sample = sample, .before = 1) %>%
        dplyr::rename(
            chr = "chrom",
            pos = "ref_position",
            # mod_strand is relative to the read, ref_mod_strand to the
            # reference
            strand = "ref_mod_strand",
            statistic = "mod_qual",
            read_name = "read_id"
        ) %>%
        dplyr::mutate(
            sample = as.factor(.data$sample),
            chr = factor(.data$chr),
            # reverse strand CG calls sit on the G, shift them onto the C of
            # the forward strand so both strands share one coordinate
            shift = .data$strand == "-" &
                is_read_cg(.data$canonical_base, .data$query_kmer),
            pos = as.integer(.data$pos) + 1L - .data$shift,
            strand = factor(.data$strand, levels = c("+", "-", "*")),
            statistic = logit(.data$statistic)
        ) %>%
        dplyr::select(
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

# Validate mod_code against the input sources and fill in the modkit default.
resolve_mod_code <- function(mod_code, methy_sources) {
    has_modkit <- any(methy_sources == "modkit")

    if (!is.null(mod_code) && !has_modkit) {
        stop(glue::glue(
            "mod_code only applies to modkit input, but no input file was ",
            "detected as modkit (detected: {paste(unique(methy_sources), collapse = ', ')}).\n",
            "Remove the mod_code argument."
        ))
    }

    if (is.null(mod_code)) {
        return(if (has_modkit) "m" else NULL)
    }

    if (!is.string(mod_code) || is.na(mod_code) || mod_code == "") {
        stop("mod_code must be a single non-empty string, e.g. \"m\" for 5mC.")
    }

    mod_code
}

#' Convert methylation calls to NanoMethViz format
#' @keywords internal
#'
#' @param input_files the files to convert
#' @param output_file the output file to write results to (must end in .bgz)
#' @param samples the names of samples corresponding to each file
#' @param mod_code the modification code to extract from modkit input. NULL
#'   uses "m" (5mC). Must be NULL unless at least one input is from modkit.
#' @param verbose TRUE if progress messages are to be printed
#'
#' @return invisibly returns the output file path, creates a tabix file (.bgz)
#'   and its index (.bgz.tbi)
convert_methy_format <- function(
    input_files,
    output_file,
    samples = extract_file_names(input_files),
    mod_code = NULL,
    verbose = TRUE
) {
    for (f in input_files) {
        assert_readable(f)
    }

    methy_sources <- purrr::map_chr(input_files, guess_methy_source)
    mod_code <- resolve_mod_code(mod_code, methy_sources)

    assert_that(
        is.character(output_file)
    )

    assert_that(is.dir(fs::path_dir(output_file)))
    file.create(path.expand(output_file))
    assert_that(is.writeable(output_file))

    for (element in vec_zip(file = input_files, sample = samples, source = methy_sources)) {
        if (verbose) {
            message(glue::glue("processing {element$file}..."))
        }
        methy_source <- element$source
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
            "modkit" = function(x, sample) reformat_modkit(x, sample, mod_code = mod_code)
        )

        # track modkit codes seen so an unmatched mod_code fails loudly
        # rather than producing an empty file
        rows_written <- 0
        codes_seen <- character()
        writer_fn <- function(x, i) {
            if (methy_source == "modkit") {
                codes_seen <<- union(codes_seen, unique(x$mod_code))
            }
            out <- reformatter(x, sample = element$sample)
            rows_written <<- rows_written + nrow(out)
            readr::write_tsv(out, file = output_file, append = TRUE)
        }
        readr::local_edition(1) # temporary fix for vroom bad value
        readr::read_tsv_chunked(
            element$file,
            col_types = col_types,
            readr::SideEffectChunkCallback$new(writer_fn)
        )

        if (methy_source == "modkit" && rows_written == 0) {
            codes_seen <- codes_seen[!is.na(codes_seen)]
            found <- if (length(codes_seen) > 0) {
                paste0("'", codes_seen, "'", collapse = ", ")
            } else {
                "none"
            }
            stop(glue::glue(
                "No calls with mod_code '{mod_code}' found in modkit file '{element$file}'.\n",
                "Mod codes found in file: {found}.\n",
                "Set mod_code to one of the codes present."
            ))
        }
    }

    invisible(output_file)
}
