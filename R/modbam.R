# get modification probability values and reference-based coordinates
modbam_to_ref_coord <- function(seq, cigar, mod_str, mod_scores, map_pos, strand) {
    # Check inputs
    assert_that(assertthat::is.string(seq))
    assert_that(assertthat::is.string(cigar))
    assert_that(assertthat::is.string(mod_str))
    assert_that(assertthat::is.string(mod_scores))
    assert_that(is.numeric(map_pos))
    assert_that(strand %in% c("+", "-"), msg = "Strand must be '+' or '-'")

    # tokenise the cigar string into a data.frame
    mod_tokens <- mod_tokeniser_cpp(mod_str, mod_scores)

    # get the mapping from strand coordinate to relative genomic coordinate
    coord_map <- get_coord_map_cpp(cigar)

    # get positions of base of interest
    if (strand == "-") {
        mod_candidate_pos <- rev(get_char_pos_cpp(seq, "G"))[mod_tokens$mod_pos]
    } else if (strand == "+") {
        mod_candidate_pos <- get_char_pos_cpp(seq, "C")[mod_tokens$mod_pos]
    }

    # calculate the relative position
    rel_pos <- coord_map[mod_candidate_pos]
    names(rel_pos) <- NULL

    # add the mapping position offset to the relative position
    if (strand == "-") {
        genome_pos <- map_pos + rel_pos - 2
    } else {
        genome_pos <- map_pos + rel_pos - 1
    }

    data.frame(
        pos = genome_pos,
        statistic = logit(mod_tokens$mod_prob)
    ) %>%
        tidyr::drop_na()
}

parse_bam_list <- function(seq, cigar, mod_str, mod_scores, map_pos, strand, mod_code) {
    # Check inputs
    assert_that(is.character(seq))
    assert_that(is.character(cigar))
    assert_that(is.character(mod_str))
    assert_that(is.character(mod_scores))
    assert_that(is.numeric(map_pos))
    assert_that(all(strand %in% c("+", "-")), msg = "Strand must be '+' or '-'")
    assert_that(is.character(mod_code) && nchar(mod_code) == 1)

    parse_bam_list_cpp(seq, cigar, mod_str, mod_scores, map_pos, as.character(strand), mod_code)
}

parse_modbam <- function(x, sample, mod_code) {
    if (is.null(x$tag$ML)) {
        mod_scores <- NULL
    } else {
        mod_scores <- x$tag$ML %>%
            purrr::map_chr(~stringr::str_c(., collapse = ","))
    }

    reads_df <- tibble::tibble(
        chr = x$rname,
        strand = x$strand,
        map_pos = x$pos,
        cigar = x$cigar,
        seq = as.character(x$seq),
        mod_string = x$tag$MM,
        mod_scores = mod_scores,
        read_name = x$qname
    ) %>%
        dplyr::filter(seq != "")

    if (nrow(reads_df) == 0) {
        return(NULL)
    }

    out <- reads_df %>%
        mutate(
            modbam_stats = parse_bam_list(
                .data$seq,
                .data$cigar,
                .data$mod_string,
                .data$mod_scores,
                .data$map_pos,
                as.character(.data$strand),
                mod_code
            )

        )

    out %>%
        dplyr::select("read_name", "chr", "strand", "modbam_stats") %>%
        tidyr::unnest("modbam_stats") %>%
        tidyr::drop_na() %>%
        dplyr::mutate(sample = sample, .before = 1)
}

read_bam <- function(bam_file, query = NULL) {
    modbam_param <- function(query = NULL) {
        if (!is.null(query)) {
            Rsamtools::ScanBamParam(
                flag = Rsamtools::scanBamFlag(isUnmappedQuery = FALSE),
                what = c("qname", "rname", "strand", "pos", "cigar", "seq"),
                tag = c("MM", "ML"),
                which = query
            )
        } else {
            Rsamtools::ScanBamParam(
                flag = Rsamtools::scanBamFlag(isUnmappedQuery = FALSE),
                what = c("qname", "rname", "strand", "pos", "cigar", "seq"),
                tag = c("MM", "ML")
            )
        }
    }

    # determine_tag <- function(records) {
    #     if (is.null(records$tag$Mm) && is.null(records$tag$Ml)) {
    #         records$tag$Mm <- NULL
    #         records$tag$Ml <- NULL
    #     } else (is.null(records$tag$MM) && is.null(records$tag$ML)) {
    #         records$tag$MM <- records$tag$Mm
    #         records$tag$ML <- records$tag$Ml
    #         records$tag$Mm <- NULL
    #         records$tag$Ml <- NULL
    #     }
    #
    #     records
    # }

    filter_modbam <- function(x) {
        tag <- x$tag
        x$tag <- NULL

        missing <- map_lgl(tag$ML, ~length(.) == 0) |
            map_lgl(x$rname, is.na) |
            map_lgl(x$strand, is.na) |
            map_lgl(x$pos, is.na) |
            map_lgl(x$cigar, is.na)

        x <- map(x, ~.[!missing])
        tag <- map(tag, ~.[!missing])

        x$tag <- tag
        x
    }

    Rsamtools::scanBam(
        bam_file,
        param = modbam_param(query = query)
    ) %>%
        map(filter_modbam)
}

read_modbam_table <- function(x, chr, start, end, sample, mod_code) {
    bam_file <- Rsamtools::BamFile(x)

    query <- make_granges(chr, start, end)

    if (length(query) == 0) {
        return(list(empty_methy_query_output()))
    }

    reads <- read_bam(bam_file, query = query)

    lapply(reads, parse_modbam, sample = sample, mod_code = mod_code)
}
