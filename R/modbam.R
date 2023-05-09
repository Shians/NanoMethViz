rev_df <- function(x) {
    assertthat::assert_that(is(x, "data.frame"))

    x[rev(1:nrow(x)), ]
}

mod_tokeniser <- function(string, scores) {
    string <- stringr::str_remove_all(string, ";")

    str_split_comma <- function(x, ...) {
        stringr::str_split(x, ",")[[1]]
    }

    mod_string_split <- str_split_comma(string)[-1] %>%
        as.numeric()
    mod_pos <- (mod_string_split + 1) %>%
        cumsum()
    mod_scores_split <- str_split_comma(scores) %>%
        as.numeric()

    data.frame(
        mod_pos = mod_pos,
        mod_prob = mod_scores_split/255
    )
}

cigar_tokeniser <- function(x) {
    non_empty_string <- function(x) { x != "" }
    state <- stringr::str_split(x, "\\d+") %>%
        unlist() %>%
        purrr::keep(non_empty_string)
    count <- stringr::str_split(x, "[[:alpha:]]") %>%
        unlist() %>%
        purrr::keep(non_empty_string)

    data.frame(state = state, count = as.numeric(count))
}

get_coord_map <- function(cigar) {
    tokens <- cigar_tokeniser_cpp(cigar)
    token_vec <- rep(tokens$state, times = tokens$count)
    seq_coord <- 0
    ref_coord <- 0

    seq_map <- numeric()
    ref_map <- numeric()

    for (tok in token_vec) {
        if (tok == "M") {
            seq_coord <- seq_coord + 1
            ref_coord <- ref_coord + 1
            seq_map <- c(seq_map, seq_coord)
            ref_map <- c(ref_map, ref_coord)
        } else if (tok %in% c("I", "S")) {
            seq_coord <- seq_coord + 1
            seq_map <- c(seq_map, seq_coord)
            ref_map <- c(ref_map, NA)
        } else if (tok %in% c("D", "N")) {
            ref_coord <- ref_coord + 1
        }
    }

    out <- setNames(ref_map, seq_map)
    out
}

modbam_to_ref_coord <- function(seq, cigar, mod_str, mod_scores, map_pos, strand) {
    mod_tokens <- mod_tokeniser_cpp(mod_str, mod_scores)
    coord_map <- get_coord_map_cpp(cigar)

    if (strand == "-") {
        mod_candidate_pos <- rev(get_char_pos_cpp(seq, "G"))[mod_tokens$mod_pos]
    } else if (strand == "+") {
        mod_candidate_pos <- get_char_pos_cpp(seq, "C")[mod_tokens$mod_pos]
    }

    rel_pos <- coord_map[mod_candidate_pos]
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

read_modbam_table <- function(x, chr, start, end, sample) {
    bam_file <- Rsamtools::BamFile(x)

    query <- GenomicRanges::GRanges(glue::glue("{chr}:{start}-{end}"))
    sb_param <- Rsamtools::ScanBamParam(
        what = c("qname", "rname", "strand", "pos", "cigar", "seq"),
        tag = c("MM", "ML"),
        which = query
    )

    reads <- Rsamtools::scanBam(bam_file, param = sb_param)

    parse_modbam <- function(x) {
        reads_df <- tibble::tibble(
            read_name = x$qname,
            chr = x$rname,
            strand = x$strand,
            map_pos = x$pos,
            cigar = x$cigar,
            seq = as.character(x$seq),
            mod_string = x$tag$MM,
            mod_scores = x$tag$ML %>% purrr::map_chr(~stringr::str_c(., collapse = ","))
        ) %>%
            dplyr::filter(seq != "")

        if (nrow(reads_df) == 0) {
            return(NULL)
        }

        out <- reads_df %>%
            mutate(
                modbam_stats = map_rows(
                    reads_df,
                    function(x) {
                        modbam_to_ref_coord(
                            x$seq,
                            x$cigar,
                            x$mod_string,
                            x$mod_scores,
                            x$map_pos,
                            x$strand
                        )
                    }
                )
           )

        out %>%
            dplyr::select("read_name", "chr", "strand", "modbam_stats") %>%
            tidyr::unnest("modbam_stats") %>%
            tidyr::drop_na() %>%
            dplyr::mutate(sample = sample, .before = 1)
    }

    lapply(reads, parse_modbam)
}
