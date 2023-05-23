# get modification probability values and reference-based coordinates
modbam_to_ref_coord <- function(seq, cigar, mod_str, mod_scores, map_pos, strand) {
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


read_modbam_table <- function(x, chr, start, end, sample) {
    bam_file <- Rsamtools::BamFile(x)

    query <- make_granges(chr, start, end)

    if (length(query) == 0) {
        return(list(empty_methy_query_output()))
    }

    sb_param <- Rsamtools::ScanBamParam(
        what = c("qname", "rname", "strand", "pos", "cigar", "seq"),
        tag = c("MM", "ML"),
        which = query
    )

    reads <- Rsamtools::scanBam(bam_file, param = sb_param)

    parse_modbam <- function(x) {
        if (is.null(x$tag$ML)) {
            mod_scores <- NULL
        } else {
            mod_scores <- x$tag$ML %>%
                purrr::map_chr(~stringr::str_c(., collapse = ","))
        }

        reads_df <- tibble::tibble(
            read_name = x$qname,
            chr = x$rname,
            strand = x$strand,
            map_pos = x$pos,
            cigar = x$cigar,
            seq = as.character(x$seq),
            mod_string = x$tag$MM,
            mod_scores = mod_scores
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
