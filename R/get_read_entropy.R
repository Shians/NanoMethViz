get_read_entropy <- function(bam_path, sample = fs::path_file(bam_path)) {
    # helper functions ----
    count_crossings <- function(numbers) {
        crossings <- sum(((numbers[-1] - 127.5) * (numbers[-length(numbers)] - 127.5)) < 0)
        return(crossings)
    }

    read_chunk <- function(bam_file) {
        Rsamtools::scanBam(bam_file, param = modbam_param())[[1]]
    }

    parse_chunk <- function(reads, sample) {
        reads <- reads[!is.null(reads$seq)]
        tibble::tibble(
            sample = factor(sample),
            read_name = reads$qname,
            chr = factor(reads$rname),
            pos = reads$pos,
            width = Biostrings::width(reads$seq),
            center = as.integer(pos + width/2),
            cg_count = purrr::map_int(as.character(reads$seq), count_cg_cpp),
            crossings = purrr::map_int(reads$tag$ML, count_crossings),
            entropy = crossings / cg_count
        )
    }

    # function body ----
    fname <- fs::path_file(bam_path)
        cli::cli_progress_bar(
            glue::glue("Parsing file: {fname}"),
            total = get_bam_total_reads(bam_path),
            format_done = paste0(
                "{.alert-success Data parsed: ", fname, " {.timestamp {cli::pb_elapsed}}}"),
            format_failed = paste0(
                "{.alert-danger Data parsing failed: ", fname, " {.timestamp {cli::pb_elapsed}}}"),
            clear = FALSE
        )

    bam_file <- Rsamtools::BamFile(bam_path, yieldSize = 15000)

    i <- 1
    df_list <- list()
    open(bam_file)
    while (Rsamtools::isIncomplete(bam_file)) {
        reads <- read_chunk(bam_file)
        df_list[[i]] <- parse_chunk(reads, sample)
        i <- i + 1
        cli::cli_progress_update(length(reads[[1]]))
    }
    close(bam_file)

    bind_rows(df_list)
}
