# UCSC Genomes ----

# exons
# the data.frame of exon information containing at least columns gene_id, chr, strand, start, end, transcript_id and symbol.

col_names <- c(
    "bin",
    "chrom",
    "chromStart",
    "chromEnd",
    "name",
    "length",
    "cpgNum",
    "gcNum",
    "perCpg",
    "perGc",
    "obsExp"
)

read_cgi_anno <- function(x) {
    x %>%
        read_tsv(col_names = col_names) %>%
        dplyr::rename(
            gene_id = name,
            chr = chrom,
            start = chromStart,
            end = chromEnd
        ) %>%
        mutate(
            # UCSC starts are 0-based, the package convention is 1-based
            start = start + 1,
            transcript_id = gene_id,
            strand = "*",
            symbol = gene_id
        )
}

download_parse_and_save <- function(genome_name, url) {
    temp_path <- tempfile()
    download.file(url, temp_path)

    anno_name <- paste0("inst/cgi_", genome_name, ".rds")
    saveRDS(read_cgi_anno(temp_path), anno_name, compress = "xz")

    fs::file_delete(temp_path)
}

# mm10 ----
download_parse_and_save(
    "mm10",
    "https://hgdownload.soe.ucsc.edu/goldenPath/mm10/database/cpgIslandExt.txt.gz"
)

# GRCm39 ----
download_parse_and_save(
    "GRCm39",
    "https://hgdownload.soe.ucsc.edu/goldenPath/mm39/database/cpgIslandExt.txt.gz"
)

# hg19 ----
download_parse_and_save(
    "hg19",
    "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/cpgIslandExt.txt.gz"
)

# hg38 ----
download_parse_and_save(
    "hg38",
    "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/cpgIslandExt.txt.gz"
)

# T2T (hs1) ----
# UCSC only distributes hs1 CpG islands as a bigBed, which uses GenBank
# accessions as sequence names and has no bin column. Rename sequences to UCSC
# names with the chromAlias table and recompute the UCSC bin so the table
# matches the other genomes' cpgIslandExt layout. import.bb() already returns
# 1-based starts.

# UCSC standard binning scheme (binFromRange in kent/src/lib/binRange.c)
ucsc_bin <- function(start, end) {
    bin_offsets <- c(512 + 64 + 8 + 1, 64 + 8 + 1, 8 + 1, 1, 0)
    start_bin <- start %/% 2^17
    end_bin <- (end - 1) %/% 2^17
    bin <- rep(NA_real_, length(start))
    for (offset in bin_offsets) {
        hit <- is.na(bin) & start_bin == end_bin
        bin[hit] <- offset + start_bin[hit]
        start_bin <- start_bin %/% 2^3
        end_bin <- end_bin %/% 2^3
    }
    bin
}

bb_path <- tempfile(fileext = ".bb")
alias_path <- tempfile(fileext = ".txt")
download.file("https://hgdownload.soe.ucsc.edu/gbdb/hs1/bbi/cpgIslandExt.bb", bb_path)
download.file("https://hgdownload.soe.ucsc.edu/goldenPath/hs1/bigZips/hs1.chromAlias.txt", alias_path)

chrom_alias <- read_tsv(alias_path, comment = "#", col_names = FALSE) %>%
    dplyr::select(ucsc = X1, genbank = X2)

cgi_anno_t2t <- rtracklayer::import.bb(bb_path) %>%
    as_tibble() %>%
    mutate(
        chr = chrom_alias$ucsc[match(as.character(seqnames), chrom_alias$genbank)],
        start = as.numeric(start),
        end = as.numeric(end),
        bin = ucsc_bin(start - 1, end),
        across(c(length, cpgNum, gcNum), as.numeric)
    ) %>%
    dplyr::select(bin, chr, start, end, gene_id = name, length, cpgNum, gcNum, perCpg, perGc, obsExp) %>%
    mutate(
        transcript_id = gene_id,
        strand = "*",
        symbol = gene_id
    )

stopifnot(!anyNA(cgi_anno_t2t$chr))
saveRDS(cgi_anno_t2t, "inst/cgi_t2t.rds", compress = "xz")

fs::file_delete(c(bb_path, alias_path))
