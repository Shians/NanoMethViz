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

# T2T ----
download_parse_and_save(
    "t2t",
    "https://hgdownload.soe.ucsc.edu/gbdb/hs1/bbi/cpgIslandExt.bb"
)

# T2T Genome ----
temp_path <- tempfile()
download.file("https://hgdownload.soe.ucsc.edu/gbdb/hs1/bbi/cpgIslandExt.bb", temp_path)

anno_t2t <- rtracklayer::import.bb(temp_path) %>%
    as_tibble()

cgi_anno_t2t <- anno_t2t %>%
    dplyr::rename(
        gene_id = name,
        chr = seqnames
    ) %>%
    mutate(
        transcript_id = gene_id,
        strand = "*",
        symbol = gene_id
    )

anno_name <- paste0("inst/cgi_t2t.rds")
saveRDS(exon_anno_t2t_formatted, anno_name, compress = "xz")

fs::file_delete(temp_path)
