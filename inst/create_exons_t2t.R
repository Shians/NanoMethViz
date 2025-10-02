temp_path <- tempfile()
download.file("https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/annotation/chm13v2.0_RefSeq_Liftoff_v5.2.gff3.gz", temp_path)

gff_col_names <- c(
    "chr",
    "source",
    "type",
    "start",
    "end",
    "score",
    "strand",
    "phase",
    "attributes"
)

anno_t2t <- read_tsv(
    temp_path,
    col_names = gff_col_names,
    comment = "#"
)

anno_t2t <- read_tsv(
    temp_path,
    col_names = gff_col_names,
    comment = "#"
)

exon_anno_t2t <- anno_t2t %>%
    filter(type == "exon") %>%
    mutate(
        transcript_id = str_extract(attributes, "Parent=(\\w+)", group = 1),
        gene_id = str_extract(attributes, "GeneID:(\\w+)", group = 1),
        symbol = str_extract(attributes, "gene=(\\w+)", group = 1)
    )

exon_anno_t2t_formatted <- exon_anno_t2t %>%
    select(
        gene_id,
        chr,
        strand,
        start,
        end,
        transcript_id,
        symbol
    ) %>%
    mutate(
        gene_id = as.character(gene_id),
        chr = as.character(chr),
        strand = as.character(strand),
        start = as.integer(start),
        end = as.integer(end),
        transcript_id = as.character(transcript_id),
        symbol = as.character(symbol)
    ) %>%
    drop_na()

anno_name <- paste0("inst/exons_t2t.rds")
saveRDS(exon_anno_t2t_formatted, anno_name, compress = "xz")

fs::file_delete(temp_path)
