get_exons_mus_musculus <- function() {
    Mus.musculus <- getFromNamespace("Mus.musculus", "Mus.musculus")
    genes <-  AnnotationDbi::keys(Mus.musculus, "GENEID")
    exon_data <- AnnotationDbi::select(
        Mus.musculus,
        keys = genes,
        keytype = "GENEID",
        columns = c(
            "GENEID",
            "TXID",
            "EXONCHROM",
            "EXONSTRAND",
            "EXONSTART",
            "EXONEND",
            "SYMBOL"
        )
    )

    exon_tibble <- as_tibble(exon_data) %>%
        dplyr::rename(
            gene_id = GENEID,
            chr = EXONCHROM,
            strand = EXONSTRAND,
            start = EXONSTART,
            end = EXONEND,
            transcript_id = TXID,
            symbol = SYMBOL
        )
}
