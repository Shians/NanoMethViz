library(NanoMethViz)
library(bsseq)
library(stringr)

bss <- nmr %>% methy_to_bsseq()
samples <- colData(bss)$sample %>%
    str_replace(".*_(\\d)_(bl6|cast)$", "\\2 \\1")
groups <- str_extract(samples, "(bl6|cast)")

ggsave(
    "img/mds.png",
    plot_mds(
        bsseq_to_log_methy_ratio(
            bss,
            regions = exons_to_genes(NanoMethViz::exons(nmr))
        ),
        labels = samples,
        groups = groups
    ) +
        ggtitle("MDS plot"),
    height = 600, width = 800, units = "px", dpi = 150
)


ggsave(
    "img/peg3_spaghetti.png",
    plot_gene(nmr, "Peg3", spaghetti = TRUE),
    height = 600, width = 800, units = "px", dpi = 150)

ggsave(
    "img/peg3_heatmap.png",
    plot_gene_heatmap(nmr, "Peg3") +
        ggtitle("Peg3"),
    height = 600, width = 800, units = "px", dpi = 150)

ggsave(
    "img/agg_genes.png",
    plot_agg_regions(
        nmr,
        regions = exons_to_genes(NanoMethViz::exons(nmr))
    ) +
        ggtitle("Aggregated profile over genes"),
    height = 600, width = 800, units = "px", dpi = 150)

