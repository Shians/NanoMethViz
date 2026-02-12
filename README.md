# NanoMethViz

<!-- badges: start -->
[![Codecov test coverage](https://codecov.io/gh/Shians/NanoMethViz/branch/master/graph/badge.svg)](https://codecov.io/gh/Shians/NanoMethViz?branch=master)
[![R-CMD-check](https://github.com/Shians/NanoMethViz/workflows/R-CMD-check/badge.svg)](https://github.com/Shians/NanoMethViz/actions)
<!-- badges: end -->

NanoMethViz is a Bioconductor package for visualising and summarising DNA methylation from Oxford Nanopore long-read sequencing.

It supports:

- tabular methylation calls (for example, Modkit, Megalodon, Nanopolish, and f5c)
- modBAM-based methylation calls

## Installation

Install from Bioconductor:

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}
BiocManager::install("NanoMethViz")
```

Install the Bioconductor devel version:

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}
BiocManager::install(version = "devel")
BiocManager::install("NanoMethViz")
```

## Quick Start

For tabix-backed methylation data, `NanoMethResult` needs three inputs:

- `methy`: path to a bgzip-compressed and tabix-indexed methylation file
- `samples`: a `data.frame` with at least a `sample` column
- `exons`: a `data.frame` with exon-level annotation

```r
library(NanoMethViz)

# Example object bundled with the package
nmr <- load_example_nanomethresult()

# Region/gene visualisation
plot_gene(nmr, "Peg3")

# Group-level summaries
plot_agg_genes(nmr)

# Dimension reduction on methylation profiles
bss <- methy_to_bsseq(nmr)
gene_regions <- exons_to_genes(exons(nmr))
lmr <- bsseq_to_log_methy_ratio(bss, gene_regions)
plot_mds(lmr)
plot_pca(lmr)
```

## Typical Workflow

1. Import methylation calls.
2. Convert to tabix format if needed (`convert_methy_format()`, `create_tabix_file()`, `modbam_to_tabix()`).
3. Prepare sample metadata and exon annotation (`get_exons_hg38()`, `get_exons_mm10()`, or custom annotation).
4. Build a result object (`NanoMethResult()` or `ModBamResult()`).
5. Query and visualise methylation (`query_methy()`, `plot_gene()`, `plot_region()`, `plot_agg_genes()`, `plot_agg_regions()`, `plot_violin()`, `plot_region_heatmap()`).
6. Export for downstream differential analysis (`methy_to_bsseq()`, `bsseq_to_edger()`).

## Common Tasks

Choose the right object type:

- Use `NanoMethResult` for tabix-backed tabular methylation calls.
- Use `ModBamResult` for modBAM-based workflows.

If you want to...

- plot one gene: use `plot_gene(x, gene)`
- plot one genomic interval: use `plot_region(x, chr, start, end)`
- show per-read heatmap in a region: use `plot_region_heatmap(x, chr, start, end)`
- aggregate over genes/regions: use `plot_agg_genes(x)` / `plot_agg_regions(x, regions)`
- get methylation values for a region: use `query_methy(x, chr, start, end)`
- run dimension reduction: use `plot_mds(log_methy_ratio_matrix)` / `plot_pca(log_methy_ratio_matrix)`
- export for differential testing: use `methy_to_bsseq(x)`, then `bsseq_to_edger(...)`

## Core Data Requirements

- Chromosome naming must match across methylation calls and annotation (`chr1` vs `1`).
- Sample names in methylation data must match `samples$sample`.
- Exon annotation should include columns: `gene_id`, `chr`, `strand`, `start`, `end`, `transcript_id`, `symbol`

## Key Global Options

- `options(NanoMethViz.site_filter = 3L)`: minimum site coverage used by query/plot functions.
- `options(NanoMethViz.highlight_col = "grey50")`: default colour for highlighted annotation regions.

## Documentation

- User Guide (Bioconductor vignette): https://www.bioconductor.org/packages/release/bioc/vignettes/NanoMethViz/inst/doc/UsersGuide.html
- User Guide (pkgdown site): https://shians.github.io/NanoMethViz/articles/UsersGuide.html
- Function reference index: https://shians.github.io/NanoMethViz/reference/index.html
- Issue tracker: https://github.com/Shians/NanoMethViz/issues

## Example Plots

### MDS plot

![](man/figures/mds.png)

### Feature aggregation

![](man/figures/agg_genes.png)

### Region methylation plot

![](man/figures/peg3_gene.png)

## License

This project is licensed under Apache License, Version 2.0.
