# Convert BSseq object to log-methylation-ratio matrix

Creates a log-methylation-ratio matrix from a BSseq object that is
useful for dimensionality reduction plots.

## Usage

``` r
bsseq_to_log_methy_ratio(
  bsseq,
  regions = NULL,
  prior_count = 2,
  drop_na = TRUE
)
```

## Arguments

- bsseq:

  the BSseq object.

- regions:

  the regions to calculate log-methylation ratios over. If left NULL,
  ratios will be calculated per site.

- prior_count:

  the prior count added to avoid taking log of 0.

- drop_na:

  whether to drop rows with all NA values.

## Value

a matrix containing log-methylation-ratios.

## Examples

``` r
nmr <- load_example_nanomethresult()
#> Successfully matched 6 samples between data and annotation.
bsseq <- methy_to_bsseq(nmr)
#> [2026-03-31 03:52:25] creating intermediate files...
#> [2026-03-31 03:52:25] parsing chr11...
#> [2026-03-31 03:52:25] parsing chr12...
#> [2026-03-31 03:52:25] parsing chr18...
#> [2026-03-31 03:52:25] parsing chr5...
#> [2026-03-31 03:52:25] parsing chr7...
#> [2026-03-31 03:52:25] parsing chrX...
#> [2026-03-31 03:52:25] samples found: B6Cast_Prom_3_cast B6Cast_Prom_3_bl6 B6Cast_Prom_2_cast B6Cast_Prom_2_bl6 B6Cast_Prom_1_cast B6Cast_Prom_1_bl6 
#> [2026-03-31 03:52:25] creating bsseq object...
#> [2026-03-31 03:52:25] reading in parsed data...
#> [2026-03-31 03:52:25] constructing matrices...
#> [2026-03-31 03:52:25] done
regions <- exons_to_genes(NanoMethViz::exons(nmr))
log_m_ratio <- bsseq_to_log_methy_ratio(bsseq, regions)
```
