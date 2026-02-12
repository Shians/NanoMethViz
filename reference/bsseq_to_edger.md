# Convert BSseq object to edgeR methylation matrix

Convert BSseq object to edgeR methylation matrix

## Usage

``` r
bsseq_to_edger(bsseq, regions = NULL)
```

## Arguments

- bsseq:

  the BSseq object.

- regions:

  the regions to calculate log-methylation ratios over. If left NULL,
  ratios will be calculated per site.

## Value

a matrix compatible with the edgeR differential methylation pipeline

## Examples

``` r
methy <- system.file("methy_subset.tsv.bgz", package = "NanoMethViz", mustWork = FALSE)
bsseq <- methy_to_bsseq(methy)
#> [2026-02-12 03:57:48] creating intermediate files...
#> [2026-02-12 03:57:48] parsing chr11...
#> [2026-02-12 03:57:48] parsing chr12...
#> [2026-02-12 03:57:48] parsing chr18...
#> [2026-02-12 03:57:48] parsing chr5...
#> [2026-02-12 03:57:48] parsing chr7...
#> [2026-02-12 03:57:48] parsing chrX...
#> [2026-02-12 03:57:48] samples found: B6Cast_Prom_3_cast B6Cast_Prom_3_bl6 B6Cast_Prom_2_cast B6Cast_Prom_2_bl6 B6Cast_Prom_1_cast B6Cast_Prom_1_bl6 
#> [2026-02-12 03:57:48] creating bsseq object...
#> [2026-02-12 03:57:48] reading in parsed data...
#> [2026-02-12 03:57:48] constructing matrices...
#> [2026-02-12 03:57:48] done
edger_mat <- bsseq_to_edger(bsseq)
```
