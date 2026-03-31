# Convert NanoMethResult object to edgeR methylation matrix

Convert NanoMethResult object to edgeR methylation matrix

## Usage

``` r
methy_to_edger(methy, regions = NULL, out_folder = tempdir(), verbose = TRUE)
```

## Arguments

- methy:

  the NanoMethResult object or path to the methylation tabix file.

- regions:

  the regions to calculate log-methylation ratios over. If left NULL,
  ratios will be calculated per site.

- out_folder:

  the folder to store intermediate files. One file is created for each
  sample and contains columns "chr", "pos", "total" and "methylated".

- verbose:

  TRUE if progress messages are to be printed

## Value

a matrix compatible with the edgeR differential methylation pipeline

## Examples

``` r
nmr <- load_example_nanomethresult()
#> Successfully matched 6 samples between data and annotation.
edger_mat <- methy_to_edger(nmr)
#> [2026-03-31 03:53:53] creating intermediate files...
#> [2026-03-31 03:53:53] parsing chr11...
#> [2026-03-31 03:53:53] parsing chr12...
#> [2026-03-31 03:53:53] parsing chr18...
#> [2026-03-31 03:53:53] parsing chr5...
#> [2026-03-31 03:53:53] parsing chr7...
#> [2026-03-31 03:53:53] parsing chrX...
#> [2026-03-31 03:53:53] samples found: B6Cast_Prom_3_cast B6Cast_Prom_3_bl6 B6Cast_Prom_2_cast B6Cast_Prom_2_bl6 B6Cast_Prom_1_cast B6Cast_Prom_1_bl6 
#> [2026-03-31 03:53:53] creating bsseq object...
#> [2026-03-31 03:53:53] reading in parsed data...
#> [2026-03-31 03:53:53] constructing matrices...
#> [2026-03-31 03:53:53] done
```
