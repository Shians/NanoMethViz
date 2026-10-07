# Create BSSeq object from methylation tabix file

Create BSSeq object from methylation tabix file

## Usage

``` r
methy_to_bsseq(methy, out_folder = tempdir(), verbose = TRUE)
```

## Arguments

- methy:

  the NanoMethResult object or path to the methylation tabix file.

- out_folder:

  the folder to store intermediate files. One file is created for each
  sample and contains columns "chr", "pos", "total" and "methylated".

- verbose:

  TRUE if progress messages are to be printed

## Value

a BSSeq object.

## Examples

``` r
nmr <- load_example_nanomethresult()
#> Successfully matched 6 samples between data and annotation.
bsseq <- methy_to_bsseq(nmr)
#> [2026-10-07 23:18:09] creating intermediate files...
#> [2026-10-07 23:18:09] parsing chr11...
#> [2026-10-07 23:18:09] parsing chr12...
#> [2026-10-07 23:18:09] parsing chr18...
#> [2026-10-07 23:18:09] parsing chr5...
#> [2026-10-07 23:18:09] parsing chr7...
#> [2026-10-07 23:18:09] parsing chrX...
#> [2026-10-07 23:18:09] samples found: B6Cast_Prom_3_cast B6Cast_Prom_3_bl6 B6Cast_Prom_2_cast B6Cast_Prom_2_bl6 B6Cast_Prom_1_cast B6Cast_Prom_1_bl6 
#> [2026-10-07 23:18:09] creating bsseq object...
#> [2026-10-07 23:18:09] reading in parsed data...
#> [2026-10-07 23:18:09] constructing matrices...
#> [2026-10-07 23:18:09] done
```
