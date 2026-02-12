# Convert exon annotation to genes

Convert exon annotation to genes

## Usage

``` r
exons_to_genes(x)
```

## Arguments

- x:

  the exon level annotation containing columns "gene_id", "chr",
  "strand" and "symbol".

## Value

the gene level annotation where each gene is taken to span the earliest
start position and latest end position of its exons.

## Examples

``` r
nmr <- load_example_nanomethresult()
#> Successfully matched 6 samples between data and annotation.
exons_to_genes(NanoMethViz::exons(nmr))
#> # A tibble: 6 × 6
#>   gene_id chr   strand symbol     start       end
#>   <chr>   <chr> <chr>  <chr>      <int>     <int>
#> 1 12189   chr11 -      Brca1  101488764 101551955
#> 2 12190   chr5  +      Brca2  150522630 150569747
#> 3 16210   chr18 +      Impact  12955852  12992950
#> 4 17263   chr12 +      Meg3   109541001 109571726
#> 5 18616   chr7  -      Peg3     6703892   6730431
#> 6 213742  chrX  -      Xist   103460366 103483254
```
