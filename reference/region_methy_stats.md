# Calculate region methylation statistics

Calculate the average methylation probability and prevalence based on
specified probability threshold.

## Usage

``` r
region_methy_stats(nmr, regions, threshold = 0.5)
```

## Arguments

- nmr:

  the NanoMethResult object.

- regions:

  the table of regions to query statistics for.

- threshold:

  the threshold to use for determining methylation calls for the
  calculation of prevalence.

## Value

table of regions with additional columns of methylation summary
statistics.

## Examples

``` r
nmr <- load_example_nanomethresult()
#> Successfully matched 6 samples between data and annotation.
gene_anno <- exons_to_genes(NanoMethViz::exons(nmr))
region_methy_stats(nmr, gene_anno)
#> # A tibble: 6 × 8
#>   gene_id chr   strand symbol     start       end mean_methy_prob prevalence
#>   <chr>   <chr> <chr>  <chr>      <int>     <int>           <dbl>      <dbl>
#> 1 12189   chr11 -      Brca1  101488764 101551955           0.590      0.600
#> 2 12190   chr5  +      Brca2  150522630 150569747           0.557      0.565
#> 3 16210   chr18 +      Impact  12955852  12992950           0.586      0.594
#> 4 17263   chr12 +      Meg3   109541001 109571726           0.524      0.525
#> 5 18616   chr7  -      Peg3     6703892   6730431           0.618      0.628
#> 6 213742  chrX  -      Xist   103460366 103483254           0.520      0.519
```
