# Nanopore Methylation Result

A NanoMethResult object stores data used for NanoMethViz visualisation.
It contains stores a path to the methylation data, sample information
and optional exon information. The object is constructed using the
NanoMethResult() constructor function described in "Usage".

## Usage

``` r
NanoMethResult(methy, samples, exons = NULL)

# S4 method for class 'NanoMethResult'
methy(object)

# S4 method for class 'NanoMethResult,ANY'
methy(object) <- value

# S4 method for class 'NanoMethResult'
samples(object)

# S4 method for class 'NanoMethResult,data.frame'
samples(object) <- value

# S4 method for class 'NanoMethResult'
exons(object)

# S4 method for class 'NanoMethResult,data.frame'
exons(object) <- value
```

## Arguments

- methy:

  the path to the methylation tabix file.

- samples:

  the data.frame of sample annotation containing at least columns sample
  and group.

- exons:

  (optional) the data.frame of exon information containing at least
  columns gene_id, chr, strand, start, end, transcript_id and symbol.

- object:

  the NanoMethResult object.

- value:

  the exon annotation.

## Value

a NanoMethResult object to be used with plotting functions

the path to the methylation data.

the sample annotation.

the exon annotation.

## Functions

- `NanoMethResult()`: Constructor

- `methy(NanoMethResult)`: methylation data path getter.

- `methy(object = NanoMethResult) <- value`: methylation data path
  setter.

- `samples(NanoMethResult)`: sample annotation getter.

- `samples(object = NanoMethResult) <- value`: sample annotation setter.

- `exons(NanoMethResult)`: exon annotation getter.

- `exons(object = NanoMethResult) <- value`: exon annotation setter.

## Slots

- `methy`:

  the path to the methylation tabix file.

- `samples`:

  the data.frame of sample annotation containing at least columns sample
  and group.

- `exons`:

  the data.frame of exon information containing at least columns
  gene_id, chr, strand, start, end, transcript_id and symbol.

## Examples

``` r
methy <- system.file(package = "NanoMethViz", "methy_subset.tsv.bgz", mustWork = FALSE)
sample <- c(
    "B6Cast_Prom_1_bl6",
    "B6Cast_Prom_1_cast",
    "B6Cast_Prom_2_bl6",
    "B6Cast_Prom_2_cast",
    "B6Cast_Prom_3_bl6",
    "B6Cast_Prom_3_cast"
)
group <- c(
    "bl6",
    "cast",
    "bl6",
    "cast",
    "bl6",
    "cast"
)
sample_anno <- data.frame(sample, group, stringsAsFactors = FALSE)
exon_tibble <- get_example_exons_mus_musculus()
NanoMethResult(methy, sample_anno, exon_tibble)
#> Successfully matched 6 samples between data and annotation.
#> An object of class "NanoMethResult"
#> Slot "methy":
#> [1] "/home/runner/work/_temp/Library/NanoMethViz/methy_subset.tsv.bgz"
#> 
#> Slot "samples":
#> # A tibble: 6 × 2
#>   sample             group
#>   <chr>              <fct>
#> 1 B6Cast_Prom_1_bl6  bl6  
#> 2 B6Cast_Prom_1_cast cast 
#> 3 B6Cast_Prom_2_bl6  bl6  
#> 4 B6Cast_Prom_2_cast cast 
#> 5 B6Cast_Prom_3_bl6  bl6  
#> 6 B6Cast_Prom_3_cast cast 
#> 
#> Slot "exons":
#> # A tibble: 303 × 7
#>    gene_id chr   strand     start       end transcript_id symbol
#>    <chr>   <chr> <chr>      <int>     <int>         <int> <chr> 
#>  1 12189   chr11 -      101551769 101551879         92234 Brca1 
#>  2 12189   chr11 -      101549015 101549113         92234 Brca1 
#>  3 12189   chr11 -      101539981 101540034         92234 Brca1 
#>  4 12189   chr11 -      101535528 101535605         92234 Brca1 
#>  5 12189   chr11 -      101533939 101534027         92234 Brca1 
#>  6 12189   chr11 -      101532020 101532159         92234 Brca1 
#>  7 12189   chr11 -      101530992 101531094         92234 Brca1 
#>  8 12189   chr11 -      101529791 101529836         92234 Brca1 
#>  9 12189   chr11 -      101528001 101528077         92234 Brca1 
#> 10 12189   chr11 -      101523325 101526639         92234 Brca1 
#> # ℹ 293 more rows
#> 

x <- load_example_nanomethresult()
#> Successfully matched 6 samples between data and annotation.
methy(x)
#> [1] "/home/runner/work/_temp/Library/NanoMethViz/methy_subset.tsv.bgz"
```
