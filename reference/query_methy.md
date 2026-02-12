# Query methylation data

Query methylation data

## Usage

``` r
query_methy(
  x,
  chr,
  start,
  end,
  simplify = TRUE,
  force = FALSE,
  truncate = TRUE,
  site_filter = getOption("NanoMethViz.site_filter", 3L)
)
```

## Arguments

- x:

  the NanoMethResult object or a path to the methylation data
  (tabix-bgzipped).

- chr:

  the vector of chromosomes

- start:

  the vector of start positions

- end:

  the vector of end positions

- simplify:

  whether returned results should be row-concatenated

- force:

  whether to force empty output when query region 'chr' does not appear
  in data. Without 'force', an empty result indicates that the requested
  'chr' appears in the data but no data overlaps with requested region,
  and an invalid 'chr' will cause an error.

- truncate:

  when querying from ModBamFiles, whether or not to truncate returned
  results to only those within the specified region. Otherwise
  methylation data for entire reads overlapping the region will be
  returned.

- site_filter:

  the minimum amount of coverage to report a site. This filters the
  queried data such that any site with less than the filter is not
  returned. The default is 1, which means that all sites are returned.
  This option can be set globally using the
  `options(NanoMethViz.site_filter = ...)` which will affect all
  plotting functions in NanoMethViz.

## Value

A table containing the data within the queried regions. If simplify is
TRUE (default) then returns all data in a single table, otherwise
returns a list of tables where each table is the data from one region.

## Details

The argument `site_filter` can be set globally using the
`options(NanoMethViz.site_filter = ...)` command. The same data entry
may appear multiple times in the output if it overlaps multiple regions.

## Examples

``` r
nmr <- load_example_nanomethresult()
#> Successfully matched 6 samples between data and annotation.
query_methy(methy(nmr), "chr7", 6703892, 6730431)
#> # A tibble: 16,473 × 7
#>    sample             chr       pos strand statistic read_name          mod_prob
#>    <chr>              <chr>   <int> <chr>      <dbl> <chr>                 <dbl>
#>  1 B6Cast_Prom_1_bl6  chr7  6704092 *          -1.05 d06f0069-470a-4ba…   0.259 
#>  2 B6Cast_Prom_1_bl6  chr7  6704092 *          -4.05 dd27e4e1-a012-494…   0.0171
#>  3 B6Cast_Prom_1_bl6  chr7  6704092 *           1.25 653d9fc1-f4a6-459…   0.777 
#>  4 B6Cast_Prom_1_bl6  chr7  6704092 *           1.39 6ea4506d-1ee3-4fc…   0.801 
#>  5 B6Cast_Prom_1_bl6  chr7  6704092 *           1.74 bf3446f8-af43-499…   0.851 
#>  6 B6Cast_Prom_1_bl6  chr7  6704092 *           2.16 d940c0fd-95c8-4e6…   0.897 
#>  7 B6Cast_Prom_1_bl6  chr7  6704092 *           4.31 bb55f1a3-015a-45b…   0.987 
#>  8 B6Cast_Prom_1_cast chr7  6704092 *          -3.83 88913f2c-cd8f-429…   0.0212
#>  9 B6Cast_Prom_1_cast chr7  6704092 *           1.12 b5f29356-1602-493…   0.754 
#> 10 B6Cast_Prom_1_cast chr7  6704092 *           1.36 2cf2b952-4091-48e…   0.796 
#> # ℹ 16,463 more rows
```
