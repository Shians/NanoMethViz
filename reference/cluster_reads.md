# Cluster reads based on methylation

Cluster reads based on methylation

## Usage

``` r
cluster_reads(x, chr, start, end, min_pts = 5)
```

## Arguments

- x:

  a ModBamResult object.

- chr:

  the chromosome name where to find the region.

- start:

  the start position of the region.

- end:

  the end position of the region.

- min_pts:

  the minimum number of points needed to form a cluster (default = 10).

## Value

A tibble with information about each read's cluster assignment and read
statistics.
