# Query exons

Query a data.frame, NanoMethResult or ModBamResult for exon annotation.

## Usage

``` r
query_exons_region(x, chr, start, end)

query_exons_gene_id(x, gene_id)

query_exons_symbol(x, symbol)
```

## Arguments

- x:

  the object to query.

- chr:

  the chromosome to query.

- start:

  the start of the query region.

- end:

  the end of the query region.

- gene_id:

  the gene_id to query.

- symbol:

  the gene_id to query.

## Value

data.frame of queried exons.

## Functions

- `query_exons_region()`: Query region.

- `query_exons_gene_id()`: Query gene ID.

- `query_exons_symbol()`: Query gene symbol.
