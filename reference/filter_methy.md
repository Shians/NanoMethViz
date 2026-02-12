# Create filtered methylation file

Create a filtered methylation file from an existing one.

## Usage

``` r
filter_methy(x, output_file, ...)
```

## Arguments

- x:

  the path to the methylation file or a NanoMethResult object.

- output_file:

  the output file to write results to (must end in .bgz).

- ...:

  filtering criteria given in dplyr syntax. Use methy_col_names() to get
  available column names.

## Value

invisibly returns 'output_file' if x is a file path, otherwise returns
NanoMethResult object with methy(x) replaced with filtered value.

## Examples

``` r
nmr <- load_example_nanomethresult()
#> Successfully matched 6 samples between data and annotation.
output_file <- paste0(tempfile(), ".tsv.bgz")
filter_methy(nmr, output_file = output_file, chr == "chrX")
#> 21,798 of 224,267 (9.72%) entries kept after filtering
#> results written to '/tmp/Rtmp2dJs99/file82ac1168b3f8.tsv.bgz' along with index file '/tmp/Rtmp2dJs99/file82ac1168b3f8.tsv.bgz.tbi'
filter_methy(methy(nmr), output_file = output_file, chr == "chrX")
#> 21,798 of 224,267 (9.72%) entries kept after filtering
#> results written to '/tmp/Rtmp2dJs99/file82ac1168b3f8.tsv.bgz' along with index file '/tmp/Rtmp2dJs99/file82ac1168b3f8.tsv.bgz.tbi'
```
