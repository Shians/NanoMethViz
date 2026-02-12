# Assert that file paths are readable

This function checks whether all provided file paths exist and are
readable. If any file paths do not exist, it throws an informative error
message.

## Usage

``` r
assert_valid_genomic_coords(chr, start, end, allow_equal = FALSE)
```

## Arguments

- x:

  A character vector of file paths to check for existence

## Value

Nothing if all files exist, otherwise throws an error
