# Validate sample annotation data frame

Validate sample annotation data frame

## Usage

``` r
assert_valid_samples(
  samples,
  required_cols = c("sample", "group"),
  context = "sample annotation"
)
```

## Arguments

- samples:

  sample annotation data.frame

- required_cols:

  required column names

- context:

  description of where this is being used for error messages
