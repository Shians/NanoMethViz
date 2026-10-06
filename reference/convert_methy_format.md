# Convert methylation calls to NanoMethViz format

Convert methylation calls to NanoMethViz format

## Usage

``` r
convert_methy_format(
  input_files,
  output_file,
  samples = extract_file_names(input_files),
  mod_code = NULL,
  verbose = TRUE
)
```

## Arguments

- input_files:

  the files to convert

- output_file:

  the output file to write results to (must end in .bgz)

- samples:

  the names of samples corresponding to each file

- mod_code:

  the modification code to extract from modkit input. NULL uses "m"
  (5mC). Must be NULL unless at least one input is from modkit.

- verbose:

  TRUE if progress messages are to be printed

## Value

invisibly returns the output file path, creates a tabix file (.bgz) and
its index (.bgz.tbi)
