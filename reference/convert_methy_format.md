# Convert methylation calls to NanoMethViz format

Convert methylation calls to NanoMethViz format

## Usage

``` r
convert_methy_format(
  input_files,
  output_file,
  samples = fs::path_ext_remove(fs::path_file(input_files)),
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

- verbose:

  TRUE if progress messages are to be printed

## Value

invisibly returns the output file path, creates a tabix file (.bgz) and
its index (.bgz.tbi)
