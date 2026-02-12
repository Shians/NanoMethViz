# Create a tabix file using methylation calls

Create a tabix file using methylation calls

## Usage

``` r
create_tabix_file(
  input_files,
  output_file,
  samples = extract_file_names(input_files),
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

## Examples

``` r
methy_calls <- system.file(package = "NanoMethViz",
    c("sample1_nanopolish.tsv.gz", "sample2_nanopolish.tsv.gz"), mustWork = FALSE)
temp_file <- paste0(tempfile(), ".tsv.bgz")

create_tabix_file(methy_calls, temp_file)
#> [2026-02-12 03:57:53] creating methylation table
#> processing /home/runner/work/_temp/Library/NanoMethViz/sample1_nanopolish.tsv.gz...
#> guessing file is produced by nanopolish...
#> processing /home/runner/work/_temp/Library/NanoMethViz/sample2_nanopolish.tsv.gz...
#> guessing file is produced by nanopolish...
#> [2026-02-12 03:57:54] sorting methylation table
#> [2026-02-12 03:57:54] compressing methylation table to tabix with index
```
