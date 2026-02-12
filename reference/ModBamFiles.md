# Constructor for a ModBamFiles object

This function creates a ModBamFiles object containing information about
the samples and file paths. This constructor checks that the files are
readable and have an index.

## Usage

``` r
ModBamFiles(samples, paths)

# S4 method for class 'ModBamFiles'
show(object)
```

## Arguments

- samples:

  a character vector with the names of the samples.

- paths:

  a character vector with the file paths for the BAM files.

- object:

  a ModBamFiles object.

## Value

A ModBamFiles object with the sample and path information.
