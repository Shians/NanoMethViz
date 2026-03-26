# Convert BAM with modifications to tabix format

The `modbam_to_tabix` function takes a ModBamResult object and converts
it into a tabix file format, which is efficient for indexing and
querying large datasets.

## Usage

``` r
modbam_to_tabix(x, out_file, mod_code = NanoMethViz::mod_code(x))
```

## Arguments

- x:

  the `ModBamResult` object.

- out_file:

  the path of the output tabix.

- mod_code:

  the modification code to use, defaults to 'm' for 5mC methylation.

## Value

invisibly returns the name of the created tabix file.

## Details

The possible tags for mod_code can be found at
<https://samtools.github.io/hts-specs/SAMtags.pdf> under the 'Base
modifications' section.

## Examples

``` r
out_file <- paste0(tempfile(), ".tsv.bgz")
mbr <- ModBamResult(
    methy = ModBamFiles(
        samples = "sample1",
        paths = system.file("peg3.bam", package = "NanoMethViz",
        mustWork = FALSE)
    ),
    samples = data.frame(
        sample = "sample1",
        group = "group1"
    )
)
#> Successfully created ModBamResult with 1 matched samples.

modbam_to_tabix(mbr, out_file)
#> ℹ Writing data to temporary file: /tmp/RtmpGrPysM/file833243de0342.tsv
#> ✔ Writing data to temporary file: /tmp/RtmpGrPysM/file833243de0342.tsv [7ms]
#> 
#> ℹ Converting data to TSV
#> ✔ Converting data to TSV [157ms]
#> 
#> ℹ Sorting data
#> ✔ Sorting data [30ms]
#> 
#> ℹ Compressing data
#> ℹ Moving data to final location: /tmp/RtmpGrPysM/file833252490891.tsv.bgz
#> ℹ Compressing data
#> ✔ Compressing data [28ms]
#> 
#> ℹ Tabix file created: /tmp/RtmpGrPysM/file833252490891.tsv.bgz
#> ✔ Tabix file created: /tmp/RtmpGrPysM/file833252490891.tsv.bgz [12ms]
#> 
```
