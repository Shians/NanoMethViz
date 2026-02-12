# modBAM methylation results

A ModBamResult object stores modBAM data used for NanoMethViz
visualisation. It contains stores a ModBamFiles object, sample
information and optional exon information. The object is constructed
using the ModBamResult() constructor function described in "Usage".

## Usage

``` r
# S4 method for class 'ModBamResult'
methy(object)

# S4 method for class 'ModBamResult,ModBamFiles'
methy(object) <- value

# S4 method for class 'ModBamResult'
samples(object)

# S4 method for class 'ModBamResult,data.frame'
samples(object) <- value

# S4 method for class 'ModBamResult'
exons(object)

# S4 method for class 'ModBamResult,data.frame'
exons(object) <- value

# S4 method for class 'ModBamResult'
mod_code(object)

# S4 method for class 'ModBamResult,character'
mod_code(object) <- value

ModBamResult(methy, samples, exons = NULL, mod_code = "m")
```

## Arguments

- object:

  the ModBamResult object.

- value:

  the mod code.

- methy:

  a ModBamFiles object.

- samples:

  the data.frame of sample annotation containing at least columns sample
  and group.

- exons:

  (optional) the data.frame of exon information containing at least
  columns gene_id, chr, strand, start, end, transcript_id and symbol.

- mod_code:

  a character with the mod code of interest. Defaults to "m" for 5mC.
  See details for other options.

## Value

a ModBamResult object to be used with plotting functions

a ModBamFiles data.frame.

the sample annotation.

the exon annotation.

the mod code.

## Details

The possible tags for mod_code can be found at
<https://samtools.github.io/hts-specs/SAMtags.pdf> under the 'Base
modifications' section.

## Functions

- `methy(ModBamResult)`: modBAM information getter.

- `methy(object = ModBamResult) <- value`: modBAM information setter.

- `samples(ModBamResult)`: sample annotation getter.

- `samples(object = ModBamResult) <- value`: sample annotation setter.

- `exons(ModBamResult)`: exon annotation getter.

- `exons(object = ModBamResult) <- value`: exon annotation setter.

- `mod_code(ModBamResult)`: mod code getter.

- `mod_code(object = ModBamResult) <- value`: mod code setter.

- `ModBamResult()`: Constructor

## Slots

- `methy`:

  a ModBamFiles data.frame specifying the samples and paths to bam
  files.

- `samples`:

  the data.frame of sample annotation containing at least columns sample
  and group.

- `exons`:

  the data.frame of exon information containing at least columns
  gene_id, chr, strand, start, end, transcript_id and symbol.

- `mod_code`:

  the modification code of interest.
