# Plot gene methylation

Plot the methylation of a gene symbol specified within the exon(x) slot.

## Usage

``` r
plot_gene(x, gene, ...)

# S4 method for class 'NanoMethResult,character'
plot_gene(
  x,
  gene,
  window_prop = 0.3,
  anno_regions = NULL,
  binary_threshold = NULL,
  avg_method = c("mean", "median"),
  spaghetti = FALSE,
  heatmap = TRUE,
  heatmap_subsample = 50,
  smoothing_window = 2000,
  gene_anno = TRUE,
  palette = ggplot2::scale_colour_brewer(palette = "Set1"),
  line_size = 1,
  mod_scale = c(0, 1),
  span = NULL
)

# S4 method for class 'ModBamResult,character'
plot_gene(
  x,
  gene,
  window_prop = 0.3,
  anno_regions = NULL,
  binary_threshold = NULL,
  avg_method = c("mean", "median"),
  spaghetti = FALSE,
  heatmap = TRUE,
  heatmap_subsample = 50,
  smoothing_window = 2000,
  gene_anno = TRUE,
  palette = ggplot2::scale_colour_brewer(palette = "Set1"),
  line_size = 1,
  mod_scale = c(0, 1),
  span = NULL
)
```

## Arguments

- x:

  the NanoMethResult or ModBamResult object.

- gene:

  the gene symbol for the gene to plot.

- ...:

  additional arguments.

- window_prop:

  the size of flanking region to plot. Can be a vector of two values for
  left and right window size. Values indicate proportion of gene length.

- anno_regions:

  the data.frame of regions to be annotated.

- binary_threshold:

  the modification probability such that calls with modification
  probability above the threshold are set to 1 and probabilities equal
  to or below the threshold are set to 0.

- avg_method:

  the average method for pre-smoothing at each genomic position. Data is
  pre-smoothed at each genomic position before the smoothed aggregate
  line is generated for performance reasons. The default is "mean" which
  corresponds to the average methylation fraction. The alternative
  "median" option is closer to an average within the more common
  methylation state.

- spaghetti:

  whether or not individual reads should be shown.

- heatmap:

  whether or not read-methylation heatmap should be shown.

- heatmap_subsample:

  how many packed rows of reads to subsample to.

- smoothing_window:

  the window size for smoothing the trend line.

- gene_anno:

  whether to show gene annotation.

- palette:

  the ggplot colour palette used for groups.

- line_size:

  the size of the lines.

- mod_scale:

  the scale range for modification probabilities. Default c(0, 1), set
  to "auto" for automatic limits.

- span:

  DEPRECATED, use smoothing_window instead. Will be removed in next
  version.

## Value

a patchwork plot containing the methylation profile in the specified
region.

## Details

This function plots the methylation data for a given gene. The main
trendline plot shows the average methylation probability across the
gene. The heatmap plot shows the methylation probability for each read
across the gene. The gene annotation plot shows the exons of the gene.
In the heatmap, each row represents one or more non-overlapping reads
where the coloured segments represent the methylation probability at
each position. Data along a read is connected by a grey line. The gene
annotation plot shows the isoforms and exons of the gene, with arrows
indicating the direction of transcription.

Since V3.0.0 NanoMethViz has changed the smoothing strategy from a loess
smoothing to a weighted moving average. This is because the loess
smoothing was too computationally expensive for large datasets and had a
span parameter that was difficult to tune. The new smoothing strategy is
controlled by the smoothing_window argument.

## Functions

- `plot_gene(x = ModBamResult, gene = character)`: S4 method for
  ModBamResult

## Examples

``` r
nmr <- load_example_nanomethresult()
#> Successfully matched 6 samples between data and annotation.
plot_gene(nmr, "Peg3")
#> `height` was translated to `width`.

```
