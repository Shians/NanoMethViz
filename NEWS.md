### Version 2.0.0
* Major changes to `plot_agg_regions()`.
  * Features of `plot_agg_regions()` and `plot_agg_regions_sample_grouped()` merged into one interface.
  * Regions now specified using single table.
* Changed `plot_regions()` default window proportion to 0.
* Changed default theme from `theme_bw()` to `theme_tufte()`.
* Added Megalodon data import instructions to "Importing Data" vignette.
* Added scico palette defaults for heatmaps. These are colourblind friendly.
* Added check for 0 length queries which would cause program to hang indefinitely.
* Added setters for NanoMethResult attributes `methy`, `samples` and `exons`.
* Added MDS and PCA plots.
* Added vignette for using external annotation and dimensionality reduction.
* Added binary thresholding for `plot_gene()`, `plot_region()` and `plot_agg_regions()`.
* Added regions argument to `bsseq_to_edger()` to calculate aggregate counts over features rather than per site.

### Version 1.1.4
* Added palette argument to aggregate plots
* Added `exons_to_genes()` function to convert exon annotation to gene annotation
* Added `plot_granges_heatmap()` function to use GRanges for plotting heatmaps

### Version 1.1.3

* Fixed group handling for list region input in `plot_agg_regions()`
* Fixed unused window size argument in `plot_region_heatmap()`
* Fixed error when reads overlap in name and position for internal function `StatLM()`

### Version 1.1.2

* Changed example dataset exon annotations from all genomic exons to just those contained in data.
* Fixed methylation heatmap to no longer be hard coded for Peg3.
* Added `plot_region_heatmap()` as analogue to `plot_region()`.
* Fixed `plot_agg_regions_sample_grouped()` to use `group` column of `NanoMethViz::samples(x)` rather than `haplotype`.
* Added unit tests.

### Version 1.1.1

* Added methylation heatmap via `plot_gene_heatmap()`.
* Fixed `gene_anno()` in `plot_gene()` for argument so FALSE actually turns off gene annotation.
* Added warning for cpp11 versions <0.2.5 which may cause memory crashes when trying to import methylation data.
* Added cpp11 version dependency to address tidyverse/readr#1145.
* Added query methylation by gene using `query_methy_gene()`.

### Version 1.0.0

* Initial Bioconductor release.
