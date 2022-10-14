### Version 2.4.0
* Fixed `plot_region_heatmap()` producing the wrong plot when a factor is used for the chromosome.
* Fixed nanopolish and f5c import positions being off by 1.
* Fixed broken `samples()` setter for NanoMethResults.
* Added `plot_agg_genes()` function as a shorthand for `plot_agg_regions(x, exons_to_genes(exons(x)))`. 
* Added the ability to interrupt `methy_to_bsseq()` calls.
* Added handling for NanoMethResults objects in `filter_methy()`. If NanoMethResult is used as input, then NanoMethResult is invisibly returned as output.
* Added black outlines to exons in annotation to distinguish contiguous segments for features like tandem repeats.
* Added `line_size` argument to `plot_gene()`, `plot_region()` and `plot_granges()` plots for adjusting line size.
* Added `subsample` argument to heatmap plots, default 50. This reduces the number of rows shown the plot to the specified amount.
* Added `get_exons_mm10()`, `get_exons_hg19()`, and `get_exons_hg38()` as replacements for `get_exons_mus_musculus()` and `get_exons_homo_sapiens()`.
* Changed heatmaps to no longer plot samples that are absent from sample annotations.
* Changed heatmap labels to appear on the right rather than on top.
* Changed heatmap alpha from 0.33 to 0.5.
* Changed arrows in exon connectors to appear in the middle as open arrow instead of at the end as closed arrow.
* Changed default X axis labels to be rescaled to appropriate SI-style. e.g. Kb, Mb, Gb.

### Version 2.2.0
* Added `heatmap` argument to `plot_gene()`, `plot_region()` and `plot_granges()`. This adds a read-heatmap to the plot.
* Added `cluster_regions()` function to perform k-means clustering on a table of genomic regions based on methylation profile.
* Added median averaging method for trends in `plot_gene()`, `plot_region()` and `plot_granges()`. This can be changed using the new `avg_method` argument, default is `mean`.
* Added `filter_methy()` function to create a filtered methylation file.
* Added `region_methy_stats()` to obtain average methylation fractions of specific regions.
* Added `methy_to_edger()` direct conversion wrapper around `methy_to_bsseq()` and `bsseq_to_edger()`.
* Added `palette` argument to `plot_gene()`, `plot_region()` and `plot_granges()` to allow custom colour palettes.
* Fixed `bsseq_to_edger()` failing when regions argument was used.
* Fixed heatmaps not staying in a single column when more than 2 groups were present.

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
