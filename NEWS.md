## Version 3.4.0
* Fixed `modbam_to_tabix()` to not delete existing directories when provided as the output file. Previously it would delete the directory and create a file with the same name, now it will fail if the output is a directory.

## Version 3.2.0
* Added preliminary support for importing modkit data using `create_tabix_file()`.

## Version 3.1.0
* Fixed parsing error for BAM files leading to an extre site to be called whenever skip width is greater than 0.
* Added functions `get_cgi_*()` to get CpG islands annotation for mm10, GRCm39, hg19, and hg38.
* Changed colour palette for heatmaps to have a less bright yellow on the low end.
* Changed default NanoMethViz.site_filter value to 3. This filters out any sites with less than 3 coverage.
* Changed BAM file parsing behaviour when parse flag "." or "?" is not set to default to "?". This is against SAM spec but follows the behaviour of IGV and older PacBio datasets.

## Version 3.0.0
* Breaking change to smoothing strategy in plot_gene(), plot_region(), and plot_granges() to use weighted moving mean instead of loess. This deprecates the `span` argument in favour of `smoothing_window` which is defaulted to 2000 bases.
  * Smoothing for various plotting functions was previously performed using loess smoothing, this performed locally weighted linear estimation to create a smoothed line, the span argument controlled the proportion of data used in this smoothing. This parameter was difficult to tune because under a fixed span, the smoothed line became flatter as the plot region grew larger. Internally, NanoMethViz dynamically calculated a span that changed inversely proportional to the width of the plotting region, decreasing the span as the plot region grew. However the calculated span was invisible to users, and it unintuitive to users how to set a span to change the appearance of the plot.
  * Changing the smoothing method to a weighted rolling mean lead to the new `smoothing_window` argument which represents the window size in bases from which data is used for smoothing around each point. This serves the same purpose as the dynamic calculation done previously, but is set more explicitly and should be more intuitive for users. The default is always 2000, and can be increased to increase smoothness and decreased to decrease smoothness.
* Breaking change to the appearance of `plot_gene()` plots, previously the isoform annotation would be restricted to only the gene of interest. It is now changed to follow the same behaviour as `plot_region()` whereby all isoforms in the region are plotted.
* Breaking change to the default plotting options for `plot_region()`, `plot_gene()` and `plot_grange()` to plot heatmap by default.
* Possible breaking change to `query_methy(simplify = FALSE)`, it will now return a list that is the same length as the number of regions queried, where it previously returned nothing if a particular sequence was missing from the tabix.
* Added `gene_anno` argument to `plot_region()` and `plot_granges()` to control whether gene annotation is plotted.
* Added plot_violin() function for creating violin plots for samples over specific regions.
* Added check to remove hard-clipped reads because they may not have matching mod strings.
* Changed gene annotation to always put label on visible isoforms, previously labels are plotted at the center of isoform.
* Changed the order of columns when querying from ModBam to be the same as when querying from Tabix, with readname at the end instead of being the second column.
* Changed `query_methy()` to always return same length output as the query when simplify = FALSE.
* Fixed memory leak in bam parsing due to out of bounds access.
* Fixed crash when CIGAR doesn't match length of SEQ.

## Version 2.6.0
* Added preliminary modbam file support.
* Changed rug plot to appear under other geoms. This helps with visibility of data when methylation values are close to 0.
* Changed heatmap alpha from 0.5 to 1, line width from 1.0 to 1.2 and line colour from black to darkgrey.
* Changed x-axis limits on plots to be controlled using coord_cartesian instead of scale_x_continuous. Plots should now accurately represent data around the boundaries.

## Version 2.4.0
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

## Version 2.2.0
* Added `heatmap` argument to `plot_gene()`, `plot_region()` and `plot_granges()`. This adds a read-heatmap to the plot.
* Added `cluster_regions()` function to perform k-means clustering on a table of genomic regions based on methylation profile.
* Added median averaging method for trends in `plot_gene()`, `plot_region()` and `plot_granges()`. This can be changed using the new `avg_method` argument, default is `mean`.
* Added `filter_methy()` function to create a filtered methylation file.
* Added `region_methy_stats()` to obtain average methylation fractions of specific regions.
* Added `methy_to_edger()` direct conversion wrapper around `methy_to_bsseq()` and `bsseq_to_edger()`.
* Added `palette` argument to `plot_gene()`, `plot_region()` and `plot_granges()` to allow custom colour palettes.
* Fixed `bsseq_to_edger()` failing when regions argument was used.
* Fixed heatmaps not staying in a single column when more than 2 groups were present.

## Version 2.0.0
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

## Version 1.1.4
* Added palette argument to aggregate plots
* Added `exons_to_genes()` function to convert exon annotation to gene annotation
* Added `plot_granges_heatmap()` function to use GRanges for plotting heatmaps

## Version 1.1.3

* Fixed group handling for list region input in `plot_agg_regions()`
* Fixed unused window size argument in `plot_region_heatmap()`
* Fixed error when reads overlap in name and position for internal function `StatLM()`

## Version 1.1.2

* Changed example dataset exon annotations from all genomic exons to just those contained in data.
* Fixed methylation heatmap to no longer be hard coded for Peg3.
* Added `plot_region_heatmap()` as analogue to `plot_region()`.
* Fixed `plot_agg_regions_sample_grouped()` to use `group` column of `NanoMethViz::samples(x)` rather than `haplotype`.
* Added unit tests.

## Version 1.1.1

* Added methylation heatmap via `plot_gene_heatmap()`.
* Fixed `gene_anno()` in `plot_gene()` for argument so FALSE actually turns off gene annotation.
* Added warning for cpp11 versions <0.2.5 which may cause memory crashes when trying to import methylation data.
* Added cpp11 version dependency to address tidyverse/readr#1145.
* Added query methylation by gene using `query_methy_gene()`.

## Version 1.0.0

* Initial Bioconductor release.
