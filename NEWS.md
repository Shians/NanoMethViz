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
