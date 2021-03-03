### Version 1.1.1

* Added methylation heatmap via `plot_gene_heatmap()`
* Fixed `gene_anno()` in `plot_gene()` for argument so FALSE actually turns off gene annotation
* Added warning for cpp11 versions <0.2.5 which may cause memory crashes when trying to import methylation data
* Added cpp11 version dependency to address tidyverse/readr#1145
* Added query methylation by gene using `query_methy_gene()`

### Version 1.0.0

* Initial Bioconductor release
