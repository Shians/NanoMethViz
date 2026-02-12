# Package index

## Package

- [`NanoMethViz`](https://shians.github.io/NanoMethViz/reference/NanoMethViz-package.md)
  [`NanoMethViz-package`](https://shians.github.io/NanoMethViz/reference/NanoMethViz-package.md)
  : NanoMethViz: Visualise methylation data from Oxford Nanopore
  sequencing

## Plotting Regions

- [`plot_gene()`](https://shians.github.io/NanoMethViz/reference/plot_gene.md)
  : Plot gene methylation
- [`plot_gene_heatmap()`](https://shians.github.io/NanoMethViz/reference/plot_gene_heatmap.md)
  : Plot gene methylation heatmap
- [`plot_region()`](https://shians.github.io/NanoMethViz/reference/plot_region.md)
  : Plot region methylation
- [`plot_region_heatmap()`](https://shians.github.io/NanoMethViz/reference/plot_region_heatmap.md)
  : Plot region methylation heatmap
- [`plot_grange()`](https://shians.github.io/NanoMethViz/reference/plot_grange.md)
  : Plot GRanges
- [`plot_grange_heatmap()`](https://shians.github.io/NanoMethViz/reference/plot_grange_heatmap.md)
  : Plot GRanges heatmap

## Other Plots

- [`plot_agg_genes()`](https://shians.github.io/NanoMethViz/reference/plot_agg_genes.md)
  : Plot gene aggregate plot
- [`plot_agg_regions()`](https://shians.github.io/NanoMethViz/reference/plot_agg_regions.md)
  : Plot aggregate regions
- [`plot_mds()`](https://shians.github.io/NanoMethViz/reference/plot_mds.md)
  : Plot MDS
- [`plot_pca()`](https://shians.github.io/NanoMethViz/reference/plot_pca.md)
  : Plot PCA
- [`plot_violin()`](https://shians.github.io/NanoMethViz/reference/plot_violin.md)
  : Plot violin for regions

## Feature Annotations

- [`get_cgi_mm10()`](https://shians.github.io/NanoMethViz/reference/get_exons.md)
  [`get_cgi_grcm39()`](https://shians.github.io/NanoMethViz/reference/get_exons.md)
  [`get_cgi_t2t()`](https://shians.github.io/NanoMethViz/reference/get_exons.md)
  [`get_cgi_hg19()`](https://shians.github.io/NanoMethViz/reference/get_exons.md)
  [`get_cgi_hg38()`](https://shians.github.io/NanoMethViz/reference/get_exons.md)
  [`get_exons_mm10()`](https://shians.github.io/NanoMethViz/reference/get_exons.md)
  [`get_exons_grcm39()`](https://shians.github.io/NanoMethViz/reference/get_exons.md)
  [`get_exons_hg19()`](https://shians.github.io/NanoMethViz/reference/get_exons.md)
  [`get_exons_hg38()`](https://shians.github.io/NanoMethViz/reference/get_exons.md)
  [`get_exons_t2t()`](https://shians.github.io/NanoMethViz/reference/get_exons.md)
  : Get exon annotations
- [`get_exons_homo_sapiens()`](https://shians.github.io/NanoMethViz/reference/get_exons_homo_sapiens.md)
  : Get exon annotations for Homo sapiens (hg19)
- [`get_exons_mus_musculus()`](https://shians.github.io/NanoMethViz/reference/get_exons_mus_musculus.md)
  : Get exon annotations for Mus musculus (mm10)
- [`exons_to_genes()`](https://shians.github.io/NanoMethViz/reference/exons_to_genes.md)
  : Convert exon annotation to genes

## Conversion

- [`bsseq_to_edger()`](https://shians.github.io/NanoMethViz/reference/bsseq_to_edger.md)
  : Convert BSseq object to edgeR methylation matrix
- [`bsseq_to_log_methy_ratio()`](https://shians.github.io/NanoMethViz/reference/bsseq_to_log_methy_ratio.md)
  : Convert BSseq object to log-methylation-ratio matrix
- [`methy_to_bsseq()`](https://shians.github.io/NanoMethViz/reference/methy_to_bsseq.md)
  : Create BSSeq object from methylation tabix file
- [`methy_to_edger()`](https://shians.github.io/NanoMethViz/reference/methy_to_edger.md)
  : Convert NanoMethResult object to edgeR methylation matrix
- [`modbam_to_tabix()`](https://shians.github.io/NanoMethViz/reference/modbam_to_tabix.md)
  : Convert BAM with modifications to tabix format

## Querying Methylation Data

- [`query_methy()`](https://shians.github.io/NanoMethViz/reference/query_methy.md)
  : Query methylation data
- [`filter_methy()`](https://shians.github.io/NanoMethViz/reference/filter_methy.md)
  : Create filtered methylation file

## Data Objects

- [`NanoMethResult()`](https://shians.github.io/NanoMethViz/reference/NanoMethResult-class.md)
  [`methy(`*`<NanoMethResult>`*`)`](https://shians.github.io/NanoMethViz/reference/NanoMethResult-class.md)
  [`` `methy<-`( ``*`<NanoMethResult>`*`,`*`<ANY>`*`)`](https://shians.github.io/NanoMethViz/reference/NanoMethResult-class.md)
  [`samples(`*`<NanoMethResult>`*`)`](https://shians.github.io/NanoMethViz/reference/NanoMethResult-class.md)
  [`` `samples<-`( ``*`<NanoMethResult>`*`,`*`<data.frame>`*`)`](https://shians.github.io/NanoMethViz/reference/NanoMethResult-class.md)
  [`exons(`*`<NanoMethResult>`*`)`](https://shians.github.io/NanoMethViz/reference/NanoMethResult-class.md)
  [`` `exons<-`( ``*`<NanoMethResult>`*`,`*`<data.frame>`*`)`](https://shians.github.io/NanoMethViz/reference/NanoMethResult-class.md)
  : Nanopore Methylation Result
- [`methy(`*`<ModBamResult>`*`)`](https://shians.github.io/NanoMethViz/reference/ModBamResult-class.md)
  [`` `methy<-`( ``*`<ModBamResult>`*`,`*`<ModBamFiles>`*`)`](https://shians.github.io/NanoMethViz/reference/ModBamResult-class.md)
  [`samples(`*`<ModBamResult>`*`)`](https://shians.github.io/NanoMethViz/reference/ModBamResult-class.md)
  [`` `samples<-`( ``*`<ModBamResult>`*`,`*`<data.frame>`*`)`](https://shians.github.io/NanoMethViz/reference/ModBamResult-class.md)
  [`exons(`*`<ModBamResult>`*`)`](https://shians.github.io/NanoMethViz/reference/ModBamResult-class.md)
  [`` `exons<-`( ``*`<ModBamResult>`*`,`*`<data.frame>`*`)`](https://shians.github.io/NanoMethViz/reference/ModBamResult-class.md)
  [`mod_code(`*`<ModBamResult>`*`)`](https://shians.github.io/NanoMethViz/reference/ModBamResult-class.md)
  [`` `mod_code<-`( ``*`<ModBamResult>`*`,`*`<character>`*`)`](https://shians.github.io/NanoMethViz/reference/ModBamResult-class.md)
  [`ModBamResult()`](https://shians.github.io/NanoMethViz/reference/ModBamResult-class.md)
  : modBAM methylation results
- [`ModBamFiles()`](https://shians.github.io/NanoMethViz/reference/ModBamFiles.md)
  [`show(`*`<ModBamFiles>`*`)`](https://shians.github.io/NanoMethViz/reference/ModBamFiles.md)
  : Constructor for a ModBamFiles object
- [`ModBamFiles-class`](https://shians.github.io/NanoMethViz/reference/ModBamFiles-class.md)
  : ModBamFiles class
- [`methy()`](https://shians.github.io/NanoMethViz/reference/methy.md) :
  Get methylation data
- [`samples()`](https://shians.github.io/NanoMethViz/reference/samples.md)
  : Get sample annotation
- [`exons()`](https://shians.github.io/NanoMethViz/reference/exons.md) :
  Get exon annotation
- [`query_exons_region()`](https://shians.github.io/NanoMethViz/reference/query_exons.md)
  [`query_exons_gene_id()`](https://shians.github.io/NanoMethViz/reference/query_exons.md)
  [`query_exons_symbol()`](https://shians.github.io/NanoMethViz/reference/query_exons.md)
  : Query exons
- [`create_tabix_file()`](https://shians.github.io/NanoMethViz/reference/create_tabix_file.md)
  : Create a tabix file using methylation calls

## Example Data

- [`load_example_modbamresult()`](https://shians.github.io/NanoMethViz/reference/load_example_modbamresult.md)
  : Load an example ModBamResult object
- [`load_example_nanomethresult()`](https://shians.github.io/NanoMethViz/reference/load_example_nanomethresult.md)
  : Load an example NanoMethResult object
- [`get_example_exons_mus_musculus()`](https://shians.github.io/NanoMethViz/reference/get_example_exons_mus_musculus.md)
  : Get example exon annotations for mus musculus (mm10)

## Miscallaneous

- [`region_methy_stats()`](https://shians.github.io/NanoMethViz/reference/region_methy_stats.md)
  : Calculate region methylation statistics
- [`cluster_regions()`](https://shians.github.io/NanoMethViz/reference/cluster_regions.md)
  : Cluster regions by K-means
- [`cluster_reads()`](https://shians.github.io/NanoMethViz/reference/cluster_reads.md)
  : Cluster reads based on methylation
- [`methy_col_names()`](https://shians.github.io/NanoMethViz/reference/methy_col_names.md)
  : Column names for methylation data
