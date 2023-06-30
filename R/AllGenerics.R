#' Get methylation data
#' @keywords internal
#' @param object the object.
#' @return the path to the methylation data.
#' @examples
#' showMethods("methy")
#'
#' @export
setGeneric("methy", valueClass = "ANY", function(object) {
    standardGeneric("methy")
})

#' Set methylation data
#' @keywords internal
#' @export
setGeneric("methy<-", function(object, value) {
    standardGeneric("methy<-")
})

#' Get sample annotation
#' @keywords internal
#' @export
setGeneric("samples", valueClass = "data.frame", function(object) {
    standardGeneric("samples")
})

#' Set sample annotation
#' @keywords internal
#' @export
setGeneric("samples<-", function(object, value) {
    standardGeneric("samples<-")
})

#' Get exon annotation
#' @keywords internal
#' @export
setGeneric("exons", valueClass = "data.frame", function(object) {
    standardGeneric("exons")
})

#' Set exon annotation
#' @keywords internal
#' @export
setGeneric("exons<-", function(object, value) {
    standardGeneric("exons<-")
})

#' Plot gene methylation
#'
#' Plot the methylation of a gene symbol specified within the exon(x) slot.
#'
#' @param x the NanoMethResult or ModBamResult object.
#' @param gene the gene symbol for the gene to plot.
#'
#' @return a patchwork plot containing the methylation profile in the specified
#'   region.
#'
#' @importFrom ggrastr rasterise
#' @export
setGeneric("plot_gene", function(x, gene, ...) {
    standardGeneric("plot_gene")
})

#' Plot gene methylation heatmap
#'
#' Plot the methylation heatmap of a gene symbol specified within the exon(x) slot.
#'
#' @param x the NanoMethResult or ModBamResult object.
#' @param gene the gene symbol for the gene to plot.
#'
#' @return a ggplot object of the heatmap
#'
#' @export
setGeneric("plot_gene_heatmap", function(x, gene, ...) {
    standardGeneric("plot_gene_heatmap")
})

#' Plot region methylation
#'
#' Plot the methylation of a genomic region.
#'
#' @param x the NanoMethResult or ModBamResult object.
#' @param chr the chromosome to plot.
#' @param start the start of the plotting region.
#' @param end the end of the plotting region.
#'
#' @return a patchwork plot containing the methylation profile in the specified
#'   region.
#'
#' @importFrom ggrastr rasterise
#' @export
setGeneric("plot_region", function(x, chr, start, end, ...) {
    standardGeneric("plot_region")
})

#' Plot region methylation heatmap
#'
#' Plot the methylation heatmap of a genomic region.
#'
#' @param x the NanoMethResult or ModBamResult object.
#' @param chr the chromosome to plot.
#' @param start the start of the plotting region.
#' @param end the end of the plotting region.
#'
#' @return a ggplot object of the heatmap.
#'
#' @export
setGeneric("plot_region_heatmap", function(x, chr, start, end, ...) {
    standardGeneric("plot_region_heatmap")
})
