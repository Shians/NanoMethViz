#' Get methylation data
#' @keywords internal
#'
#' @param object the object.
#'
#' @return the path to the methylation data.
#'
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
#'
#' @param object the object.
#'
#' @return the sample annotation.
#'
#' @examples
#' showMethods("samples")
#'
#' @export
#'
#' @keywords internal
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
#'
#' @param object the object.
#'
#' @return the exon annotation.
#'
#' @examples
#' showMethods("exons")
#'
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

#' Plot gene
#'
#' @param x the NanoMethResult object.
#' @param gene the gene symbol for the gene to plot.
#' @param ... additional arguments
#'
#' @return a patchwork plot containing the methylation profile in the specified
#'   region.
#'
#' @examples
#' nmr <- load_example_nanomethresult()
#' plot_gene(nmr, "Peg3")
#'
#' @importFrom ggrastr rasterise
#' @export
setGeneric("plot_gene", function(x, gene, ...) {
    standardGeneric("plot_gene")
})

#' Plot gene methylation heatmap
#'
#' @param x the NanoMethResult object.
#' @param gene the gene symbol for the gene to plot.
#' @param ... additional arguments
#'
#' @return a ggplot object of the heatmap
#'
#' @examples
#' nmr <- load_example_nanomethresult()
#' plot_gene_heatmap(nmr, "Peg3")
#'
#' @export
setGeneric("plot_gene_heatmap", function(x, gene, ...) {
    standardGeneric("plot_gene_heatmap")
})

#' Plot region
#'
#' @param x the NanoMethResult object.
#' @param chr the chromosome to plot.
#' @param start the start of the plotting region.
#' @param end the end of the plotting region.
#' @param ... additional arguments.
#'
#' @return a patchwork plot containing the methylation profile in the specified
#'   region.
#'
#' @examples
#' nmr <- load_example_nanomethresult()
#' plot_region(nmr, "chr7", 6703892, 6730431)
#'
#' @importFrom ggrastr rasterise
#' @export
setGeneric("plot_region", function(x, chr, start, end, ...) {
    standardGeneric("plot_region")
})

#' Plot region methylation heatmap
#'
#' @param x the NanoMethResult object.
#' @param chr the chromosome to plot.
#' @param start the start of the plotting region.
#' @param end the end of the plotting region.
#' @param ... additional arguments.
#'
#' @return a ggplot object of the heatmap.
#'
#' @examples
#' nmr <- load_example_nanomethresult()
#' plot_region_heatmap(nmr, "chr7", 6703892, 6730431)
#'
#' @export
setGeneric("plot_region_heatmap", function(x, chr, start, end, ...) {
    standardGeneric("plot_region_heatmap")
})
