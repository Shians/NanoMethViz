% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/convert_bsseq_to_edger.R
\name{bsseq_to_edger}
\alias{bsseq_to_edger}
\title{Convert BSseq object to edgeR methylation matrix}
\usage{
bsseq_to_edger(bsseq)
}
\arguments{
\item{bsseq}{the BSseq object.}
}
\value{
a matrix compatible with the edgeR differential methylation pipeline
}
\description{
Convert BSseq object to edgeR methylation matrix
}
\examples{
methy <- system.file("methy_subset.tsv.bgz", package = "NanoMethViz")
bsseq <- methy_to_bsseq(methy)
edger_mat <- bsseq_to_edger(bsseq)
}