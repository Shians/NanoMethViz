% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/query_exons.R
\name{query_exons}
\alias{query_exons}
\alias{query_exons_region}
\alias{query_exons_gene_id}
\alias{query_exons_symbol}
\title{Query exons}
\usage{
query_exons_region(x, chr, start, end)

query_exons_gene_id(x, gene_id)

query_exons_symbol(x, symbol)
}
\arguments{
\item{x}{the object to query.}

\item{chr}{the chromosome to query.}

\item{start}{the start of the query region.}

\item{end}{the end of the query region.}

\item{gene_id}{the gene_id to query.}

\item{symbol}{the gene_id to query.}
}
\value{
data.frame of queried exons.
}
\description{
Query a data.frame, NanoMethResult or ModBamResult for exon annotation.
}
\section{Functions}{
\itemize{
\item \code{query_exons_region()}: Query region.

\item \code{query_exons_gene_id()}: Query gene ID.

\item \code{query_exons_symbol()}: Query gene symbol.

}}
