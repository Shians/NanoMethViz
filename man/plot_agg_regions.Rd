% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_agg_regions.R
\name{plot_agg_regions}
\alias{plot_agg_regions}
\title{Plot aggregate regions}
\usage{
plot_agg_regions(
  x,
  regions,
  group_col = NULL,
  flank = 2000,
  stranded = TRUE,
  span = 0.05,
  palette = ggplot2::scale_colour_brewer(palette = "Set1")
)
}
\arguments{
\item{x}{the NanoMethResult object.}

\item{regions}{a table of regions containing at least columns chr, strand,
start and end. Any additiona columns can be used for grouping.}

\item{group_col}{the column to group aggregated trends by. This column can
be in from the regions table or samples(x).}

\item{flank}{the number of flanking bases to add to each side of each region.}

\item{stranded}{TRUE if negative strand features should have coordinates
flipped to reflect features like transcription start sites.}

\item{span}{the span for loess smoothing.}

\item{palette}{the colour palette used for groups.}
}
\value{
a ggplot object containing the aggregate methylation trend.
}
\description{
Plot aggregate regions
}
\examples{
nmr <- load_example_nanomethresult()
gene_anno <- exons_to_genes(NanoMethViz::exons(nmr))
plot_agg_regions(nmr, gene_anno)
plot_agg_regions(nmr, gene_anno, group_col = "sample")
plot_agg_regions(nmr, gene_anno, group_col = "group")

}
