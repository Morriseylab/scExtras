% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SlingShot.R
\name{plotLineageHeatMap}
\alias{plotLineageHeatMap}
\title{Plot heatmap of the pseudotime data}
\usage{
plotLineageHeatMap(
  object,
  features,
  lineage = "lineage1",
  col,
  group.by = NULL,
  features.callout = NULL,
  show_row_names = TRUE,
  show_heatmap_legend = TRUE
)
}
\arguments{
\item{object}{Seurat object}

\item{features}{Vector of genes to be plotted}

\item{lineage}{THe linage to be plotted such as lineage1}

\item{col}{Color palette, this vector needs to have the names be the cell types or cluster names}

\item{group.by}{Cluster variable used in slingshot}

\item{features.callout}{vector of gene names to be called out}

\item{show_row_names}{T/F to show row names}

\item{show_heatmap_legend}{T/F show legend}
}
\description{
Plot heatmap of the pseudotime data
}
