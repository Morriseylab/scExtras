% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Plotting.R
\name{BiGenePlot}
\alias{BiGenePlot}
\title{bi_getValues utility function.}
\usage{
BiGenePlot(
  object,
  feature1,
  feature2,
  feature1.min = 1,
  feature2.min = 1,
  dims = 1:2,
  reduction = "umap",
  pt.size = 0.1,
  cols = c("#E41A1C", "#377EB8", "#4DAF4A", "grey75"),
  plotLineage = FALSE,
  title = NULL
)
}
\arguments{
\item{object}{Seurat object}

\item{feature1}{Feature 1}

\item{feature2}{Feature 2}

\item{feature1.min}{Min UMI for a cell to be postive for feature. Set higher for more a restrictive selection}

\item{feature2.min}{Min UMI for a cell to be postive for feature. Set higher for more a restrictive selection}

\item{dims}{Dimensions to plot, must be a two-length numeric vector specifying x- and y-dimensions}

\item{reduction}{Which dimensionality reduction to use}

\item{pt.size}{Adjust point size}

\item{cols}{Vector of colors for Feature1+/Feature2+, feature1+,feature2+ and feature1-/feature2-}

\item{title}{Plot Title}
}
\description{
bi_getValues utility function.
}
