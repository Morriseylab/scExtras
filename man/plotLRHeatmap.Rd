% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotLRHeatmap.R
\name{plotLRHeatmap}
\alias{plotLRHeatmap}
\title{Create ligand receptor heatmap plot}
\usage{
plotLRHeatmap(result, clusterby)
}
\arguments{
\item{result}{dataframe containing the fields- interacting_pair,Receptor,Ligand,Receptor_cluster,Ligand_cluster}

\item{clusterby}{clusterby option for heatmap. Choose from "both", "row","column" and "none"}
}
\value{
Returns heatmap with ligand-receptor pair data
}
\description{
Create ligand receptor heatmap plot
}
