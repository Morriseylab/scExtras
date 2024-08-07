% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ligrec.R
\name{RunLigRec}
\alias{RunLigRec}
\title{Get ligand receptor pairs between all clusters}
\usage{
RunLigRec(object, group.by = "ident", org, perc = 30, assay = "RNA")
}
\arguments{
\item{object}{Seurat object}

\item{group.by}{grouping variable}

\item{org}{organism}

\item{perc}{minimum percentage of cells that should express the ligand or receptor}
}
\value{
Returns the object with object@misc[[ligrecres]] containing the matrix of lignad-receptor pair data
}
\description{
Get ligand receptor pairs between all clusters
}
