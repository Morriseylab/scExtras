#'Run PCA and Jackstaw and plots.
#' @param object Seurat object
#' @param npcs Compute N PC dims
#' @param jackstraw Run jackstraw
#' @param plotdir Where to export plots.
#' @export
PCATools <- function(object,npcs=50,jackstraw=T,plotdir='./'){
  object <- RunPCA(object = object, npcs = npcs, verbose = FALSE)
  if(jackstraw==TRUE){
    object <- JackStraw(object = object, num.replicate = 100,dims = npcs)
    object <- ScoreJackStraw(object = object,dims=1:npcs)
  }
  object <- LogSeuratCommand(object = object)
  object
}

#'Run all Dimension Reduction methods and find marker genes between clusters
#' @param object Seurat object
#' @param k n number of neighbors for umap
#' @param DM Run diffusion map
#' @param UMAP Run UMAP
#' @param TSNE Run TSNE
#' @param finfallnarkers T/F run findallmarkers
#' @param resolution Resolution param for FindCluster
#' @import dplyr tidyr Seurat
#' @export
ClusterDR <-function(object,k=30, dims=dim,DM=F,UMAP=T,TSNE=T,findallamrkers=T,resolution=0.5,n.componets=2){

  if(TSNE==TRUE){
     object <- RunTSNE(object = object, reduction = "pca",dims = dims)
  }
  if(UMAP==TRUE){
    object <- RunUMAP(object = object, reduction = "pca", n.neighbors = k,n.components = n.components,dims = dims)
  }
  if(DM==TRUE){
    object <- RunDiffusion(object = object,dims=dims)
  }
  object <- FindNeighbors(object = object,dims=dims,k.param = k)
  object <- FindClusters(object = object,resolution=resolution)
  object$var_cluster <- object@active.ident
  object@misc[["findallmarkers"]] <- FindAllMarkers(object = object, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  object
}
