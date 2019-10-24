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
#' @param n.neighbors n number of neighbors for umap
#' @param dim
#' @param DM Run diffusion map
#' @param UMAP Run UMAP
#' @param TSNE Run TSNE
#' @param findallmarkers T/F run findallmarkers
#' @param resolution Resolution param for FindCluster
#' @param n.components How many components to compute for Dim reductions (UMAP,DM,etc)
#' @param min.dist minimum distance parameter for RunUMAP function. Controls how tightly to compress the umap points together
#' @import dplyr tidyr Seurat
#' @export

ClusterDR <-function(object,n.neighbors=30, dims,DM=F,UMAP=T,TSNE=T,findallmarkers=T,resolution=0.5,n.components=2,min.dist=0.3){

  if(TSNE==TRUE){
     object <- RunTSNE(object = object, reduction = "pca",dims = dims)
  }
  if(UMAP==TRUE){
    object <- RunUMAP(object = object, reduction = "pca", n.neighbors = n.neighbors,n.components = n.components,dims = dims,min.dist=min.dist)
  }
  if(DM==TRUE){
    object <- RunDiffusion(object = object,dims=dims)
  }
  object <- FindNeighbors(object = object,dims=dims,k.param = k)
  object <- FindClusters(object = object,resolution=resolution)
  object$var_cluster <- object@active.ident
  if(findallmarkers==TRUE){
      object@misc[["findallmarkers"]] <- FindAllMarkers(object = object, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  }

  object
}
