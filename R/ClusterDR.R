#'Run all Dimension Reduction methods and find marker genes between clusters
#' @param object Seurat object
#' @param npcs number of PC's to run
#' @param maxdim maximum number of dimensions. Defaults to npcs specified
#' @param k n number of neighbors for umap
#' @import dplyr tidyr Seurat
#' @export

ClusterDR <-function(object,npcs=50, maxdim='auto',k=30, DM=F,UMAP=T,TSNE=T,resolution=0.5,){
  object <- RunPCA(object = object, npcs = npcs, verbose = FALSE)

  if(maxdim=='auto'){
    object <- JackStraw(object = object, num.replicate = 100,dims = npcs)
    object <- ScoreJackStraw(object = object,dims=1:npcs)
    dim <- object@reductions$pca@jackstraw$overall.p.values %>%
      as.data.frame(.) %>%
      mutate(adj = p.adjust(Score,method='bonferroni')) %>%
      filter(adj <0.05) %>%
      summarise(max=max(PC)) %>%
      pull(max)

  } else {
    dim<-maxdim
  }
  print(dim)
  if(TSNE==TRUE){
     object <- RunTSNE(object = object, reduction = "pca",dims = 1:dim)
  }
  if(UMAP==TRUE){
    object <- RunUMAP(object = object, reduction = "pca", n.neighbors = k,n.components = 3,dims = 1:dim)
  }
  if(DM==TRUE){
    object <- RunDiffusion(object = object,dims=1:dim)
  }
  object <- FindNeighbors(object = object,dims=1:dim,k.param = k)
  object <- FindClusters(object = object,resolution=resolution)
  object$var_cluster <- object@active.ident
  object@misc[["findallmarkers"]] <- FindAllMarkers(object = object, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  object
}
