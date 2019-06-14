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

#'Plot Pairwise PCA plots.
#' @param object Seurat object
#' @param npcs Compute N PC dims
#' @export
PCAPwPlot <- function(object,npcs=50)
  plist <- list()
for (i in seq(1,npcs,by=2)){
  plist[[as.character(i)]] <- FeatureScatter(scrna,paste0('PC_',i),paste0('PC_',i+1)) + geom_point(color='gray') +
    theme(legend.position="none",
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank()
    ) +
    ggtitle(paste0('PC',i,' vs PC',i+1))
}
  plot_grid(plotlist = plist,nrow = 5)
}


#'Run all Dimension Reduction methods and find marker genes between clusters
#' @param object Seurat object
#' @param k n number of neighbors for umap
#' @param DM Run diffusion map
#' @param UMAP Run UMAP
#' @param TSNE Run TSNE
#' @param resolution Resolution param for FindCluster
#' @import dplyr tidyr Seurat
#' @export
ClusterDR <-function(object,k=30, dims=dim,DM=F,UMAP=T,TSNE=T,resolution=0.5){

  if(TSNE==TRUE){
     object <- RunTSNE(object = object, reduction = "pca",dims = dims)
  }
  if(UMAP==TRUE){
    object <- RunUMAP(object = object, reduction = "pca", n.neighbors = k,n.components = 3,dims = dims)
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
