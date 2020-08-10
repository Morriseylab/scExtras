#'Function to run doublet finder
#' @param object Seurat object
#' @param dims Number of dimensions to use
#' @param group variable to group the cells by
#' @param doublet.formrate Doublet formation rate (Refer to Cellranger documentation. Depends on number of cells fed in)
#' @param sct logical (T/F)
#' @import DoubletFinder
#' @export
RunDoubletfinder <- function(object, dims,group="var_cluster",doublet.formrate =0.75, sct=T){

  ### pK Identification (no ground-truth)
  sweep.res <- paramSweep_v3(object, PCs = dims, sct = sct)
  sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  pK=as.numeric(as.character(bcmvn$pK))
  BCmetric=bcmvn$BCmetric
  pK_choose = pK[which(BCmetric %in% max(BCmetric))]

  ## Homotypic Doublet Proportion Estimate
  annotations=eval(parse(text= paste("object@meta.data$",group,sep="")))
  homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
  nExp_poi <- round(doublet.formrate*length(rownames(object)))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

  ## Run DoubletFinder
  object <- doubletFinder_v3(object, PCs = dims, pN = 0.25, pK = pK_choose, nExp = nExp_poi, reuse.pANN = FALSE, sct = sct)
  colnames(object@meta.data)[ncol(object@meta.data)]="var_doubletfinder_res"

  return(object)
}
