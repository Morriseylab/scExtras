#'run Slingshot
#' @param object Seurat object
#' @param groups groups
#' @param reduction Reduction method. Default is 'dm'
#' @param start.clus start cluster
#' @param end.clus end cluster
#' @import dplyr tidyr Seurat slingshot
#' @export
#'
runSlingshot  <- function(object,reduction='dm',groups=NULL, start.clus=NULL,end.clus=NULL){
  rd <- Embeddings(object,reduction)
  cl <- Idents(object = object)
  object@misc[['sds']] <-  list("dr"=reduction,"data"=slingshot(rd,cl,start.clus=start.clus,end.clus=end.clus))
  ps <- slingPseudotime(object@misc[['sds']]$data)
  object@meta.data[,colnames(ps)] <- ps
  return(object)
}
