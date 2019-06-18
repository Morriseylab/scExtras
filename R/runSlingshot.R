#'run Slingshot
#' @param object Seurat object
#' @param group,by groups
#' @param reduction Reduction method. Default is 'dm'
#' @param start.clus start cluster
#' @param end.clus end cluster
#' @import dplyr tidyr Seurat slingshot
#' @export
#'
runSlingshot  <- function(object,reduction='dm',group.by=NULL, start.clus=NULL,end.clus=NULL){
  rd <- Embeddings(object,reduction)
  cl <- Idents(object = object)
  group.by <- group.by %||% 'ident'
  cl <- FetchData(object = object, vars = group.by) %>% pull(`group.by`)

  object@misc[['sds']] <-  list("dr"=reduction,"data"=slingshot(rd,cl,start.clus=start.clus,end.clus=end.clus))
  ps <- slingPseudotime(object@misc[['sds']]$data)
  object@meta.data[,colnames(ps)] <- ps
  return(object)
}
