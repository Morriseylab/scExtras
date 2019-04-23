#'Get marker genes for the given cluster
#' @param object Seurat object
#' @param cluster name of cluster you want markers for
#' @import dplyr tidyr Seurat
#' @export
getClusterMarkers <- function(object,cluster=0){

  object@misc[['findallmarkers']] %>% filter(cluster==!!cluster)

}
