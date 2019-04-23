#'Plot heatmap of top n genes
#' @param object Seurat object
#' @param nfeatures number of features to plot. Default is 10
#' @import dplyr tidyr Seurat
#' @export
#'
HeatMapTopGenes <- function(object,nfeatures=10){
  if(nfeatures > 50){
    print('Too many features, choose < 50')
    return()
  }

  object@misc[['findallmarkers']] %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC) %>% pull(gene) %>%
    DoHeatmap(object = object, features = .) + NoLegend()

}
