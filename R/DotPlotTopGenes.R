#'Plot dotplot of top n genes
#' @param object Seurat object
#' @param nfeatures number of genes to plot. Number has to be between 1 and 10.
#' @import dplyr tidyr Seurat
#' @export
DotPlotTopGenes <- function(object,nfeatures=3){
  if(nfeatures > 50){
    print('Too many features, choose < 10')
    return()
  }

  object@misc[['findallmarkers']] %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC) %>% pull(gene) %>%
    DotPlot(object = object, features = .)

}
