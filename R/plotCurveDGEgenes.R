#'Create feature plots of genes
#' @param object Seurat object
#' @param curve curve to be plotted
#' @param reduction.use Reduction method. Default is 'dm'
#' @param n number of genes to plot. Default is 25
#' @import dplyr tidyr Seurat
#' @export
#'
plotCurveDGEgenes <- function(object=NULL,curve=NULL,n=25,reduction.use='dm'){
  genes = object@misc$sds$dge[[curve]] %>% arrange(p.value) %>% head(n) %>% pull(gene)
  plot_grid(  plotlist = FeaturePlot(scrna,genes,reduction = reduction.use,cols = c('grey','purple')))

}
