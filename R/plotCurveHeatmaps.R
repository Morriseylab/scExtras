#'Plot heatmap of the pseudotime data
#' @param object Seurat object
#' @param curve curve to be plotted
#' @param filename output filename
#' @param n number of genes to plot. Default is 25
#' @import dplyr tidyr Seurat
#' @export
#'
plotCurveHeatmaps <- function(object=NULL,curve=NULL,filename='heatmap.png',n=25){
  c = sym(curve)
  cells = object@meta.data %>% tibble::rownames_to_column('cellid') %>% arrange(desc(!!c)) %>% filter(!is.na(!!c))
  genes = object@misc$sds$dge[[curve]] %>% arrange(p.value) %>% head(n) %>% pull(gene)
  FetchData(object,genes,cells$cellid) %>% t(.) %>%
    NMF::aheatmap(.,Colv=NA,distfun='pearson',scale='row',annCol=cells$var_celltype,annColors = list(X1=cpallette),
                  filename=filename)

}
