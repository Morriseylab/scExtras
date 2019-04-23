#'Create a 3D plot
#' @param object Seurat object
#' @param groupby grouping variable to color by. Has to be a column name present in metadata of the seurat object
#' @param reduction Reduction method. Default is 'dm'
#' @param colors color palette
#' @import dplyr tidyr Seurat
#' @export
#'
make3dPlot <- function(object,groupby,reduction='dm',colors=NULL){
  dims=1:3
  dims <- paste0(Key(object = object[[reduction]]), dims)
  data <- FetchData(object = object, vars = c(dims,groupby))

  if(is.factor(data[,groupby])){
    colors=cpallette
  }
  plot_ly(data, x=~get(dims[1]), y=~get(dims[2]), z=~get(dims[3]),colors=colors,color=~get(groupby),size=.5 ) %>%
    add_markers()
}
