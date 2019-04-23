#'Get maximum number of PC dimensions
#' @param object Seurat object
#' @import dplyr tidyr Seurat
#' @export
#'
getMaxDim <- function(object){

  object@reductions$pca@jackstraw$overall.p.values %>%
    as.data.frame(.) %>%
    mutate(adj = p.adjust(Score,method='bonferroni')) %>%
    filter(adj <0.05) %>%
    summarise(max=max(PC)) %>%
    pull(max)


}
