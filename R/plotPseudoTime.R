#'run Pseudotime
#' @param object Seurat object
#' @param groupby variable to group by
#' @param reduction Reduction method. Default is 'dm'
#' @param dims number of dimensions
#' @import dplyr tidyr Seurat ggplot2
#' @export
#'
plotPseudoTime = function(object,groupby,reduction='dm',dims=1:2){

  curved <- bind_rows(lapply(names(object@misc$sds$data@curves), function(x){c <- slingCurves(object@misc$sds$data)[[x]]
  d <- as.data.frame(c$s[c$ord,seq_len(2)])
  d$curve<-x
  return(d)
  })
  )

  dims <- paste0(Key(object = object[[reduction]]), dims)

  p=FetchData(object = object, vars = c(dims,groupby)) %>%
    ggplot(.,aes_string(x=dims[1],y=dims[2]))+geom_point(aes(color=!!sym(groupby))) +
    theme(legend.position="top") +
    guides(col = guide_legend(nrow = 2)) +
    geom_path(aes_string(dims[1], dims[2],linetype="curve"),curved,size=1)
  p



}
