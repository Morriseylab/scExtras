#'run Pseudotime
#' @param object Seurat object
#' @param group.by variable to group by
#' @param reduction Which dimensionality reduction to use
#' @param dims Dimensions to plot, must be a two-length numeric vector specifying x- and y-dimensions
#' @import dplyr tidyr Seurat ggplot2
#' @export
#'
plotPseudoTime = function(object,group.by=NULL,reduction='dm',dims=1:2){
  object[['ident']] <- Idents(object = object)
  group.by <- group.by %||% 'ident'
  curved <- bind_rows(lapply(names(object@misc$sds$data@curves), function(x){c <- slingCurves(object@misc$sds$data)[[x]]
    d <- as.data.frame(c$s[c$ord,seq_len(2)])
    d$curve<-x
    return(d)
    })
  )

  dims <- paste0(Key(object = object[[reduction]]), dims)

  p=FetchData(object = object, vars = c(dims,group.by)) %>%
    ggplot(.,aes_string(x=dims[1],y=dims[2]))+geom_point(aes(color=!!sym(group.by))) +
    theme(legend.position="top") +
    guides(col = guide_legend(nrow = 2)) +
    geom_path(aes_string(dims[1], dims[2],linetype="curve"),curved,size=1)
  p
}
