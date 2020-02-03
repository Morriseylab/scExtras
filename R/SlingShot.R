#'run Slingshot
#' @param object Seurat object
#' @param sds.name Name in hte Misc slot to save the Slingshot object
#' @param group,by groups
#' @param reduction Reduction method. Default is 'dm'
#' @param start.clus start cluster
#' @param end.clus end cluster
#' @param approx_points
#' @param allow.breaks logical, determines whether curves that branch very close to the origin should be allowed to have different starting points.
#' @param extend acter, how to handle root and leaf clusters of lineages when constructing the initial, piece-wise linear curve. Accepted values are 'y' (default), 'n', and 'pc1'. See 'Details' for more.
#' @param stretch numeric factor by which curves can be extrapolated beyond endpoints. Default is 2, see principal_curve.
#' @import dplyr tidyr Seurat slingshot
#' @export
#'
runSlingshot  <- function(object,sds.name='sds',reduction='dm',group.by=NULL, start.clus=NULL,end.clus=NULL, approx_points = FALSE, allow.breaks=TRUE, extend='n',stretch=0){
  rd <- Embeddings(object,reduction)
  cl <- Idents(object = object)
  group.by <- group.by %||% 'ident'
  cl <- FetchData(object = object, vars = group.by) %>% pull(`group.by`)

  object@misc[[sds.name]] <-  list("dr"=reduction,"data"=slingshot(rd,cl,start.clus=start.clus,end.clus=end.clus,approx_points=approx_points,allow.breaks=allow.breaks,extend=extend,stretch=stretch))
  #ps <- slingPseudotime(object@misc[['sds']]$data)
  #object@meta.data[,colnames(ps)] <- as.data.frame(ps)
  object <- LogSeuratCommand(object = object)
  return(object)
}

#'run Psuedotime Diff Expression
#' @param object Seurat object
#' @import dplyr tidyr Seurat gam
#' @export
#'
runPseudoTimeDGE <- function(object){
  var <- VariableFeatures(scrna)

  DGE <- list()
  for(c in names(object@misc$sds$data@curves)){
    object@misc$sds$dge[[c]] <- FetchData(object,append(var, c,0)) %>% tidyr::gather(gene,signal, -one_of(c)) %>% dplyr::rename(curve = 1) %>%
      tidyr::nest(-gene) %>%
      mutate(
        fit = purrr::map(data, ~ gam::gam(signal ~ lo(curve), data = .x)),
        tidied = purrr::map(fit, tidy)
      ) %>%
      tidyr::unnest(tidied) %>%
      filter(term !='Residuals')
  }
  object <- LogSeuratCommand(object = object)
  return(object)

}

#'run Pseudotime
#' @param object Seurat object
#' @param group.by variable to group by
#' @param reduction Which dimensionality reduction to use
#' @param dims Dimensions to plot, must be a two-length numeric vector specifying x- and y-dimensions
#' @import dplyr tidyr Seurat ggplot2
#' @export
#'
plotPseudoTime = function(object,
                          sds=NULL,
                          group.by = NULL,
                          reduction = 'dm',
                          dims = 1:2
                          ) {
  object[['ident']] <- Idents(object = object)
  group.by <- group.by %||% 'ident'
  dims <- paste0(Key(object = object[[reduction]]), dims)

  curved <-
    bind_rows(lapply(names(getCurves(sds)), function(x) {
      c <- slingCurves(sds)[[x]]
      d <- as.data.frame(c$s[c$ord, dims])
      d$curve <- x
      return(d)
    }))


  DimPlot(sub,cols=cpallette,label = T,group.by = group.by) +
    geom_path(aes_string(dims[1], dims[2], linetype = "curve"), curved, size =1)




  # data <- FetchData(object = object, vars = c(dims, group.by))
  # p <- ggplot(data, aes_string(x = dims[1], y = dims[2])) +
  #   geom_point(aes(color =!!sym(group.by))) +
  #   theme(legend.position = "top")
  #
  # if (is.character(data[[group.by]]) | is.factor(data[[group.by]])){
  #   p <- p + guides(col = guide_legend(nrow = 2))
  # } else {
  #   p <- p + scale_color_distiller(palette = "RdYlBu", na.value = 'grey90')
  # }
  # p <- p +geom_path(aes_string(dims[1], dims[2], linetype = "curve"), curved, size =1)
  # p


}


#'Plot heatmap of the pseudotime data
#' @param object Seurat object
#' @param curve curve to be plotted
#' @param filename output filename
#' @param n number of genes to plot. Default is 25
#' @import dplyr tidyr Seurat
#' @export
#'
plotCurveHeatmaps <-
  function(object = NULL,
           curve = NULL,
           annCol= 'var_cluster',
           filename = 'heatmap.png',
           n = 25) {
    cells <- object@meta.data %>% tibble::rownames_to_column('cellid')  %>% dplyr::arrange(!!sym(curve)) %>% filter(!is.na(!!sym(curve)))
    genes <- object@misc$sds$dge[[curve]] %>% dplyr::arrange(p.value) %>% head(n) %>% pull(gene)
    FetchData(object=object, vars=genes, cells=cells$cellid) %>% t(.) %>%
      NMF::aheatmap(
        .,
        Colv = NA,
        distfun = 'pearson',
        scale = 'row',
        annCol = cells[[annCol]],
        annColors = list(X1 = cpallette),
        filename = filename
      )

  }


#'Create feature plots of genes
#' @param object Seurat object
#' @param curve curve to be plotted
#' @param reduction Reduction method. Default is 'dm'
#' @param n number of genes to plot. Default is 25
#' @import dplyr tidyr Seurat
#' @export
#'
plotCurveDGEgenes <- function(object=NULL,curve=NULL,n=25,reduction='dm'){
  genes = object@misc$sds$dge[[curve]] %>% arrange(p.value) %>% head(n) %>% pull(gene)
  plot_grid(  plotlist = FeaturePlot(scrna,genes,reduction = reduction,cols = c('grey','purple')))

}

