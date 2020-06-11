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

  #Set ident to groupby option
  if(!is.null(group.by)){
      Idents(object)<- group.by
  }

  cl <- Idents(object = object)
  ### Need to fix this not working
  #group.by <- group.by %||% 'ident'
  #cl <- FetchData(object = object, vars = group.by) %>% pull(`group.by`)

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

#'CurvePlot Plot Slingshot Curves
#' @param object Seurat object
#' @param sds Slingshot Data object
#' @param group.by variable to group by
#' @param reduction Which dimensionality reduction to use, default UMAP
#' @param dims Dimensions to plot, must be a two-length numeric vector specifying x- and y-dimensions
#' @param cols Color palette
#' @param label Label plots
#' @import dplyr tidyr Seurat ggplot2
#' @export
#'
CurvePlot = function(object,
                          sds=NULL,
                          group.by = NULL,
                          reduction = 'umap',
                          dims = 1:2,
                          cols=NULL,
                          label=T
                          ) {
  object[['ident']] <- Idents(object = object)
  group.by <- group.by %||% 'ident'
  dims <- paste0(Key(object = object[[reduction]]), dims)

  curved <-
    bind_rows(lapply(names(slingCurves(sds)), function(x) {
      c <- slingCurves(sds)[[x]]
      d <- as.data.frame(c$s[c$ord, dims])
      d$curve <- x
      return(d)
    }))


  DimPlot(object,cols=cols,label = label,group.by = group.by,reduction = reduction) +
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
#' @param sdsname Name of Slingshot object stored in Seurat Object
#' @param features Vector of genes to be plotted
#' @param lineage THe linage to be plotted such as lineage1
#' @param col Color palette, this vector needs to have the names be the cell types or cluster names
#' @import dplyr tidyr Seurat ComplexHeatmap
#' @export
#'
plotlineageHeatMap <- function(object,sdsname,features,lineage='lineage1',col){

  ## Maybe add curve is user puts in integer

  sds <- object@misc[[sdsname]]$data

  ### Should add a check for lineages in model

  qlineage <- quo(lineage)

  cells <- inner_join(
    slingCurveWeights(sds,as.probs=T) %>% as.data.frame() %>%
      setNames(gsub('curve','lineage',names(.))) %>%
      rownames_to_column('cellid') %>%
      gather(curve,w,-cellid) %>%
      group_by(cellid) %>%
      top_n(n=1,wt=w) %>%
      filter(curve==!!lineage),
    slingPseudotime(sub@misc$umap_cl3$data) %>% as.data.frame %>%
      setNames(gsub('curve','lineage',names(.))) %>%
      rownames_to_column('cellid') %>%
      select(cellid,!!qlineage) %>%
      dplyr::rename('time'=lineage)
  ) %>% arrange(time) %>%
    inner_join(., object@meta.data %>% rownames_to_column('cellid') )

  data <- FetchData(object=object, vars=features,cells=cells$cellid) %>% t(.)
  mat_scaled = t(scale(t(data)))

  f1=circlize::colorRamp2(c(-2,0,2), c('skyblue1', "grey10","yellow"))


  col_fun =circlize::colorRamp2(c(0, 20), c("blue", "red"))

  ha = HeatmapAnnotation(
    pseudotime=anno_lines(cells$time),
    celltype = cells$var_celltype,
    col = list(celltype = cellcolorpal
    ))

  ht1 <- Heatmap(mat_scaled,
                 col=f1,
                 show_row_dend = F,
                 row_names_side='left',
                 show_column_names = F,
                 cluster_columns = F,
                 cluster_rows = F,
                 top_annotation = ha
  )
  ht1



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

