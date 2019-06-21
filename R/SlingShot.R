#'run Slingshot
#' @param object Seurat object
#' @param group,by groups
#' @param reduction Reduction method. Default is 'dm'
#' @param start.clus start cluster
#' @param end.clus end cluster
#' @import dplyr tidyr Seurat slingshot
#' @export
#'
runSlingshot  <- function(object,reduction='dm',group.by=NULL, start.clus=NULL,end.clus=NULL){
  rd <- Embeddings(object,reduction)
  cl <- Idents(object = object)
  group.by <- group.by %||% 'ident'
  cl <- FetchData(object = object, vars = group.by) %>% pull(`group.by`)

  object@misc[['sds']] <-  list("dr"=reduction,"data"=slingshot(rd,cl,start.clus=start.clus,end.clus=end.clus))
  ps <- slingPseudotime(object@misc[['sds']]$data)
  object@meta.data[,colnames(ps)] <- ps
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

