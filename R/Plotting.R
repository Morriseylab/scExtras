
#'Plot a stacked bar graph of
#' @param object Seurat object
#' @param group.by Variable for x axis typical cluster or celltype
#' @param color.by Variable to split each group by such as library
#' @param cols Color palette for grouping variable
#' @export

sampleBarGraph <- function(object,group.by=NULL,color.by=NULL,col=NULL){
  if(is.null(group.by)){
    print("group.by cannot be null")
    exit()
  }
  if(is.null(color.by)){
    print("color.by cannot be NULL")
    exit()
  }

  p <- object@meta.data %>%
    group_by(!!sym(color.by),!!sym(group.by)) %>%
    summarise(n = n()) %>%
    mutate(pct = n / sum(n))  %>%
    ggplot(.,aes(x=!!sym(group.by),y=pct,fill=!!sym(color.by))) +
    geom_bar(position="fill", stat="identity")  +
    scale_y_continuous(labels = scales::percent_format()) +
    theme_bw() + theme(legend.position = 'bottom')

  if(!is.null(col)){
    p <- p + scale_fill_manual(values=col)
  }
  p


}



#'Plot Pairwise PCA plots.
#' @param object Seurat object
#' @param dims Compute N PC dims
#' @export
PCAPwPlot <- function(object,dims=1:50){
  plist <- list()

  for (i in seq(min(dims),max(dims),by=2)){
    plist[[as.character(i)]] <- FeatureScatter(object,paste0('PC_',i),paste0('PC_',i+1)) + geom_point(color='gray') +
      theme(legend.position="none",
            axis.text.x = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.x = element_blank(),
            axis.ticks.y = element_blank(),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            plot.title = element_text(size=8)
      ) +
      ggtitle(paste0('PC',i,' vs PC',i+1))
  }
  plot_grid(plotlist = plist,nrow = 5)
}

#'Plot Pairwise PCA plots.
#' @param object Seurat object
#' @param dims Compute N PC dims
#' @export
DimPCAPlot <- function(object,dims=1:10,feature){
  plist <- list()
  for (i in dims){
    plist[[as.character(i)]] <- FeatureScatter(object,paste0('PC_',i),feature) + geom_point(color='gray') +
      theme(legend.position="none",
            axis.text.x = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.x = element_blank(),
            axis.ticks.y = element_blank(),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            plot.title = element_text(size=8)
      ) +
      ggtitle(paste0('PC',i,' vs ',feature))
  }
  plot_grid(plotlist = plist,ncol = 5)
}


#'bi_getValues utility function.
#' @param object Seurat object
#' @param feature1 Feature 1
#' @param feature2 Feature 2
#' @param feature1.min Min UMI for a cell to be postive for feature. Set higher for more a restrictive selection
#' @param feature2.min Min UMI for a cell to be postive for feature. Set higher for more a restrictive selection
#' @param dims Dimensions to plot, must be a two-length numeric vector specifying x- and y-dimensions
#' @param reduction Which dimensionality reduction to use
#' @param pt.size Adjust point size
#' @param cols Vector of colors for Feature1+/Feature2+, feature1+,feature2+ and feature1-/feature2-
#' @param plotLineage
#' @param title Plot Title
#' @export

BiGenePlot <-
  function(object,
           feature1,
           feature2,
           feature1.min = 1,
           feature2.min = 1,
           dims = 1:2,
           reduction = "umap",
           pt.size = 0.1,
           cols = c("#E41A1C", "#377EB8", "#4DAF4A", 'grey75'),
           plotLineage = FALSE,
           title = NULL
  )
  {
    ###

    if(length(cols)!=4){
      stop("Please input a vector of 4 colors in the following oerder: Feature1+/Feature2+, feature1+,feature2+ and feature1-/feature2-")
    }



    feature1.name <- paste0(feature1, '+')
    feature2.name <- paste0(feature2, '+')
    feature.both.name <- paste0(feature1, '+/', feature2, '+')
    feature.none.name <- paste0(feature1, '-/', feature2, '-')

    dims <- paste0(Key(object = object[[reduction]]), dims)
    data <-
      FetchData(object = object,
                vars = c(dims, 'ident', feature1, feature2)) %>%
      mutate(value = case_when(
        (!!sym(feature1) >= feature1.min & !!sym(feature2) >= feature2.min) ~ feature.both.name,
        (!!sym(feature1) >= feature1.min & !!sym(feature2) < feature2.min) ~ feature1.name,
        (!!sym(feature1) < feature1.min & !!sym(feature2) >= feature2.min) ~ feature2.name,
        (!!sym(feature1) < feature1.min & !!sym(feature2) < feature2.min)  ~ feature.none.name
      )
      ) %>%
      mutate(value = factor(
        value,
        levels = c(
          feature.both.name,
          feature1.name ,
          feature2.name ,
          feature.none.name
        )
      )) %>%
      arrange(desc(value))

    p <- ggplot(data, aes_string(x = dims[1], y = dims[2])) +
      geom_point(aes(color = value)) +
      scale_color_manual(values = cols,drop = F) +
      theme_void() + coord_equal() +
      theme(
        plot.title = element_text(hjust = 0.5),
        legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
      )

    if(plotLineage){
      curved <-
        bind_rows(lapply(names(object@misc$sds$data@curves), function(x) {
          c <- slingCurves(object@misc$sds$data)[[x]]
          d <- as.data.frame(c$s[c$ord, seq_len(2)])
          d$curve <- x
          return(d)
        }))
      p <- p +geom_path(aes_string(dims[1], dims[2], linetype = "curve"), curved, size =1)
    }
p

  }

#'Function to plot multiple genes in a single Violin plot
#' @param object Seurat object
#' @param features genes to plot
#' @param group.by variable to group the cells by
#' @param cols Colors
#' @param orientation orientation to order the plots horizontally (single row) or vertically (single column)
#' @export
celltypeVlnPlot <- function(object, features,group.by='var_cluster',cols,orientation="vertical"){

  d <- FetchData(object,c(features,group.by)) %>% tidyr::gather(gene,signal,-`group.by`)
  d$gene <- factor(d$gene,levels=features)
  gp=eval(parse(text = paste0("d$",group.by,sep="")))
  d[,1] <- factor(gp,levels= sort(unique(gp)))
  if (orientation=="horizontal"){
    ggplot(d,aes_string(x=group.by,y="signal",color=group.by,fill=group.by)) +
      geom_violin() + facet_wrap(~gene,scales='free_x',nrow = 1) +
      theme_base() + coord_flip() +
      scale_fill_manual(values=cols) + scale_color_manual(values=cols) +
      scale_x_discrete(limits = levels(gp)) +
      theme(legend.position="none",
            axis.title.y=element_blank(),
            axis.title.x=element_blank(),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank(),
            strip.text.x = element_text(face = "bold.italic")
      )
  }else if(orientation =="vertical"){
    ggplot(d,aes_string(x=group.by,y="signal",color=group.by,fill=group.by)) +
      geom_violin() + facet_wrap(~gene,scales='free_x',ncol = 1) +
      theme_base() +
      scale_fill_manual(values=cols) + scale_color_manual(values=cols) +
      scale_x_discrete(limits = levels(gp)) +
      theme(legend.position="none",
            axis.title.y=element_blank(),
            axis.title.x=element_blank(),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank(),
            strip.text.x = element_text(face = "bold.italic")
      )
  }
}



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

#'Function to plot
#' @param object Seurat object
#' @param feature Feature to split on, ie sample
#' @param reduction Reduction used for the plot, umap default
#' @param ncol How many columns to use for the plot
#' @param cols Color palette
#' @return patchwork object
#' @import patchwork
#' @import purr
#' @import dplyr
#' @export




SplitMetaPlot <- function(object,reduction='umap',feature=NULL,ncol=4,col=NULL){
  if(is.null(col)){
    stop("Please provide color palette")
  }

  groups <- object@meta.data %>% pull(!!feature) %>% unique()
  if(length(groups) > 35) {
    stop("Too Many Groups")
  }
  plots <- groups %>% map(~DimPlot(object = object, reduction = "umap",label=F,pt.size = .1,cells.highlight = object@meta.data %>% filter(!!sym(feature)==!!.x) %>% rownames())  +
                            coord_equal() +
                            theme_void() +
                            ggtitle(.x) +
                            NoLegend()
  )

  if(!is.null(col)){
    plots <- map2(plots,col[1:length(groups)], ~.x + scale_color_manual(values=c('grey',.y)))
  }
  patchwork::wrap_plots(plots,ncol=ncol)
}


