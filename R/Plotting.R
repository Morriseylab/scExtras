
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
#' @param features Compute N PC dims
#' @export

bigene_getValues <- function(object,features,limita,limitb){
  if(features !=2){
    stop("Please enter 2 features")
  }

  gene_values=FetchData(object=object,vars=features)
  colnames(gene_values) <- c('genea','geneb')
  # gene_values=as.data.frame(gene_values) %>%
  #   mutate(value =ifelse(genea>=limita[1] & geneb <limitb[1], gene_probes[1], ifelse(genea<limita[1] & geneb >=limitb[1],
  #                   gene_probes[2],ifelse(genea>=limita[1] & geneb >=limitb[1],"DoublePos","NULL"))))
  as.data.frame(gene_values) %>%
    mutate(value = ifelse(genea>=limita[1] & geneb>=limitb[1],
                          'both',
                          ifelse(genea>=limita[1] & geneb<limitb[1],
                                 gene_probes[1],
                                 ifelse(genea<=limita[1] & geneb>=limitb[1],
                                        gene_probes[2],
                                        'none')))
    ) #%>% select(value)
}


#'bi_getValues utility function.
#' @param object Seurat object
#' @param features 2 genes
#' @param x
#' @param y
#' @param limita
#' @param limitb
#' @param pt.size Adjust point size
#' @param title Plot Tttle
#' @param reduction Default is umap
#' @export

BiGenePlot <- function (object, features, x=1,y=2, limita=c(1,100), limitb=c(1,100), pt.size = 0.1,
                         title = NULL,reduction="umap")
{

  gene_values <- bigene_getValues(object,gene_probes,limita,limitb)
  projection=as.data.frame(eval(parse(text=paste("Embeddings(object, reduction =\"",type,"\")",sep=""))))
  colnames(projection) <- c("Component.1", "Component.2")
  proj_gene <- data.frame(cbind(projection, gene_values))
  #proj_gene$value = factor(proj_gene$value,levels=unique(proj_gene$value))
  proj_gene$value = factor(proj_gene$value,levels=c('both',gene_probes[1],gene_probes[2],'none'))
  proj_gene <- arrange(proj_gene, desc(value))

  dims <- paste0(Key(object = object[[reduction]]), dims)

  p=FetchData(object = object, vars = c(dims,groupby)) %>%
    ggplot(.,aes_string(x=dims[1],y=dims[2]))+geom_point(aes(color=!!sym(groupby))) +




  p <- ggplot(proj_gene, aes(Component.1, Component.2)) +
    geom_point(aes(colour = value), size = marker_size) +
    scale_color_manual(values=c("#E41A1C","#377EB8","#4DAF4A", 'grey75'),drop=F) +
    theme(legend.key.size = unit(10,"point")) + xlab(paste("Component", x)) +
    ylab(paste("Component", y))


  if (!is.null(title)) {
    p <- p + ggtitle(title)
  }
  p <- p + monocle_theme_opts() + theme(plot.title = element_text(hjust = 0.5),
                                        legend.position="bottom",
                                        legend.title=element_blank(),
                                        legend.text=element_text(size=14),
                                        panel.grid.major = element_blank(),
                                        panel.grid.minor = element_blank())

  return(p)
}














#' Visualize 'features' on a dimensional reduction plot
#'
#' Colors single cells on a dimensional reduction plot according to a 'feature'
#' (i.e. gene expression, PC scores, number of genes detected, etc.)
#'
#' @param order Boolean determining whether to plot cells in order of expression. Can be useful if
#' cells expressing given feature are getting buried.
#' @param features Vector of features to plot.
#' @param cols The two colors to form the gradient over. Provide as string vector with
#' the first color corresponding to low values, the second to high. Also accepts a Brewer
#' color scale or vector of colors. Note: this will bin the data into number of colors provided.
#' @param min.cutoff,max.cutoff Vector of minimum and maximum cutoff values for each feature,
#'  may specify quantile in the form of 'q##' where '##' is the quantile (eg, 'q1', 'q10')
#' @param split.by A factor in object metadata to split the feature plot by, pass 'ident'
#'  to split by cell identity'; similar to the old \code{FeatureHeatmap}
#' @param blend Scale and blend expression values to visualize coexpression of two features
#' @param blend.threshold The color cutoff from weak signal to strong signal; ranges from 0 to 1.
#' @param ncol Number of columns to combine multiple feature plots to, ignored if \code{split.by} is not \code{NULL}
#' @param combine Combine plots into a single gg object; note that if TRUE; themeing will not work when plotting multiple features
#' @param coord.fixed Plot cartesian coordinates with fixed aspect ratio
#' @param by.col If splitting by a factor, plot the splits per column with the features as rows.

#'
#' @return A ggplot object
#'
#' @importFrom grDevices rgb
#' @importFrom cowplot theme_cowplot
#' @importFrom RColorBrewer brewer.pal.info
#' @importFrom ggplot2 labs scale_x_continuous scale_y_continuous theme element_rect dup_axis guides
#' element_blank element_text margin scale_color_brewer scale_color_gradientn scale_color_manual coord_fixed
#' ggtitle
#' @examples
#' FeaturePlot(object = pbmc_small, features = 'PC_1')
#' @export
FeaturePlot2 <- function(
  object,
  features,
  dims = c(1, 2),
  cells = NULL,
  cols = c("lightgrey",  "blue"),
  pt.size = NULL,
  order = TRUE,
  min.cutoff = NA,
  max.cutoff = NA,
  reduction = NULL,
  split.by = NULL,
  shape.by = NULL,
  blend = FALSE,
  blend.threshold = 0.5,
  label = FALSE,
  label.size = 4,
  ncol = NULL,
  combine = TRUE,
  coord.fixed = FALSE,
  by.col = TRUE
) {
  no.right <- theme(
    axis.line.y.right = element_blank(),
    axis.ticks.y.right = element_blank(),
    axis.text.y.right = element_blank(),
    axis.title.y.right = element_text(
      face = "bold",
      size = 14,
      margin = margin(r = 7)
    )
  )
  if (is.null(reduction)) {
    default_order <- c('umap', 'tsne', 'pca')
    reducs <- which(default_order %in% names(object@reductions))
    reduction <- default_order[reducs[1]]
  }
  if (length(x = dims) != 2 || !is.numeric(x = dims)) {
    stop("'dims' must be a two-length integer vector")
  }
  if (blend && length(x = features) != 2) {
    stop("Blending feature plots only works with two features")
  }
  dims <- paste0(Key(object = object[[reduction]]), dims)
  cells <- cells %||% colnames(x = object)
  data <- FetchData(object = object, vars = c(dims, features), cells = cells)
  features <- colnames(x = data)[3:ncol(x = data)]
  min.cutoff <- mapply(
    FUN = function(cutoff, feature) {
      return(ifelse(
        test = is.na(x = cutoff),
        yes = min(data[, feature]),
        no = cutoff
      ))
    },
    cutoff = min.cutoff,
    feature = features
  )
  max.cutoff <- mapply(
    FUN = function(cutoff, feature) {
      return(ifelse(
        test = is.na(x = cutoff),
        yes = max(data[, feature]),
        no = cutoff
      ))
    },
    cutoff = max.cutoff,
    feature = features
  )
  check.lengths <- unique(x = vapply(
    X = list(features, min.cutoff, max.cutoff),
    FUN = length,
    FUN.VALUE = numeric(length = 1)
  ))
  if (length(x = check.lengths) != 1) {
    stop("There must be the same number of minimum and maximum cuttoffs as there are features")
  }
  brewer.gran <- ifelse(
    test = length(x = cols) == 1,
    yes = brewer.pal.info[cols, ]$maxcolors,
    no = length(x = cols)
  )
  data[, 3:ncol(x = data)] <- sapply(
    X = 3:ncol(x = data),
    FUN = function(index) {
      data.feature <- as.vector(x = data[, index])
      min.use <- SetQuantile(cutoff = min.cutoff[index - 2], data.feature)
      max.use <- SetQuantile(cutoff = max.cutoff[index - 2], data.feature)
      data.feature[data.feature < min.use] <- min.use
      data.feature[data.feature > max.use] <- max.use
      if (brewer.gran == 2) {
        return(data.feature)
      }
      data.cut <- if (all(data.feature == 0)) {
        0
      }
      else {
        as.numeric(x = as.factor(x = cut(
          x = as.numeric(x = data.feature),
          breaks = brewer.gran
        )))
      }
      return(data.cut)
    }
  )
  colnames(x = data)[3:ncol(x = data)] <- features
  rownames(x = data) <- cells
  data$split <- if (is.null(x = split.by)) {
    RandomName()
  } else {
    switch(
      EXPR = split.by,
      ident = Idents(object = object)[cells],
      object[[split.by, drop = TRUE]][cells]
    )
  }
  if (!is.factor(x = data$split)) {
    data$split <- factor(x = data$split)
  }
  if (!is.null(x = shape.by)) {
    data[, shape.by] <- object[[shape.by, drop = TRUE]]
  }
  plots <- vector(
    mode = "list",
    length = ifelse(
      test = blend,
      yes = 4,
      no = length(x = features) * length(x = levels(x = data$split))
    )
  )
  xlims <- c(floor(x = min(data[, dims[1]])), ceiling(x = max(data[, dims[1]])))
  ylims <- c(floor(min(data[, dims[2]])), ceiling(x = max(data[, dims[2]])))
  if (blend) {
    ncol <- 4
    color.matrix <- BlendMatrix(col.threshold = blend.threshold)
    colors <- list(
      color.matrix[, 1],
      color.matrix[1, ],
      as.vector(x = color.matrix)
    )
  }
  for (i in 1:length(x = levels(x = data$split))) {
    ident <- levels(x = data$split)[i]
    data.plot <- data[as.character(x = data$split) == ident, , drop = FALSE]
    if (blend) {
      data.plot <- cbind(data.plot[, dims], BlendExpression(data = data.plot[, features[1:2]]))
      features <- colnames(x = data.plot)[3:ncol(x = data.plot)]
    }
    for (j in 1:length(x = features)) {
      feature <- features[j]
      if (blend) {
        cols.use <- as.numeric(x = as.character(x = data.plot[, feature])) + 1
        cols.use <- colors[[j]][sort(x = unique(x = cols.use))]
      } else {
        cols.use <- NULL
      }
      plot <- SingleDimPlot(
        data = data.plot[, c(dims, feature, shape.by)],
        dims = dims,
        col.by = feature,
        order = order,
        pt.size = pt.size,
        cols = cols.use,
        shape.by = shape.by,
        label = label,
        label.size = label.size
      ) +
        #scale_x_continuous(limits = xlims) +
        #scale_y_continuous(limits = ylims) +
        theme_cowplot()
      if (length(x = levels(x = data$split)) > 1) {
        plot <- plot + theme(panel.border = element_rect(fill = NA, colour = 'black'))
        plot <- plot + if (i == 1) {
          labs(title = feature)
        } else {
          labs(title = NULL)
        }
        if (j == length(x = features) && !blend) {
          suppressMessages(
            expr = plot <- plot +
              scale_y_continuous(sec.axis = dup_axis(name = ident)) +
              no.right
          )
        }
        if (j != 1) {
          plot <- plot + theme(
            axis.line.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.text.y = element_blank(),
            axis.title.y.left = element_blank()
          )
        }
        if (i != length(x = levels(x = data$split))) {
          plot <- plot + theme(
            axis.line.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.text.x = element_blank(),
            axis.title.x = element_blank()
          )
        }
      } else {
        plot <- plot + labs(title = feature)
      }
      if (!blend) {
        plot <- plot + guides(color = NULL)
        if (length(x = cols) == 1) {
          plot <- plot + scale_color_brewer(palette = cols)
        } else if (length(x = cols) > 1) {
          plot <- suppressMessages(
            expr = plot + scale_color_gradientn(
              colors = cols,
              guide = "colorbar"
            )
          )
        }
      }
      if (coord.fixed) {
        plot <- plot + coord_fixed()
      }
      plot <- plot
      plots[[(length(x = features) * (i - 1)) + j]] <- plot
    }
  }
  if (blend) {
    blend.legend <- BlendMap(color.matrix = color.matrix)
    for (i in 1:length(x = levels(x = data$split))) {
      suppressMessages(expr = plots <- append(
        x = plots,
        values = list(
          blend.legend +
            scale_y_continuous(
              sec.axis = dup_axis(name = ifelse(
                test = length(x = levels(x = data$split)) > 1,
                yes = levels(x = data$split)[i],
                no = ''
              )),
              expand = c(0, 0)
            ) +
            labs(
              x = features[1],
              y = features[2],
              title = if (i == 1) {
                paste('Color threshold:', blend.threshold)
              } else {
                NULL
              }
            ) +
            no.right
        ),
        after = 4 * i - 1
      ))
    }
  }
  plots <- Filter(f = Negate(f = is.null), x = plots)
  if (combine) {
    if (is.null(x = ncol)) {
      ncol <- 2
      if (length(x = features) == 1) {
        ncol <- 1
      }
      if (length(x = features) > 6) {
        ncol <- 3
      }
      if (length(x = features) > 9) {
        ncol <- 4
      }
    }
    ncol <- ifelse(
      test = is.null(x = split.by) || blend,
      yes = ncol,
      no = length(x = features)
    )
    legend <- if (blend) {
      'none'
    } else {
      split.by %iff% 'none'
    }
    if (by.col & !is.null(x = split.by)) {
      plots <- lapply(X = plots, FUN = function(x) {
        suppressMessages(x +
                           theme_cowplot() + ggtitle("") +
                           scale_y_continuous(sec.axis = dup_axis(name = "")) + no.right
        )
      })
      nsplits <- length(x = levels(x = data$split))
      idx <- 1
      for (i in (length(x = features) * (nsplits - 1) + 1):(length(x = features) * nsplits)) {
        plots[[i]] <- suppressMessages(plots[[i]] + scale_y_continuous(sec.axis = dup_axis(name = features[[idx]])) + no.right)
        idx <- idx + 1
      }
      idx <- 1
      for (i in which(x = 1:length(x = plots) %% length(x = features) == 1)) {
        plots[[i]] <- plots[[i]] + ggtitle(levels(x = data$split)[[idx]])
        idx <- idx + 1
      }
      idx <- 1
      if (length(x = features) == 1) {
        for (i in 1:length(x = plots)) {
          plots[[i]] <- plots[[i]] + ggtitle(levels(x = data$split)[[idx]])
          idx <- idx + 1
        }
      }
      plots <- plots[c(do.call(
        what = rbind,
        args = split(x = 1:length(x = plots), f = ceiling(x = seq_along(along.with = 1:length(x = plots))/length(x = features)))
      ))]
      plots <- CombinePlots(
        plots = plots,
        ncol = nsplits,
        legend = legend
      )
    }
    else {
      plots <- CombinePlots(
        plots = plots,
        ncol = ncol,
        legend = legend,
        nrow = split.by %iff% length(x = levels(x = data$split))
      )
    }
  }
  return(plots)
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





