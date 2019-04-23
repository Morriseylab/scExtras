#'Run Diffusion map
#' @param object Seurat object
#' @param dims number of dimensions
#' @param reduction reductionm method. Defaults to pca
#' @param features vector of gene names
#' @param assay Assay name. Defaults to 'RNA'
#' @param max.dim maximum dimensions
#' @param q.use q. value cutoff
#' @param reduction.name Dimension Reduction method
#' @param reduction.key Dimension Reduction key
#' @import dplyr tidyr Seurat broom
#' @export
RunDiffusion <- function(
  object,
  dims = 1:5,
  reduction = 'pca',
  features = NULL,
  assay = 'RNA',
  max.dim = 3L,
  q.use = 0.01,
  reduction.name = "dm",
  reduction.key = "DM_",
  ...
) {
  if (!is.null(x = dims) || is.null(x = features)) {
    data.use <- Embeddings(object[[reduction]])[, dims]
    assay <- DefaultAssay(object = object[[reduction]])
  } else {
    data.use <- t(x = GetAssayData(object = object, slot = 'data', assay = assay)[features, ])
  }

  data.dist <- parallelDist::parDist(data.use)
  data.diffusion <- data.frame(destiny::DiffusionMap(data = as.matrix(data.dist)+1,n_eigs = max.dim)@eigenvectors)

  colnames(x = data.diffusion) <- paste0(reduction.key, 1:ncol(x = data.diffusion))
  rownames(x = data.diffusion) <-  rownames(data.use)
  # for (i in 1:max.dim) {
  #    x <- data.diffusion[, i]
  #   x <- MinMax(data = x, min = quantile(x = x, probs = q.use),
  #              quantile(x = x, probs = 1 - q.use))
  #  data.diffusion[, i] <- x
  #}

  assay <- DefaultAssay(object = object[[reduction]])

  dm.reduction <- CreateDimReducObject(
    embeddings = as.matrix(data.diffusion),
    key = reduction.key,
    assay = assay
  )


  object[[reduction.name]] <- dm.reduction



  #object <- LogSeuratCommand(object = object)
  return(object)
}
