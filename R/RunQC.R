#'Perform QC on single cell RNA-seq data with seurat, store all stats and return seurat object
#' @param dir Directory containing the cellranger files for input
#' @param org Organism
#' @param ccscale Specify whether to scale for cell cycle genes or not
#' @param filter Specify whether to filter data
#' @param LowerFeatureCutoff Lower value for Feature cutoff
#' @param UpperfeatureCutoff Upper bound, can be an integer or if set to MAD will be 3 median + 3*MAD
#' @param UpperMitoCutoff Upper Mito Pecent cuoff default is 0.05
#' @param sc.transform Define whether or not to apply scTransform function in seurat. T/F
#' @param doubletdetection Logical to specify if you want to run scds doublet detection
#' @return Seurat object
#' @import dplyr tidyr Seurat scds
#' @export
#' @examples
#' scrna = processExper(dir=prjdir,'test_prj',org='mouse',files=c("filtered_gene_bc_matrices/mm10"),ccscale=F,filter = T,LowerFeatureCutoff=200,UpperFeatureCutoff="MAD",UpperMitoCutoff=0.05)

RunQC <- function(dir,
                  name,
                  org='mouse',
                  files,
                  filter = T,
                  LowerFeatureCutoff=200,
                  UpperFeatureCutoff="MAD",
                  UpperMitoCutoff=5,
                  doubletdetection = F
){
  try(if(length(files)==0) stop("No files"))

  if(UpperFeatureCutoff!="MAD" & !is.integer(UpperFeatureCutoff)) {
    stop("Please use MAD and numeric cutoff for UpperFeatureCount")
  }


  if(length(files)==1){
    # Load the dataset
    if(dir.exists(files[1])){
      inputdata <- Read10X(data.dir =files[1])
    }else{
      inputdata <- Read10X_h5(filename =files[1])
    }
    colnames(inputdata) <- paste0(colnames(inputdata), '-',name)
    # Initialize the Seurat object with the raw (non-normalized data).
    object <- CreateSeuratObject(counts= inputdata, min.cells = 10, min.features = 200,project = name)
  }else{
    #Initialize the first object with the raw (non-normalized data) and add rest of the data
    if(dir.exists(files[1])){
      inputdata <- Read10X(data.dir =files[1])
    }else{
      inputdata <- Read10X_h5(filename =files[1])
    }
    colnames(inputdata) <- paste0(colnames(inputdata), '-',name, '-rep1')
    object <- CreateSeuratObject(counts= inputdata, min.cells = 10, min.features = 200, project = name)
    #cat('Rep1', length(object@cell.names), "\n")
    for(i in 2:length(files)){
      if(dir.exists(files[i])){
        tmp.data <- Read10X(data.dir =files[1])
      }else{
        tmp.data <- Read10X_h5(filename =files[i])
      }

      colnames(tmp.data) <- paste0(colnames(tmp.data), '-',name, '-rep',i)

      tmp.object <- CreateSeuratObject(counts= tmp.data, min.cells = 10, min.features = 200, project = name)
      #cat('Rep', i, ": ", length(tmp.object@cell.names), "\n", sep="")
      object <- merge(object, tmp.object, do.normalize = FALSE, min.cells = 0, min.features = 0)
    }

    # cat("merged: ", length(object@cell.names), "\n", sep="")
  }
  if(doubletdetection ==T){
    sce= as.SingleCellExperiment(object)
    sce = cxds(sce)
    sce = bcds(sce)
    sce = cxds_bcds_hybrid(sce)
    object=as.Seurat(sce)
  }
  if(org=='mouse'){
    object[["percent.mito"]] <- PercentageFeatureSet(object, pattern = "^mt-")
  }else{
    object[["percent.mito"]] <- PercentageFeatureSet(object, pattern = "^MT-")
  }

  write.csv(object@meta.data, file=paste(dir,"/percentmito.csv",sep=""))
  png(paste(dir,"/QC_Vlnplot.png",sep=""), width=10, height=6, units="in", res=300)
  print({VlnPlot(object = object, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"),
                 ncol = 3)})
  dev.off()

  object@misc[["filterstats"]] <- list()
  object@misc[["filterstats"]][['TotalCellsbeforefilteration']] <- dim(object)[2]

  if(filter){
    #Using a median + 3 MAD cutoff for high genes.
    if(UpperFeatureCutoff=="MAD"){
      UpperFeatureCutoff <- median(object$nFeature_RNA) + 3*mad(object$nFeature_RNA)
    }

    object@misc[["filterstats"]][['TotalSamples']] <- dim(object[[]][1])[1]

    cells.use <- colnames(object)[which(object[[]]['percent.mito'] < UpperMitoCutoff)]
    object@misc[["filterstats"]][['Mitofilter']] <- dim(object[[]][1])[1] - length(cells.use)
    object <- subset(object, cells = cells.use)

    cells.use <- colnames(object)[which(object[[]]['nFeature_RNA'] > LowerFeatureCutoff)]
    object@misc[["filterstats"]][['LowFeatureFilter']] <- dim(object[[]][1])[1] - length(cells.use)
    object <- subset(object, cells = cells.use)

    cells.use <- colnames(object)[which(object[[]]['nFeature_RNA'] < UpperFeatureCutoff)]
    object@misc[["filterstats"]][['HighFeatureFilter']] <- dim(object[[]][1])[1] - length(cells.use)
    object <- subset(object, cells = cells.use)

  }
  return(object)
}
