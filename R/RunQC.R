#'Perform QC on single cell RNA-seq data with seurat, store all stats and return seurat object
#' @param object Seurat object
#' @param org Organism (mouse/human)
#' @param filter Specify whether to filter data
#' @param LowerFeatureCutoff Lower value for Feature cutoff
#' @param UpperfeatureCutoff Upper bound, can be an integer or if set to MAD will be 3 median + 3*MAD
#' @param UpperMitoCutoff Upper Mito Pecent cuoff default is 0.05
#' @param doubletdetection Logical to specify if you want to run scds doublet detection
#' @param dir Output directory for QC files
#' @return Seurat object
#' @import dplyr tidyr Seurat scds
#' @export
#' @examples
#' scrna = RunQC(object=scrna,org="mouse",filter = T,LowerFeatureCutoff=200,UpperFeatureCutoff="MAD",UpperMitoCutoff=0.05,doubletdetection = T)

RunQC <- function(object,
                  org='mouse',
                  filter = T,
                  LowerFeatureCutoff=200,
                  UpperFeatureCutoff="MAD",
                  UpperMitoCutoff=5,
                  doubletdetection = T,
                  dir
){

  if(UpperFeatureCutoff!="MAD" & !is.integer(UpperFeatureCutoff)) {
    stop("Please use MAD and numeric cutoff for UpperFeatureCount")
  }

  if(doubletdetection ==T){
    sce= as.SingleCellExperiment(object)
    sce = cxds(sce,estNdbl=T)
    sce = bcds(sce,estNdbl=T)
    sce = cxds_bcds_hybrid(sce,estNdbl=T)
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
