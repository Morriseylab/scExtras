#'Preprocess single cell RNA-seq data with seurat and return seurat object
#' @param dir Directory containing the cellranger files for input
#' @param name Project name
#' @param org Organism
#' @param files path to the dir having the barcodes.tsv, genes.tsv and matrix.mtx files
#' @param ccscale Specify whether to scale for cell cycle genes or not
#' @param filter Specify whether to filter data
#' @param LowerFeatureCutoff Lower value for Feature cutoff [200]
#' @param UpperfeatureCutoff Upper bound, can be an integer or if set to MAD will be 3 median + 3*MAD
#' @return Seurat object
#' @import dplyr tidyr Seurat
#' @export
#' @examples
#' scrna = processExper(dir=prjdir,'test_prj',org='mouse',files=c("filtered_gene_bc_matrices/mm10"),ccscale=F,filter = T)

processExper <- function(dir,name,
                         org='mouse',
                         files,
                         ccscale=F,
                         filter = T,
                         LowerFeatureCutoff=200,
                         UpperFeatureCutoff="MAD"

                         ){
  try(if(length(files)==0) stop("No files"))

  try(if(UpperFeatureCutoff=="MAD" || !(is.numeric(UpperFeatureCutoff))) stop("Please use MAD and numeric cutoff for UpperFeatureCount"))



  if(length(files)==1){
    # Load the dataset
    inputdata <- Read10X(data.dir =files[1])
    colnames(inputdata) <- paste0(colnames(inputdata), '-',name)
    # Initialize the Seurat object with the raw (non-normalized data).
    object <- CreateSeuratObject(counts= inputdata, min.cells = 10, min.features = 200,project = name)
  }else{
    #Initialize the first object with the raw (non-normalized data) and add rest of the data
    inputdata <- Read10X(data.dir =files[1])
    colnames(inputdata) <- paste0(colnames(inputdata), '-',name, '-rep1')
    object <- CreateSeuratObject(counts= inputdata, min.cells = 10, min.features = 200, project = name)
    #cat('Rep1', length(object@cell.names), "\n")
    for(i in 2:length(files)){
      tmp.data <- Read10X(data.dir =files[i])
      colnames(tmp.data) <- paste0(colnames(tmp.data), '-',name, '-rep',i)

      tmp.object <- CreateSeuratObject(counts= tmp.data, min.cells = 10, min.features = 200, project = name)
      # cat('Rep', i, ": ", length(tmp.object@cell.names), "\n", sep="")
      object <- merge(object, tmp.object, do.normalize = FALSE, min.cells = 0, min.features = 0)
    }
    #cat("merged: ", length(object@cell.names), "\n", sep="")
  }



  if(org=='mouse'){
    mito.features <- grep(pattern = "^mt-", x = rownames(x = object), value = TRUE)
  }else{
    mito.features <- grep(pattern = "^MT-", x = rownames(x = object), value = TRUE)
  }


  percent.mito <- Matrix::colSums(x = GetAssayData(object = object, slot = "counts")[mito.features,
                                                                                     ])/Matrix::colSums(x = GetAssayData(object = object, slot = "counts"))

  # The [[ operator can add columns to object metadata, and is a great place
  # to stash QC stats
  object[["percent.mito"]] <- percent.mito
  VlnPlot(object = object, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"),
          ncol = 3)

  if(filter){
    #Using a median + 3 MAD cutoff for high genes.
    if(UpperFeatureCutoff=="MAD"){
      UpperFeatureCutoff <- median(object$nFeature_RNA) + 3*mad(object$nFeature_RNA)
    }

    object <- subset(object, subset = LowerFeatureCutoff > 200 & percent.mito < 0.05 & nFeature_RNA < UpperFeatureCutoff)
    object@misc[[stats]] <- list(featurecutofflow = 200,
                                   featurecutoffhigh = featureCutoff )

  }

  #normalize data
  object <- NormalizeData(object = object)

  #detection of variable genes
  #calculates the average expression and dispersion for each gene, places these genes into bins, and then calculates a z-score for dispersion within each bin
  object <-FindVariableFeatures(object = object,
                                selection.method = "vst", nfeatures = 2000, verbose = FALSE)


  if(ccscale==T){
    if(org=='human'){
      #Assign scores in the CellCycleScoring function.Stores S and G2/M scores in object@meta.data, along with the predicted classification of each cell in either G2M, S or G1 phase
      object <- CellCycleScoring(object = object, s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes)
    }else{
      m2h <- readr::read_csv(mouseorthologfile)
      cc.genes$s.genes <- m2h %>% filter(human_name %in% cc.genes$s.genes) %>% pull(mouse_name)
      cc.genes$g2m.genes <- m2h %>% filter(human_name %in% cc.genes$g2m.genes) %>% pull(mouse_name)
      object <- CellCycleScoring(object = object, s.features  = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes)
    }
    #Scaling the data and removing unwanted sources of variation
    object <- ScaleData(object = object, vars.to.regress = c("nUMI", "percent.mito","S.Score", "G2M.Score"))
  }else{
    object <- ScaleData(object = object, vars.to.regress = c("nUMI", "percent.mito"))
  }
  return(object)
}
