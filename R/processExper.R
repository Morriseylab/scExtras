#'Preprocess single cell RNA-seq data with seurat and return seurat object
#' @param object Seurat object
#' @param ccscale Specify whether to scale for cell cycle genes or not
#' @param sc.transform Define whether or not to apply scTransform function in seurat. T/F
#' @return Seurat object
#' @import dplyr tidyr Seurat
#' @export
#' @examples
#' scrna = processExper(object=scrna,ccscale=F,sc.transform =F)

processExper <- function(object,ccscale=F,sc.transform=F){
  if(sc.transform==T){

    if(ccscale==T){
      object <- NormalizeData(object = object)

      #detection of variable genes
      #calculates the average expression and dispersion for each gene, places these genes into bins, and then calculates a z-score for dispersion within each bin
      object <-FindVariableFeatures(object = object,
                                    selection.method = "vst", nfeatures = 2000, verbose = FALSE)

      if(org=='human'){
        #Assign scores in the CellCycleScoring function.Stores S and G2/M scores in object@meta.data, along with the predicted classification of each cell in either G2M, S or G1 phase
        object <- CellCycleScoring(object = object, s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes)
      }else{
        m2h <- readr::read_csv(mouseorthologfile)
        cc.genes$s.genes <- m2h %>% filter(human_name %in% cc.genes$s.genes) %>% pull(mouse_name)
        cc.genes$g2m.genes <- m2h %>% filter(human_name %in% cc.genes$g2m.genes) %>% pull(mouse_name)
        object <- CellCycleScoring(object = object, s.features  = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes)
      }
      object <- SCTransform(object, verbose = FALSE, vars.to.regress=c("nCount_RNA","percent.mito","S.Score", "G2M.Score"))
    }else{
      object <- SCTransform(object, vars.to.regress = c("nCount_RNA","percent.mito"), verbose = FALSE)
    }
  }else{
    #normalize data
    object <- NormalizeData(object = object)

    #detection of variable genes
    #calculates the average expression and dispersion for each gene, places these genes into bins, and then calculates a z-score for dispersion within each bin
    object <-FindVariableFeatures(object = object,
                                  selection.method = "vst", nfeatures = 2000, verbose = FALSE)

    if(org=='human'){
      #Assign scores in the CellCycleScoring function.Stores S and G2/M scores in object@meta.data, along with the predicted classification of each cell in either G2M, S or G1 phase
      object <- CellCycleScoring(object = object, s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes)
    }else{
      m2h <- readr::read_csv(mouseorthologfile)
      cc.genes$s.genes <- m2h %>% filter(human_name %in% cc.genes$s.genes) %>% pull(mouse_name)
      cc.genes$g2m.genes <- m2h %>% filter(human_name %in% cc.genes$g2m.genes) %>% pull(mouse_name)
      object <- CellCycleScoring(object = object, s.features  = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes)
    }

    if(ccscale==T){
      #Scaling the data and removing unwanted sources of variation
      object <- ScaleData(object = object, vars.to.regress = c("nCount_RNA", "percent.mito","S.Score", "G2M.Score"))
    }else{
      object <- ScaleData(object = object, vars.to.regress = c("nCount_RNA", "percent.mito"))
    }
  }
  object <- LogSeuratCommand(object = object)
}
