#'Create Seurat object from 10x output
#' @param files input files
#' @name projectname Project name
#' @return Seurat object
#' @import dplyr tidyr Seurat scds
#' @export
#' @examples
#' scrna = CreateSeuratObj(dir=prjdir,name='test_prj',files=c("10x/Gene/filtered"))

CreateSeuratObj <- function(files,name
){
  try(if(length(files)==0) stop("No files"))

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
    for(i in 2:length(files)){
      if(dir.exists(files[i])){
        tmp.data <- Read10X(data.dir =files[1])
      }else{
        tmp.data <- Read10X_h5(filename =files[i])
      }

      colnames(tmp.data) <- paste0(colnames(tmp.data), '-',name, '-rep',i)

      tmp.object <- CreateSeuratObject(counts= tmp.data, min.cells = 10, min.features = 200, project = name)
      object <- merge(object, tmp.object, do.normalize = FALSE, min.cells = 0, min.features = 0)
    }
  }
  return(object)
}
