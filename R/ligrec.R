#'Get ligand receptor pairs between all clusters
#' @param object Seurat object
#' @param group.by grouping variable
#' @param org organism
#' @param perc minimum percentage of cells that should express the ligand or receptor
#' @import dplyr tidyr Seurat
#' @export
#' @return Returns the object with object@misc[[ligrecres]] containing the matrix of lignad-receptor pair data
RunLigRec <- function(object,group.by='ident',org,perc=30){
  #get grouping variable
  var=as.character(group.by)
  #genes=fread("data/ligrecgenes.txt",header = TRUE)
  if(org=="mouse"){
    data('Mm_PairsLigRec',package="ligrec")
    rl = mm
  }else if(org=="human"){
    data('Hs_PairsLigRec',package="ligrec")
    rl=hs
  }

  genes <- intersect(rownames(GetAssayData(object = object, slot = "counts",assay='RNA')), unique(c(as.character(rl$ligand),as.character(rl$receptor))))

  #For all unique genes in the ligrec list, get their expression value for all cells and the groups the cells belong to
  da=DefaultAssay(object)
  if(da == "SCT"){
    DefaultAssay(object) <- "RNA"
    my.data <- cbind(FetchData(object,c(var)), FetchData(object = object,genes,slot="counts"))
    DefaultAssay(object) <- "SCT"
  }else{
    my.data <- cbind(FetchData(object,c(var)), FetchData(object = object,genes,slot="counts"))
  }
  colnames(my.data)[1]= "clust"
  perc=perc/100
  result=data.frame()
  res=data.frame()
  my.data$clust= factor(my.data$clust, levels= unique(my.data$clust))
  #loop over each cluster to find pairs
  for(i in 1:(length(levels(my.data$clust)))){
    for(j in 1:(length(levels(my.data$clust)))){
      #from the large martix, subselect receptor and lig subgoups (if i=1 and j=2, keep cells in grps 1 and 2)
      test=my.data[my.data$clust==levels(my.data$clust)[i] | my.data$clust==levels(my.data$clust)[j],]
      #Subselect genes in receptor list in cells in rec subgroup (say 1)
      R_c1=test[test$clust==levels(my.data$clust)[i] ,(colnames(test) %in% rl$receptor)]
      #Subselect genes in ligand list in cells in lig subgroup (say 2)
      L_c2=test[test$clust==levels(my.data$clust)[j] , (colnames(test) %in% rl$ligand)]
      if(nrow(R_c1)!=0 &nrow(L_c2)!=0){
        #keep genes that are expressed in more than user-input percent of the cells
        keep1 = colSums(R_c1>0)>=perc*dim(R_c1)[1]
        keep2 = colSums(L_c2>0)>=perc*dim(L_c2)[1]
        R_c1=R_c1[,keep1]
        L_c2=L_c2[,keep2]
        #get list of lig-rec pairs
        res=rl[(rl$ligand %in% colnames(L_c2)) & (rl$receptor %in% colnames(R_c1)),]
      }else{}
      if(nrow(res)!=0){
        res$Receptor_cluster=levels(my.data$clust)[i]
        res$Lig_cluster=levels(my.data$clust)[j]
        result=rbind(result,res)
      }else{result=result}
    }
  }
  # get final list of all lig-rec pairs
  object@misc[['ligrecres']] <- result
  object
}
