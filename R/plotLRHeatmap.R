#'Create ligand receptor heatmap plot
#' @param result dataframe containing the fields- interacting_pair,Receptor,Ligand,Receptor_cluster,Ligand_cluster
#' @param clusterby clusterby option for heatmap. Choose from "both", "row","column" and "none"
#' @import dplyr tidyr RColorBrewer NMF tibble
#' @export
#' @return Returns heatmap with ligand-receptor pair data

plotLRHeatmap <- function(result,clusterby){
  dist2 <- function(x, ...) {as.dist(1-cor(t(x), method="pearson"))}
  result=result %>% dplyr::select(interacting_pair,Receptor,Ligand,Receptor_cluster:Ligand_cluster)
  tab=table(result[,4:5])
  tab=as.data.frame(tab)
  tab = tab %>% spread(Receptor_cluster,Freq)
  rownames(tab)=tab$Ligand_cluster
  tab= tab %>% dplyr::select(-Ligand_cluster)
  tab=Filter(var, tab)
  if(clusterby=="both"){
    row=TRUE
    column=TRUE
  }else if(clusterby=="row"){
    row=TRUE
    column=FALSE
  }else if(clusterby=="column"){
    row=FALSE
    column=TRUE
  }else if(clusterby=="none"){
    row=NA
    column=NA
  }
  aheatmap(as.matrix(tab),distfun=dist2,Rowv=row,Colv=column,col = colorRampPalette(brewer.pal(n = 9,"YlGnBu"))(30),main= "Receptor genes (x) vs Ligand genes (y)")
}
