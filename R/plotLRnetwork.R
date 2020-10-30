#'Create ligand receptor network plot
#' @param results dataframe containing the fields- interacting_pair,Receptor,Ligand,Receptor_cluster,Ligand_cluster
#' @param filtermin minimum number of interactions between the nodes
#' @param filtermax maximum number of interactions between the nodes
#' @import dplyr tidyr igraph tibble
#' @export
#' @return Returns network plot with ligand-receptor pair data

plotLRnetwork <- function(results, filtermin=0,filtermax=0){
  edges=results %>% dplyr::select(Receptor_cluster,Ligand_cluster)
  e2=as.data.frame(table(edges[,1:2]))
  filtermin=ifelse(filtermin==0,min(e2$Freq),filtermin)
  filtermax=ifelse(filtermax==0,max(e2$Freq),filtermax)
  e2=e2[e2$Freq>= filtermin & e2$Freq<= filtermax,]
  e2$pair=paste(e2$Receptor_cluster,"_",e2$Ligand_cluster,sep="")
  results$pair=paste(results$Receptor_cluster,"_",results$Ligand_cluster,sep="")
  results=results[results$pair %in% e2$pair,]

  rec <- results %>% distinct(Receptor_cluster) %>% rename(label = Receptor_cluster)
  lig <- results %>% distinct(Ligand_cluster) %>% rename(label = Ligand_cluster)
  nodes <- full_join(rec,lig, by = "label")
  nodes <- nodes %>% rowid_to_column("id")
  col=cpallette[1:nrow(nodes)]
  nodes$color=col
  perpair <- results %>% group_by(Receptor_cluster, Ligand_cluster) %>% summarise(freq = n()) %>% ungroup()
  edges <- perpair %>%  left_join(nodes, by = c("Receptor_cluster" = "label")) %>% rename(to = id)
  edges <- edges %>% left_join(nodes, by = c("Ligand_cluster" = "label")) %>% rename(from = id)
  edges <- dplyr::select(edges, from, to, freq)
  edges=left_join(edges,nodes,by=c("from"="id")) %>% dplyr::select(-label)
  edge.col=edges$color
  edge.lab=as.character(edges$freq)
  OldRange = (max(edges$freq) - min(edges$freq))
  NewRange = 8-2
  edges$width = (((edges$freq - min(edges$freq)) * NewRange) / OldRange) + 1.5
  width=edges$width

  routes_igraph <- graph_from_data_frame(d = edges, vertices = nodes, directed = TRUE)
  curves <-autocurve.edges2(routes_igraph)
  plot(routes_igraph, edge.arrow.size = 0.2,vertex.label.color="black",edge.label.color="black",vertex.color=col,edge.color=edge.col,edge.width=width,edge.arrow.width=9.5,edge.arrow.size=9.5,edge.curved=curves)
}
