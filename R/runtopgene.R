#'Run topGene API and convert it to enrichResult
#' @param genelist vector of genelist
#' @import httr jsonlite dplyr tidyr multienrichjam
#' @export

#g=c("FLDB","APOE","ENSG00000113196","ENSMUSG00000020287")

runtopgene <- function(genelist){
  pc_json <- list(symbols=genelist)
  res <- POST("https://toppgene.cchmc.org/API/lookup"
              , body = pc_json
              , encode = "json")
  appData <- content(res)
  genes=appData$Genes
  genes=as.data.frame(do.call(rbind, genes))
  entrez=genes$Entrez
  pc_json2 <- list(Genes = entrez)
  res2 <- POST("https://toppgene.cchmc.org/API/enrich"
               , body = pc_json2
               , encode = "json")
  appData2 <- content(res2)
  df=appData2$Annotations
  df2=as.data.frame(do.call(rbind, df))
  df3=df2 %>% tidyr::unnest(Genes)
  df4=df3 %>% tidyr::unnest_wider(Genes)
  df5=df4 %>% group_by(ID) %>% mutate(geneID = paste(Entrez, collapse="/"))
  df5=df5 %>% select(-Entrez) %>% group_by(ID) %>% mutate(Symbol = paste(Symbol, collapse=","))

  df5=df5[!duplicated(df5),]
  df5$GeneRatio=paste0(df5$GenesInTermInQuery,"/",df5$GenesInTerm,sep="")
  df5$BgRatio=paste0(df5$GenesInQuery,"/",df5$TotalGenes,sep="")
  df5 = df5 %>% rename('Description'='Name','pvalue'='PValue','qvalue'='QValueBonferroni','Count'='GenesInTermInQuery')
  df6=data.frame(lapply(df5, function(x) unlist(x)))
  rownames(df6)=df6$ID
  edo=enrichDF2enrichResult(df6,pAdjustMethod = "BH",keyColname = "ID", geneColname = "geneID", pvalueColname = "QValueFDRBH", descriptionColname = "Description", pvalueCutoff = 0.05)
  return(edo)
}
