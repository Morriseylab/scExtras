library(Seurat)
library(dplyr,quietly = T)
library(cowplot,quietly = T)
library(scExtras)
library(readxl)
library(ggplot2)

outdir<-'C2C_0D/Seurat' # All results will be saved here plots, RData file
projectname<-'C2C_0D' # specify project name,this will also be the Rdata file name
input10x <- c('STARSolo/C2C_0DSolo.out/Gene/filtered/') # dir(s) of the 10x output files, genes.tsv,barcodes.tsv

org<-'mouse'
npcs<-40 #How many inital PC dimensions to compute.
k=30 #This for nearest neighbors, 30 is default

dir.create(outdir,recursive = T,showWarnings = F)
plotdir <- paste0(outdir,'/plots')
dir.create(plotdir,showWarnings = F)
qcdir <-paste0(outdir,'/qc')
dir.create(qcdir,showWarnings = F)

#Create seurat object
scrna = CreateSeuratObj(files=input10x,name=projectname)

#Filter data and run double detector
scrna = RunQC(scrna,
              org=org,
              filter = T,
              LowerFeatureCutoff=200,
              UpperFeatureCutoff="MAD",
              UpperMitoCutoff=5,
              doubletdetection = T,
              dir=qcdir)

#normalize and scale data
 ##### Replace this section if running sctransform
scrna = processExper(scrna,ccscale=F,return_var_genes = F,org=org)
 ####

#Run PCA
scrna= PCATools(scrna, npcs=npcs, jackstraw=T, plotdir = qcdir)

#Run Clustering,Tsne and umap
npcs= 15 #set number of PC's to use based on elbow and jackstraw plots
scrna <- ClusterDR(scrna,dims=1:npcs,n.neighbors =20,min.dist=0.2, findallmarkers =F,res=0.4 )

#Find ligand receptors
scrna <- RunLigRec(scrna,org=org)

#Save seurat object as RDS file
saveRDS(scrna,paste0(outdir,'/',projectname,'.RDS'))

