# scExtras


Provides additional functions for Seurat v3  

* Pipeline tools
* Diffusionmap
* Slingshot
* Ligand-Receptor Analysis

## Requirements
```
devtools::install_github('chris-mcginnis-ucsf/DoubletFinder')
devtools::install_github("jokergoo/ComplexHeatmap")
devtools::install_github("Morriseylab/ligrec")

install.packages(c("Seurat","NMF","data.table","broom","quantreg","gam","parallelDist"))

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("kstreet13/slingshot")
BiocManager::install(c("destiny","scds","Biobase"))

```

## Install 
```
devtools::install_github('Morriseylab/scExtras')
```


## Run Pipeline 
Set all input parameters

```
outdir<-'Seurat' 
projectname<-'HumanLung' 
input10x <- c('STARSolo/HumanLungSolo.out/Gene/filtered/') 

org<-'mouse'
npcs<-40 
k=30

dir.create(outdir,recursive = T,showWarnings = F)
plotdir <- paste0(outdir,'/plots')
dir.create(plotdir,showWarnings = F)
qcdir <-paste0(outdir,'/qc')
dir.create(qcdir,showWarnings = F)
```
Read in 10x files, this can be H5 files or directories containg the mtx, barcode and gene files and create Seurat Object.

```
scrna = CreateSeuratObj(files=input10x,name=projectname)
```
Filter the data and run doublet detector.Doublet detection is performed using [scds](https://www.bioconductor.org/packages/release/bioc/html/scds.html) and [DoubletFinder](https://github.com/chris-mcginnis-ucsf/DoubletFinder)
```
scrna = RunQC(scrna,org=org,filter = T,LowerFeatureCutoff=200,UpperFeatureCutoff="MAD",UpperMitoCutoff=5,doubletdetection = T,dir=qcdir)
```
Normalize data and scale data, ccscale=T will regress out cell cycle effects. To run scTransform, use the [vignette](https://github.com/ChristophH/sctransform)  
```
scrna = processExper(scrna,ccscale=F,return_var_genes = F,org=org)
```
Perform PCA and save plots in QC dir. 
```
scrna= PCATools(scrna, npcs=npcs, jackstraw=T, plotdir = qcdir)
```
Perform Louvain Clustering and UMAP reduction 
```
npcs= 15 #set number of PC's to use based on elbow and jackstraw plots
scrna <- ClusterDR(scrna,dims=1:npcs,n.neighbors =20,min.dist=0.2, findallmarkers =F,res=0.4 )

```

## Diffusion Map
Seruat v3 removed the Diffusionmap dimension reduction routine. 
```
scrna <- RunDiffusion(scrna, dims=1:20)
```

## Trajectory Analysis using slingshot
 ```
 scrna <- runSlingshot(mes,reduction='umap',approx_points = 200,extend= "n",stretch=0)
 lineageDimPlot(scrna,reduction = "umap",group.by = "var_cluster",lineage = "all")
 
 

```

## Basic Ligand-Receptor Analysis 

```
scrna <- RunLigRec(scrna,org=org)
```

##Save seurat object as RDS file

```
saveRDS(scrna,paste0(outdir,'/',projectname,'.RDS'))
```


