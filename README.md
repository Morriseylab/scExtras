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

Read in 10x files, this can be H5 files or directories containg the mtx, barcode and gene files. Doublet detection is performed using [scds](https://www.bioconductor.org/packages/release/bioc/html/scds.html) and [DoubletFinder](https://github.com/chris-mcginnis-ucsf/DoubletFinder)

```
input10x <- c('Rep1.H5','Rep2.H5')
outdir<-'Seurat'
qcdir <- paste0(outdir,'/QC')
org<-'human'
projectname<- 'HumanLung1'

scrna= RunQC(dir=outdir,org=org,name=projectname,files=input10x ,filter=T, doubletdetection = T,UpperMitoCutoff=10)
```
Normalize data and scale data, cscale=T will regress out cell cycle effects and sc.transform will use [SCTransform](https://github.com/ChristophH/sctransform)  
```
scrna = processExper(scrna ,ccscale = T, sc.transform = T)
```
Perform PCA and save plots in QC dir. 
```
scrna= PCATools(scrna, npcs=30, jackstraw=F, plotdir = qcdir)
```
Perform Louvain Clustering and UMAP redcution 
```
scrna <- ClusterDR(scrna,dims=1:npcs,n.neighbors =k)
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



