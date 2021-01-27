# scExtras


Provides additional functions for Seurat v3  

* Pipeline tools
* Diffusionmap
* Slingshot
* Ligand-Receptor Analysis

## Requirements
```
devtools::install_github('chris-mcginnis-ucsf/DoubletFinder')

BiocManager::install("kstreet13/slingshot")
```

## Install 
```
devtools::install_github('Morriseylab/scExtras')
```


## Run Pipeline 

Read in 10x files, this can be H5 files or directories containg the mtx, barcode and gene files. 

```
input10x <- c('Rep1.H5','Rep2.H5')
outdir<-'Seurat'
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
scrna= PCATools(scrna, npcs=30, jackstraw=F, plotdir = 'QC')
```
Perform Louvain Clustering and UMAP redcution 
```
scrna <- ClusterDR(scrna,dims=1:npcs,n.neighbors =k)
```


## Diffusionmap
Seruat v3 removed the Diffusionmap dimension reduction routine. 

## Trajectory Analysis using slingshot


## Basic Ligand-Receptor Analysis 





