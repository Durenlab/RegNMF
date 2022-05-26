# Integrate with Seurat

Here we demonstrate how to use reduced Dimension (H) from RegNMF in Seurat, using dataset of 10K human 
peripheral blood mononuclear cells (PBMCs) freely available from 10X Genomics. 

In this tutorial we will do following steps:
* create a seuart Object contain RNA-seq and ATAC-seq data
* add reduced dimension matrix H into this seurat object
* run umap on the reduced dimenson 
* clustering on the reduced dimension
* plot a umap visualization
* find the cluster specific differential expressing genes
### Load in packages
```R
library(Seurat)
library(ggplot2)
```
### Input the multiome data

Change the ~ to the directory of your data stored

Here we used a subset of the pbmc10k cells, *this is not required. 

* the barcode of cells used in this example is provided in barcode_use.txt
```R
#Read in h5 file
inputdata.10x <- Read10X_h5("~/pbmc_10k/pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.h5")
#split the data
rna_counts <- inputdata.10x$`Gene Expression`
atac_counts <- inputdata.10x$Peaks
#flifer by the barcode
barcode_use<-read.table("~/pbmc_10k/barcode_use.txt")
#filter rna
barcode<-rna_counts@Dimnames[[2]]
idx<-match(barcode_use$V1,barcode)
rna_counts<-rna_counts[,idx]
#filter atac
atac_counts <- inputdata.10x$Peaks
barcode<-atac_counts@Dimnames[[2]]
idx<-match(barcode_use$V1,barcode)
atac_counts<-atac_counts[,idx]
```

### Create a Seurat object 
```R
#create seurat object separately
pbmc <- CreateSeuratObject(counts = rna_counts,assay = "RNA")
pbmc.atac <- CreateSeuratObject(counts = atac_counts,assay = "Peaks")
# combine into 1 seurat object with 2 assays
pbmc@assays$Peaks<-pbmc.atac$Peaks
pbmc@meta.data$nCount_Peaks<-pbmc.atac@meta.data$nCount_Peaks
pbmc@meta.data$nFeature_Peaks<-pbmc.atac@meta.data$nFeature_Peaks
```
### Preparing scREG reduced dimension --- H

Here I used pre-saved result, you can just use the H matrix form RegNMF results
```R
H<-read.table("~/pbmc_10k/H.txt",header = FALSE, sep = "\t")
H=as.matrix(t(H))
rownames(H)<-barcode_use$V1
Dimnames2<-matrix(0,dim(H)[2],1)
for (i in (1:dim(H)[2])){
Dimnames2[i]<-paste0('Dim',i)
}
colnames<-as.character(Dimnames2)
colnames(H)<-colnames
```
### Put H into the seurat object
```R
cell.embeddings=H
reduction.data<- CreateDimReducObject(embeddings = cell.embeddings,key = "Dim",assay = "RNA")
pbmc@reductions$RegNMF<-reduction.data
```
### Now you should be able to run Umap and clustering & the downstream analysis
```R
pbmc <- RunUMAP(pbmc, reduction = "RegNMF", dims = 1:100, reduction.name = "umap.RegNMF")
pbmc <- FindNeighbors(pbmc, reduction = "RegNMF", dims = 1:100)###
pbmc <- FindClusters(pbmc,graph.name = "RNA_snn", resolution = 0.5)
```
### Plot Umap
```R
DimPlot(pbmc, reduction = "umap.RegNMF", group.by = "seurat_clusters", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("scREG")
```

You can of course run PCA and UMAP on asssay RNA or assay Peaks, using seurat,and plot the Umap for each.

### Find all markers of cluster 2
```R
cluster2.markers <- FindMarkers(pbmc, ident.1 = 2, min.pct = 0.25)
head(cluster2.markers, n = 5)
```
