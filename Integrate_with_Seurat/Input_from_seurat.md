
# Make input for RegNMF from seurat object
## Libraries
```r
library(Seurat)
library(GeneBreak)
library(Matrix, warn.conflicts = FALSE, quietly = TRUE)
library(dplyr)
library(tidyr)
library(stringr)
library(gtools)
library(RegNMF)
```
- I assume you are using paired multiome data and you have seurat object ready.
- You can conduct all data preprocessing steps that you need using seurat, then come to this part
## load seurat object
```r
seuratObj<-readRDS("~/seuratObj.rds") # replace by the path to your seurat object
#first extract the RNA and ATAC respectively
rna=seuratObj@assays$RNA@data # change data to count if you wish to use raw count
atac=seuratObj@assays$ATAC@data
# filter the cells based on the barcode list if needed
barcode_use<-read.table("~/barcode_use.txt")
#filter rna
barcode_rna<-rna@Dimnames[[2]]
idx<-match(barcode_use$V1,barcode_rna)
rna<-rna[,idx]
#rna<-rna[,!idx] #if barcode given are cells to trim out
#filter atac
barcode_atac<-atac@Dimnames[[2]]
idx<-match(barcode_use$V1,barcode_atac)
atac<-atac[,idx]
#atac<-atac[,!idx] #if barcode given are cells to trim out
```
## prepare the inputs from seurat object
```r
#E,O
E=as.matrix(rna)
O=as.matrix(atac)
#symbol,peak
symbol=data.frame(symbol=rna@Dimnames[[1]])
peakname=data.frame(symbol=atac@Dimnames[[1]])
#symbol location
data(ens.gene.ann.hg38)
gene_location=data.frame(cbind(ens.gene.ann.hg38$Gene,ens.gene.ann.hg38$Chromosome,ens.gene.ann.hg38$Start,ens.gene.ann.hg38$End))
gene_location$X2<- sub("^", "chr",gene_location$X2)
A=match(symbol,gene_location$X1)
keepid=!is.na(A)
symbol_location=gene_location[A[keepid],c(2,3)]
chr_idx=match(symbol_location[,1],chr)
chr_idx[is.na(chr_idx)]=0
symbol_location=cbind(as.numeric(chr_idx),as.numeric(symbol_location$X3))
#peak location
peak_loc=data.frame(peak_location=atac@Dimnames[[1]])
peak_location=peak_loc %>% separate(peak_location, c("chr","start","end"), sep ="-")
peak_location=peak_location[,c(1,2)]
chr=unique(peak_location[,1])
chr=chr[grep("^chr",chr)]
chr=chr[mixedorder(chr)]
chr_idx=match(peak_location[,1],chr)
chr_idx[is.na(chr_idx)]=0
peak_location=cbind(as.numeric(chr_idx),as.numeric(peak_location$start))
```
Run the RegNMF
```r
W123H=RegNMF(E=E, 
             O=O, 
             Symbol=symbol, 
             PeakName=peakname, 
             Symbol_location=symbol_location, 
             Peak_location=peak_location)
```
