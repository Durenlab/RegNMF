read_ATAC_GEX.default<-function(foldername){
  if(typeof(foldername)!="character"){
    stop("foldername should be character\n")
  }

  else{

    filename=paste0(foldername,"matrix.mtx")
    a=read.table(filename,sep = "",header = T,comment.char = '%')
    filename=paste0(foldername,"features.tsv")
    C=read.table(filename,sep='\t')
    chr=unique(C[,4])
    chr=chr[grep("^chr",chr)]
    isATAC=match(C[,3],"Peaks")
    isATAC[is.na(isATAC)]=0
    isATAC=as.logical(isATAC)
    PeakName=C[isATAC,2]
    Symbol=C[!isATAC,2]
    features=C[,2]

    chr_idx=match(C[,4],chr)
    chr_idx[is.na(chr_idx)]=0


    Symbol_location=matrix(chr_idx[!isATAC],nrow = length(chr_idx[!isATAC]))
    Symbol_location=cbind(Symbol_location,C[!isATAC,5])
    PeakName_location=matrix(chr_idx[isATAC],nrow = length(chr_idx[isATAC]))
    PeakName_location=cbind(PeakName_location,C[isATAC,5])
    filename=paste0(foldername,"barcodes.tsv")
    barcode=read.table(filename,sep="\t")

    ##rna
    f=match(features,Symbol)
    d=f;
    d[is.na(d)]=0
    d=as.logical(d)

    tmp=match(a[,1],which(d==TRUE))
    tmp[is.na(tmp)]=0
    tmp=as.logical(tmp)
    a_rna=a[tmp,]
    a_rna[,1]=f[a_rna[,1]]
    E=sparseMatrix(a_rna[,1],a_rna[,2],x=log2(1+a_rna[,3]),dims = c(length(Symbol),length(barcode[,1])))

    ##atac
    f=match(features,PeakName)
    d=f;
    d[is.na(d)]=0
    d=as.logical(d)

    tmp=match(a[,1],which(d==TRUE))
    tmp[is.na(tmp)]=0
    tmp=as.logical(tmp)
    a_atac=a[tmp,]
    a_atac[,1]=f[a_atac[,1]]
    O=sparseMatrix(a_atac[,1],a_atac[,2],x=log10(1+a_atac[,3]),dims = c(length(PeakName),length(barcode[,1])))

    ##
    idx=Symbol_location[,1]>0
    Symbol=Symbol[idx];
    Symbol_location=Symbol_location[idx,]
    E=E[idx,]
    idx=PeakName_location[,1]>0
    PeakName=PeakName[idx]
    PeakName_location=PeakName_location[idx,]
    O=O[idx,]
    ##non-zero 10 cell
    idx=rowSums(as.matrix(E)>0 )>0;
    Symbol=Symbol[idx];
    Symbol_location=Symbol_location[idx,]
    E=E[idx,]
    idx=rowSums(as.matrix(O)>0)>10;
    PeakName=PeakName[idx]
    PeakName_location=PeakName_location[idx,]
    O=O[idx,]

    return(list(E=E,
                O=O,
                PeakName=PeakName,
                barcode=barcode,
                Symbol_location=Symbol_location,
                Symbol=Symbol,
                Peak_location=PeakName_location,
                chr=chr))
  }
}
