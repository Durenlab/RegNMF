SplitGroup.default<-function(foldername,barcord,W3,H,Reg_symbol_name,Reg_peak_name,cluster){
  clustern=length(unique(cluster))
  barcord_cluster=data.frame(barcord=barcord,cluster=cluster)
  bfilename=paste0(foldername,"barcord_cluster.bed")
  write.table(barcord_cluster,bfilename,sep="\t",col.names = F,row.names = F,quote = FALSE)
  chr=c()
  peaks=c()
  peake=c()
  for (i in 1:length(Reg_peak_name)) {
    a=strsplit(Reg_peak_name[i],':')
    b=strsplit(a[[1]][2],'-')
    chr[i]=a[[1]][1]
    peaks[i]=b[[1]][1]
    peake[i]=b[[1]][2]
  }
  df=data.frame(chr=chr,peaks=peaks,peake=peake,symbolName=Reg_symbol_name)

  H_norm=H/sqrt(rowSums(H*H))
  W3_norm=t(t(W3)*sqrt(rowSums(H*H)))
  H_w=matrix(nrow = nrow(H_norm),ncol = clustern)

  for (i in 1:clustern) {
    H_w[,i]=rowMeans(H_norm[,cluster==i])

  }

  W3_cluster=W3_norm%*%H_w
  RegFolderName=paste0(foldername,"old_Reg_cluster/")
  if(!dir.exists(RegFolderName)){
    dir.create(RegFolderName,recursive = TRUE )
  }
  for (i in 1:clustern) {
    topk=order(W3_cluster[,i])[1:10000]
    outdf=df[topk,]
    outdf$Reg=W3_cluster[topk,i]
    filename=paste0(RegFolderName,"Reg_cluster",i,".bed")
    write.table(outdf,filename,sep="\t",col.names = F,row.names = F,quote = FALSE)
  }


  return(a=c(barcordFileName=bfilename,RegFolderName=RegFolderName))
}

