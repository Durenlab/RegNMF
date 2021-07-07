
#Reg_symbol_name=c()
#peak_location=c()
#for (i in 1:1000000) {
#  a=paste0("Symbol",i);
#  Reg_symbol_name[i]=a;
#  a=paste0("chr?:",i,"-",i+100);
#  peak_location[i]=a;
#}
#c=rep(1:12,300)
#d=matrix(c(c[1:3012],1:3012),3012,2)
#d=data.frame(d)
SplitGroup.default<-function(foldername,barcord,W3,H,Reg_symbol_name,Reg_peak_name,cluster){
  clustern=length(unique(cluster))
  barcord_cluster=data.frame(barcord=barcord,cluster=cluster)
  bfilename=paste0(foldername,"barcord_cluster.bed")
  write.table(barcord_cluster,filename,col.names = F,row.names = F,quote = FALSE)
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
  if(!dir.exists(name))dir.create(name,recursive = TRUE )
  dir.create("",recursive = TRUE )
  for (i in 1:clustern) {
    topk=order(W3_cluster[,i])[1:10000]
    outdf=df[topk,]
    outdf$Reg=W3_cluster[topk,i]
    filename=paste0(RegFolderName,"Reg_cluster",i,".bed")
    write.table(outdf,filename,col.names = F,row.names = F,quote = FALSE)
  }


  return(a=c(barcordFileName=bfilename,RegFolderName=RegFolderName))
}

