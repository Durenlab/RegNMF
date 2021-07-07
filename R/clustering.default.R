clustering.default<-function(H,cell_num_smooth=sqrt(dim(H)[2])){
  HH=t(t(H)/colSums(H));
  S=matrix(ncol = 3);
  S=na.omit(S);
  S_tmp=Fold_sim_based_clustering2.default(HH,cell_num_smooth);
  #  S[,1]=S_tmp[1,];
  #  S_tmp=Fold_sim_based_clustering2.default(HH,50);
  #  S[,2]=S_tmp[1,];
  #  S_tmp=Fold_sim_based_clustering2.default(HH,150);
  #  S[,3]=S_tmp[1,];
  #cluster level W3
  #  H_norm=W123H[4]/sqrt(rowSums(W123H$H*W123H$H));
  #  W3_norm=W123H[3]/sqrt(rowSums(W123H$H*W123H$H));
  #  H_w=matrix(ncol = length(unique(S[,1])));
  #  H_w=na.omit(H_w);
  #  for (i in 1:length(unique(S[,1]))) {
  #    H_w[,i]=rowMeans(H_norm[,S[,1]==i]);
  #  }
  #  W3_cluster=W3_norm%*%H_w;
  T=Rtsne::Rtsne(t(HH));
  tsne_plot <- data.frame(x=T$Y[,1],y=T$Y[,2],col=S_tmp[1,]);
  ans=ggplot2::ggplot(tsne_plot)+ggplot2::geom_point(ggplot2::aes(x=x,y=y,colour=factor(col)));

  return(list(HH=HH,
              S=S_tmp,
              T=T,
              plot=ans))
}
