RegNMF.default <- function(E,
                   O,
                   Symbol,
                   PeakName,
                   Symbol_location,
                   Peak_location,
                   feature_cut_perc=0.01,
                   corr_cut_k=100000,
                   cell_num_smooth=sqrt(dim(E)[2]),
                   core=8){

  K=100;
  #####normlize and filter
  E1=normlize_scalling.default(E);
  O1=normlize_scalling.default(O);
  numCut=feature_cut_perc*dim(E1)[2];
  a=rowSums(as.matrix(E1)>0);
  a1=rowSums(as.matrix(O1)>0);
  E11=E1[a>numCut,];
  O11=O1[a1>numCut,];
  Symbol_location=Symbol_location[a>numCut,];
  Peak_location=Peak_location[a1>numCut,];
  Symbol=Symbol[a>numCut];
  PeakName=PeakName[a1>numCut];

  R=Fold_RE_TG_multiAdjust.default(E11,O11,Symbol_location,Peak_location)
  R=as(R,"sparseMatrix");
  mk=sort(R[R>0],decreasing = TRUE);
  corr_cut_k=min(length(mk),corr_cut_k)
  mk=mk[1:corr_cut_k];

  corr_cut=mk[length(mk)];
  a=colSums((as.matrix(R)>=corr_cut));
  O11=O11[a>0,];
  a2=rowSums((as.matrix(R)>=corr_cut));
  E11=E11[a2>0,];
  R=R[a2>0,a>0];
  Symbol_location=Symbol_location[a2>0,];
  Peak_location=Peak_location[a>0,]
  Symbol=Symbol[a2>0];
  PeakName=PeakName[a>0];
  c=which(((as.matrix(R)>=corr_cut)==1),arr.ind = TRUE);
  Reg=O11[c[,2],]+E11[c[,1],];
  Reg_gene_location=Symbol_location[c[,1],];
  Reg_peak_location=Peak_location[c[,2],];
  Reg_dis=abs(Reg_gene_location[,2]-Reg_peak_location[,2]);
  d0=2*10^5;
  Reg_w=exp(-Reg_dis/d0);
  Reg_adj=Reg*Reg_w;
  O12=tfidf_trans.default(O11);
  Reg_info=cbind(c[,1],c[,2]);
  Reg_info=cbind(Reg_info,Reg_dis);

  W10=matrix(runif(dim(O12)[1]*K,min = 0,max=1),dim(O12)[1],K);
  W20=matrix(runif(dim(E11)[1]*K,min = 0,max=1),dim(E11)[1],K);
  W30=matrix(runif(dim(Reg_adj)[1]*K,min = 0,max=1),dim(Reg_adj)[1],K);
  H0=matrix(runif(dim(O12)[2]*K),K,dim(O12)[2]);
  wh1=CNmf(as.matrix(O12),100,100,W10,H0,core);
  wh2=CNmf(as.matrix(E11),100,100,W20,H0,core);
  wh3=CNmf(as.matrix(Reg_adj),100,100,W30,H0,core);
  wh1$W=t(t(wh1$W)*sqrt(rowSums(wh1$H*wh1$H)))
  wh2$W=t(t(wh2$W)*sqrt(rowSums(wh2$H*wh2$H)))
  wh3$W=t(t(wh3$W)*sqrt(rowSums(wh3$H*wh3$H)))
  lambda=defaultpar_CoupledNMF.default(O12,wh1$W,E11,wh2$W,Reg_adj,wh3$W,2,1);

  W10=matrix(runif(dim(O12)[1]*K,min = 0,max=1),dim(O12)[1],K);
  W20=matrix(runif(dim(E11)[1]*K,min = 0,max=1),dim(E11)[1],K);
  W30=matrix(runif(dim(Reg_adj)[1]*K,min = 0,max=1),dim(Reg_adj)[1],K);
  H0=matrix(runif(dim(O12)[2]*K),K,dim(O12)[2]);
  W123H=CPPNMF_cluster_joint_cross_domain_try(as.matrix(O12),as.matrix(E11),as.matrix(Reg_adj),K,300,lambda[[1]],lambda[[2]],W10,W20,W30,H0,c[,1],c[,2],Reg_w,core);


  W123H[["Reg_gene_name"]]=Symbol[c[,1]];
  W123H[["Reg_peak_name"]]=PeakName[c[,2]];
  W123H[["Gene_name"]]=Symbol;
  W123H[["peak_name"]]=PeakName;



  return(W123H);


}
