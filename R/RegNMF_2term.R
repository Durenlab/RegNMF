RegNMF_2term <- function(A,
                         E,
                         feature_cut_perc=0.01,
                         cell_num_smooth=sqrt(dim(E)[2]),
                         core=8){

  K=100;
  #####Normlize and filter
  A1=normlize_scalling.default(A);
  E1=normlize_scalling.default(E);
  numCut=feature_cut_perc*dim(E1)[2];
  a=rowSums(as.matrix(E1)>0);
  a1=rowSums(as.matrix(A1)>0);
  E11=E1[a>numCut,];
  A11=A1[a1>numCut,];

###Making W10,W20,H0, which are used in initialize lambda
  W10=matrix(runif(dim(A11)[1]*K,min = 0,max=1),dim(A11)[1],K);
  W20=matrix(runif(dim(E11)[1]*K,min = 0,max=1),dim(E11)[1],K);
  H0=matrix(runif(dim(A11)[2]*K),K,dim(A11)[2]);
  wh1=CNmf(as.matrix(A11),100,100,W10,H0,core);
  wh2=CNmf(as.matrix(E11),100,100,W20,H0,core);
  wh1$W=t(t(wh1$W)*sqrt(rowSums(wh1$H*wh1$H)))
  wh2$W=t(t(wh2$W)*sqrt(rowSums(wh2$H*wh2$H)))
  beta=1
  arfa=1
  lambda=beta*mean(t(wh1$W)%*%A11)/mean(t(wh2$W)%*%E11);

  W10=matrix(runif(dim(A11)[1]*K,min = 0,max=1),dim(A11)[1],K);
  W20=matrix(runif(dim(E11)[1]*K,min = 0,max=1),dim(E11)[1],K);
  H0=matrix(runif(dim(A11)[2]*K),K,dim(A11)[2]);
  W12H=CPPNMF_cluster_joint_cross_domain_2term(as.matrix(A11),as.matrix(E11),K,300,lambda,W10,W20,H0,core);

  return(W12H);

}
