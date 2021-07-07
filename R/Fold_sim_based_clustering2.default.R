Fold_sim_based_clustering2.default <- function(O,k=20){
  sim=cosine(O);
  KK=colSums(sim);
  twom=sum(KK);
  sim_norm=sim-(as.matrix(KK) %*% (t(as.matrix(KK/twom))))
  f=matrix(0,dim(sim_norm)[1],dim(sim_norm)[2]);
  for (i in (1:dim(sim_norm)[2])) {
    f[,i]=order(sim_norm[,i],decreasing = TRUE);
  }
  print(1);

  KNN=matrix(0,dim(sim_norm)[2],dim(sim_norm)[2]);

  for (i in (1:dim(sim_norm)[2])) {
    KNN[f[2:(k+1),i],i]=1;
  }
  KNN2=t(KNN)%*%KNN;
  A=Cjaccard(KNN2);
  A[is.na(A)] <- 0;
  KK=colSums(as.matrix(A));
  twom=sum(KK);
  B=A-((as.matrix(KK)%*%t(as.matrix(KK)))/twom);
  a=netZooR::alpaca.genlouvain(B);
  return(a);
}
