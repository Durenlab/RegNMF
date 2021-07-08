defaultpar_CoupledNMF.default <- function(PeakO,W1,X,W2,Reg,W3,beta=1,arfa=1){
  lambda1=beta*mean(t(W1)%*%PeakO)/mean(t(W2)%*%X);
  lambda2=arfa*mean(t(W1)%*%PeakO)/mean(t(W3)%*%Reg);
  return(list(lambda1,lambda2));
}