normlize_scalling.default <- function (A){
  B=as.matrix(A);
  B[B==0]=NaN;
  b=apply(B,2,median,na.rm=TRUE);
  return(as(t(t(as.matrix(A))/b*mean(b)),"sparseMatrix"));
}
