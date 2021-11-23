normlize_scalling.default <- function (A){
 listCols<-function(m){
    #converts a sparse Matrix into a list of its columns
    res<-split(m@x, findInterval(seq_len(Matrix::nnzero(m)), m@p, left.open=TRUE))
    names(res)<-colnames(m)
    res
}
  list0=listCols(A)
b=c()
for (i in 1:length(list0)){
b=c(b,median(list0[i][[1]],na.rm=TRUE))
}
return(Matrix::t(Matrix::t(E)/b*mean(b)))
}
