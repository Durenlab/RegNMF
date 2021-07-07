tfidf_trans.default <-function(O){
  O=as.matrix((O>0));
  tf1=t(t(O)/log(1+colSums(O)));
  idf=log(1+dim(O)[2]/(1+rowSums(O>0)));
  O1=tf1*idf;
  O1[is.na(O1)] <- 0;
  return(as(O1,"sparseMatrix"))
}
