Fold_RE_TG_multiAdjust.default<- function (E2,O2,Symbol_location,Peak_location){
  E2=as.matrix(E2);
  O2=as.matrix(O2);
  print(typeof(E2));
  print(typeof(O2));
  print(typeof(Symbol_location));
  print(typeof(Peak_location));
  P_1=Fold_RE_TG_MultiAdjustCore(E2, O2, Symbol_location, Peak_location);

  #T_List=Fold_RE_TG_MultiAdjustCore(P_1E2, O2, Symbol_location, Peak_location);
  #T_1 <- sparseMatrix(T_List[1],T_List[2],T_List[3]);
  #T_1[is.na(T_1)] <- 0;
  return(P_1);
}


#Fold_RE_TG_MultiAdjustCore(P_1, E11, O11, Symbol_location, Peak_location);
