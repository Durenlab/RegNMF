#include "../inst/include/util.h"


// [[Rcpp::export]]
double eps(double a){
  int i=1;
  int b =(int)abs(a);
  double c=b;
  double epsilon=DBL_EPSILON;
  while((c/2)>1){
    i++;
    c/=2;
    epsilon*=2;
  }
  return epsilon;
}



Eigen::MatrixXd CppoperationMA_demo(Eigen::MatrixXd M, NumericVector A,int type){

  Eigen::MatrixXd ans(M.rows(),M.cols());

  switch (type)
  {
  case 0:
    for (int i = 0; i < M.rows(); i++){
      for (int j = 0; j < M.cols(); j++){
        ans(i,j)=A[i]+M(i,j);
      }
    }
    break;
  case 1:
    for (int i = 0; i < M.rows(); i++){
      for (int j = 0; j < M.cols(); j++){
        ans(i,j)=M(i,j)-A[i];
      }
    }
    break;
  case 2:
    for (int i = 0; i < M.rows(); i++){
      for (int j = 0; j < M.cols(); j++){
        ans(i,j)=A[i]*M(i,j);
      }
    }
    break;
  case 3:
    for (int i = 0; i < M.rows(); i++){
      for (int j = 0; j < M.cols(); j++){
        ans(i,j)=M(i,j)/A[i];
      }
    }
    break;

  }
  return ans;

}



Eigen::MatrixXd ifzeroMCPP(Eigen::MatrixXd M){
  Eigen::MatrixXd ans(M.rows(),M.cols());
  for (int i = 0; i < M.rows(); i++)
  {
    for (int j = 0; j < M.cols(); j++)
    {
      if(M(i,j)>0)    ans(i,j)=M(i,j);
      else    ans(i,j)=0;
    }

  }
  return ans;
}


Eigen::MatrixXd chooesVinMCPP(Eigen::MatrixXd M, NumericVector A ,int type){
  if(type==0){
    Eigen::MatrixXd ans(A.length(),M.cols());
    for(int i=0;i<A.length();i++){
      for(int j=0;j<ans.cols();j++){
        ans(i,j)=M((A[i]-1),j);
      }
    }
    return ans;
  }

  if(type==1){
    Eigen::MatrixXd ans(M.rows(),A.length());
    for(int i=0;i<ans.rows();i++){
      for(int j=0;j<A.length();j++){
        ans(i,j)=M(i,(A[j]-1));
      }
    }
    return ans;
  }

}

Eigen::MatrixXd epsMCpp(Eigen::MatrixXd M){
  Eigen::MatrixXd ans(M.rows(),M.cols());
  for (int i = 0; i < M.rows(); i++)
  {
    for (int j = 0; j < M.cols(); j++)
    {
      ans(i,j)=eps(M(i,j));
    }

  }
  return ans;
}


double Cvar(NumericVector AA, double Mean, int length){
  double Sum=sum(AA)/(length-1);
  Sum-=pow(Mean,2.0)*length/(length-1);
  return Sum;
}




double Cttest(NumericVector A, NumericVector B,NumericVector AA,NumericVector BB,int lengthA, int lengthB){
  double meanA=sum(A)/lengthA;
  double meanB=sum(B)/lengthB;
  double t=(meanA-meanB)/sqrt((Cvar(AA,meanA,lengthA)/lengthA)+(Cvar(BB,meanB,lengthB)/lengthB));
  return t;
}

// [[Rcpp::export]]
NumericMatrix Cjaccard(NumericMatrix MM){
  NumericMatrix a (MM.rows(),MM.cols());

  for(int i=0;i<MM.rows();i++){
    for(int j=0;j<MM.cols();j++){
      a(i,j)=MM(i,j)/(MM(i,i)+MM(j,j)-MM(i,j));
    }
  }
  return a;
}
