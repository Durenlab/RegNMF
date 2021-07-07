#ifndef UTIL_H
#define UTIL_H

#include <Rcpp.h>
#include<RcppEigen.h>
#include<math.h>
#include<float.h>
#include<stdlib.h>
#include<vector>

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::plugins(cpp11)]]
//[[Rcpp::plugins(openmp)]]
using namespace Rcpp;




double eps(double a);
Eigen::MatrixXd CppoperationMA_demo(Eigen::MatrixXd M, NumericVector A,int type);
Eigen::MatrixXd ifzeroMCPP(Eigen::MatrixXd M);
Eigen::MatrixXd chooesVinMCPP(Eigen::MatrixXd M, NumericVector A ,int type);
Eigen::MatrixXd epsMCpp(Eigen::MatrixXd M);
double Cvar(NumericVector AA, double Mean, int length);
double Cttest(NumericVector A, NumericVector B,NumericVector AA,NumericVector BB,int lengthA, int lengthB);
NumericMatrix Cjaccard(NumericMatrix MM);

#endif
