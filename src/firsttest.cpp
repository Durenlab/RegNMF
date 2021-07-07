#include "../inst/include/util.h"

//[[Rcpp::export]]
Rcpp::List CNmf(Eigen::Map<Eigen::MatrixXd> V, int K, int maxiter, Eigen::Map<Eigen::MatrixXd> W0,Eigen::Map<Eigen::MatrixXd> H0,int core){
  Eigen::setNbThreads(core);
  int n = Eigen::nbThreads();
  Rprintf("Core=%d\n",n);

  Eigen::MatrixXd W(W0.rows(),W0.cols());
  Eigen::MatrixXd H(H0.rows(),H0.cols());


  for(int iter=0;iter<maxiter;iter++){
    Rcpp::checkUserInterrupt();
    Eigen::MatrixXd W0TV=W0.transpose()*V;
    Eigen::MatrixXd W0TW0H=W0.transpose()*W0*H0;
    H=ifzeroMCPP((H0.array()*(W0TV.array()/(W0TW0H+epsMCpp(W0TV)).array())).matrix());

    Eigen::MatrixXd VHT=V*H.transpose();
    Eigen::MatrixXd W0HHT=W0*H*H.transpose();

    W=ifzeroMCPP((W0.array()*(VHT.array()/(W0HHT+epsMCpp(VHT)).array())).matrix());

    W0=W;
    H0=H;
    Rprintf("%d\n",iter);
  }

  return Rcpp::List::create(Named("W") = wrap(W),
                            Named("H") = wrap(H));
}



// [[Rcpp::export]]
Rcpp::List CPPNMF_cluster_joint_cross_domain_try(Eigen::Map<Eigen::MatrixXd> PeakO,
                                                 Eigen::Map<Eigen::MatrixXd> X,
                                                 Eigen::Map<Eigen::MatrixXd> Reg,
                                                 int K,
                                                 int maxiter,
                                                 double lambda1,
                                                 double lambda2,
                                                 Eigen::Map<Eigen::MatrixXd> W10,
                                                 Eigen::Map<Eigen::MatrixXd> W20,
                                                 Eigen::Map<Eigen::MatrixXd> W30,
                                                 Eigen::Map<Eigen::MatrixXd> H0,
                                                 NumericVector c1,
                                                 NumericVector c2,
                                                 NumericVector Reg_w,
                                                 int core){


  Eigen::setNbThreads(core);
  int n = Eigen::nbThreads();
  Rprintf("Core=%d\n",n);
  double tolx=1e-4;
  double tolfun=1e-6;
  double sqrteps = sqrt(DBL_EPSILON);
  double dnorm,dnorm0;
  Eigen::MatrixXd HH, numerO,numerX, numerR ,W1, W2, W3, H, Tmp1, Tmp2,Tmp3,HT;





  //dnorm0=(PeakO-(W10*H0)).squaredNorm()+(lambda1*(X-(W20*H0)).squaredNorm())+(lambda2*(Reg-(W30*H0)).squaredNorm());

  for(int iter=0;iter<maxiter;iter++){
    Rcpp::checkUserInterrupt();
    Eigen::MatrixXd W10T=W10.transpose();
    Eigen::MatrixXd W20T=W20.transpose();
    Eigen::MatrixXd W30T=W30.transpose();
    Eigen::MatrixXd H0T=H0.transpose();
    numerO=(W10T*PeakO)+(lambda1*(W20T*X)) +(lambda2*(W30T*Reg));

    H=ifzeroMCPP((H0.array()*(numerO.array()/(((W10T*W10)+(lambda1*W20T*W20)+(lambda2*W30T*W30))*H0+epsMCpp(numerO)).array())).matrix());

    HT=H.transpose();
    HH=H*HT;


    numerO=PeakO*HT;
    W1=ifzeroMCPP((W10.array()*(numerO.array()/(W10*HH+epsMCpp(numerO)).array())).matrix());


    numerX=X*HT;

    W2=ifzeroMCPP((W20.array()*(numerX.array()/(W20*HH+epsMCpp(numerX)).array())).matrix());


    Tmp1=chooesVinMCPP(W1,c2,0);
    Tmp2=chooesVinMCPP(W2,c1,0);
    numerR=Tmp1+Tmp2;
    numerR=CppoperationMA_demo(numerR,Reg_w,2);
    Tmp1=Reg*HT;
    numerR+=Tmp1;

    W3=ifzeroMCPP((W30.array()*(numerR.array()/(W30*HH+W30+epsMCpp(numerR)).array())).matrix());

    if(iter>300){
      if((dnorm0-dnorm)<=tolfun||(dnorm0-dnorm)<=dnorm0){

        Rprintf("dnorm0-dnorm %f is small", dnorm0-dnorm);
        break;
      }
      else if(iter==maxiter)
        break;
    }

    W10 = W1;
    H0 = H;
    W20 = W2;
    W30 = W3;

    if(iter%20==0){
    /*
      Tmp1=-(W1*H);
      Tmp1+=PeakO;
      dnorm=Tmp1.squaredNorm();
      Tmp1=-(W2*H);
      Tmp1+=X;
      dnorm+=lambda1*Tmp1.squaredNorm();
      Tmp1=-(W3*H);
      Tmp1+=Reg;
      dnorm+=lambda2*Tmp1.squaredNorm();
      */
      Rprintf("iteration %d\n",iter);
    }


}

return Rcpp::List::create(Named("W1") = wrap(W1),
                          Named("W2") = wrap(W2),
                          Named("W3") = wrap(W3),
                          Named("H") = wrap(H));


}



// [[Rcpp::export]]
NumericMatrix Fold_RE_TG_MultiAdjustCore(NumericMatrix E2,
                                         NumericMatrix O2,
                                         NumericMatrix Symbol_location,
                                         NumericMatrix Peak_location){
  LogicalVector location1,location2,preid;
  NumericMatrix E2sqrt (E2.rows(),E2.cols());
  NumericMatrix P_1 (E2.rows(),O2.rows());
  IntegerVector id,rowid,colid;

  for(int i=0;i<(E2.length());i++){
    if(E2[i]!=0){E2sqrt[i]=pow(E2[i],2);}
    }
  Function w("which");
  Rprintf("check2\n");
  for(int i=0;i<O2.rows();i++){
    Rcpp::checkUserInterrupt();
    location1=Symbol_location(_,0)==Peak_location(i,0);
    location2=abs(Symbol_location(_,1)-Peak_location(i,1))<1000000;
    preid=location1&location2;
    id=w(preid==TRUE);
    id=id-1;
    int tsum=0;
    for(int j=0;j<O2.cols();j++){
      if(O2(i,j)>0) tsum++;
    }
    for(int j=0;j<id.length();j++){
      int cnt0=0;
      int cnt1=0;

      for(int k=0;k<O2.cols();k++){
        if(E2(id(j),k)!=0){
          if(O2(i,k)>0){
            cnt1++;
          }

          else if (O2(i,k)==0) {
            cnt0++;
          }
        }

      }

      NumericVector set0 (cnt0);
      NumericVector set0sqrt (cnt0);
      NumericVector set1 (cnt1);
      NumericVector set1sqrt (cnt1);

      for(int k=0;k<O2.cols();k++){
      if(E2(id(j),k)!=0){
        if(O2(i,k)>0){
          cnt1--;
          set1(cnt1)=E2(id(j),k);
          set1sqrt(cnt1)=E2sqrt(id(j),k);
          }

        else if (O2(i,k)==0) {
          cnt0--;
          set0(cnt0)=E2(id(j),k);
          set0sqrt(cnt0)=E2sqrt(id(j),k);

          }
      }

      }
      P_1(id(j),i)=Cttest(set1,set0,set1sqrt,set0sqrt,tsum,O2.cols()-tsum);
    }
     if(i%5000==0) {Rprintf("%d\n",i);
     }

  }
  return wrap(P_1);
}

