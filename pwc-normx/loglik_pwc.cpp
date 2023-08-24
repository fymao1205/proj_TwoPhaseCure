//#include <Rcpp.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include <Rmath.h>
#include <iterator>
#include <iostream>
#include <vector>
#include <algorithm>

#include "commonf.h"

using namespace Rcpp;
using namespace arma;


// Obtain environment containing function
//Rcpp::Environment package_env("package:mnormt"); 

// Make function callable from C++
//Rcpp::Function bipmvnrm = package_env["pmnorm"];    

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

//[[Rcpp::export()]]
NumericMatrix gaussleg(int n, double x1, double x2){
  
  double EPS = 3e-14;
  double m = 0.5*(n+1);
  //double xm, xl;
  double xm = 0.5 * (x2 + x1);
  double xl = 0.5 * (x2 - x1);
  
  //NumericVector x(n), w(n);
  //x = rep(0.0, n);
  //w = rep(0.0, n);
  NumericVector x=no_init(n);
  NumericVector w=no_init(n);
  double z;
  //double tmp;
  for (int i=0; i<m; i++) {
    
    double tmp = M_PI * ( 1.0*(i+1) - 0.25)/(n*1.0 + 0.5);
    z = cos(tmp);
    
    double tol = 9999;
    //double p1, p2, p3;
    double pp;
    //double z1;
    while (tol > EPS) {
      double p1 = 1.0;
      double p2 = 0.0;
      for (int j=0; j<n; j++) {
        double p3 = p2;
        p2 = p1;
        p1 = ((2 * j + 1) * z * p2 - j  * p3)/(j+1);
      }
      
      pp = n * (z * p1 - p2)/((z-1)*(1+z));
      double z1 = z;
      z = z1 - p1/pp;
      tol = abs(z - z1);
      //z -= p1/pp;
      //tol= abs(p1/pp);
    }
    
    int s=n -1 - i;
    x[i] = xm - xl * z;
    x[s] = xm + xl * z;
    w[i] = (2 * xl)/((1-z)*(1+z) * pp * pp);
    w[s] = w[i];
  }
  
  NumericMatrix res=no_init_matrix(n, 2);
  res(_,0) = x;
  res(_,1) =w;
  
  //double tmp1 = cos(tmp);
  //Rcout << "tmp is " << tmp << std::endl;
  //Rcout << "tmp1 is " << tmp1 << std::endl;
  return(res);
}


double logL_cal_pwc_norm(int del, double t, double x1, double eta1, double beta1, 
                         double muZ, double muT, double muX1, double sigma2, NumericVector alp, 
                    //NumericVector eta, NumericVector alp, NumericVector gam, NumericVector beta, 
                    NumericVector brks){
  double expcovZ=exp(x1*eta1+muZ); //exp(eta[0]+x1*eta[1]+x2*eta[2]); 
  //double hf = alp[0]*alp[1]*pow(alp[0]*t, alp[1]-1)*exp(beta[0]*x1+beta[1]*x2); 
  //double Hf = pow(alp[0]*t, alp[1])*exp(beta[0]*x1+beta[1]*x2);
  NumericVector levs = alp*exp(beta1*x1+muT); 
  double hf = hpwc_double(t, brks, levs, 0); 
  double Hf = Hpwc_double(t, brks, levs, 0);
  double Sf = exp(-Hf); 
  double piZ = expcovZ/(1+expcovZ);
  
  double logL1 = del*(log(hf)-Hf+log(piZ)) + (1-del)*log(1-(1-Sf)*piZ);
  //double mu1=gam[0]+gam[1]*x2;
  double logL2 = R::dnorm4(x1, muX1, sigma2, 1);//x1*(gam[0]+gam[1]*x2) - log(1+exp(gam[0]+gam[1]*x2)); //x1*(covX-log(1+expcovX))+(1-x1)*(-log(1+expcovX));
  
  double res=logL1+logL2;
  
  return res;
}


double logL_obs_pwc_norm(int r, int del, double t, double x1, double eta1, double beta1, 
                         double muZ, double muT, double muX1, double sigma2, NumericVector alp, 
                    //NumericVector eta, NumericVector alp, NumericVector gam, NumericVector beta, 
                    NumericVector brks){
  double res=0.0;
  if(r==1){
    res = logL_cal_pwc_norm(del, t, x1, eta1, beta1, muZ, muT, muX1, sigma2, alp, brks);
  }else{
    
    NumericMatrix gau_a = gaussleg(20, -1, 1.0);
    NumericVector u = gau_a(_,0);
    NumericVector w = gau_a(_,1);
    
    double tmp=0.0;
    for(int i=0; i<20; i++){
      tmp += (w[i]*(exp(logL_cal_pwc_norm(del, t, u[i], eta1, beta1, muZ, muT, muX1, sigma2, alp, brks))+ exp(logL_cal_pwc_norm(del, t, 1/u[i], eta1, beta1, muZ, muT, muX1, sigma2, alp, brks))*pow(u[i],-2) ));
    }
    
    //double res0 = logL_cal_pwc_norm(del, t, 0.0, x2, eta, alp, gam, beta, brks);
    //double res1 = logL_cal_pwc(del, t, 1.0, x2, eta, alp, gam, beta, brks);
    res = log(tmp);//log(exp(res0)+exp(res1));
  }
  
  return res;
}

// [[Rcpp::export()]]
NumericVector logLn_obs_pwc_norm(IntegerVector r, IntegerVector del, NumericVector t, NumericVector x1, 
                                 double beta1, double eta1, NumericVector muZ, NumericVector muT, 
                                 double muX1, double sigma2, NumericVector alp,
                                 NumericVector brks){
  
  NumericVector res=no_init(del.size());
  
  for(int i=0; i<del.size(); i++){
    res[i] = logL_obs_pwc_norm(r[i], del[i], t[i], x1[i], eta1, beta1, muZ[i], muT[i], muX1, sigma2, alp, brks); //wgt[i]*(temp);
  }
  
  return res;
}


double logL_cal_pwc_norm_H0(int del, double t, double muZ, double muT, NumericVector alp,
                            NumericVector brks){
  double expcovZ=exp(muZ); //exp(eta[0]);//
  //double hf = alp[0]*alp[1]*pow(alp[0]*t, alp[1]-1)*exp(beta[0]*x1+beta[1]*x2); 
  //double Hf = pow(alp[0]*t, alp[1])*exp(beta[0]*x1+beta[1]*x2);
  NumericVector levs = alp*exp(muT); 
  double hf = hpwc_double(t, brks, levs, 0); 
  double Hf = Hpwc_double(t, brks, levs, 0);
  double Sf = exp(-Hf); 
  double piZ = expcovZ/(1+expcovZ);
  
  double logL1 = del*(log(hf)-Hf+log(piZ)) + (1-del)*log(1-(1-Sf)*piZ);
  //double logL2 = x1*(gam[0]+gam[1]*x2) - log(1+exp(gam[0]+gam[1]*x2)); //x1*(covX-log(1+expcovX))+(1-x1)*(-log(1+expcovX));
  
  double res=logL1;
  
  return res;
}

// [[Rcpp::export()]]
NumericVector logLn_obs_pwc_norm_H0(IntegerVector del, NumericVector t, NumericVector muZ, NumericVector muT, 
                                    NumericVector alp, NumericVector brks){
  
  NumericVector res=no_init(del.size());
  
  for(int i=0; i<del.size(); i++){
    res[i] = logL_cal_pwc_norm_H0(del[i], t[i], muZ[i], muT[i], alp, brks); //wgt[i]*(temp);
  }
  
  return res;
}

// [[Rcpp::export()]]
NumericVector vHpwc_H0(NumericVector t, NumericVector muT, 
                       NumericVector alp, NumericVector brks){
  
  NumericVector res=no_init(t.size());
  for(int i=0; i<t.size(); i++){
    
    //double piX = expit.f(eta[0]+eta[1]*x2[i]);
    NumericVector levs = alp*exp(muT[i]);
    res[i] = Hpwc_double(t[i], brks, levs, 0); 
  }
  
  return res;
  
}

// [[Rcpp::export()]]
NumericVector vppwc_H0(NumericVector t, NumericVector muT, 
                       NumericVector alp, NumericVector brks){
  
  NumericVector res=no_init(t.size());
  for(int i=0; i<t.size(); i++){
    
    //double piX = expit.f(eta[0]+eta[1]*x2[i]);
    NumericVector levs = alp*exp(muT[i]);
    res[i] = ppwc_double(t[i], brks, levs, 1, 0); 
  }
  
  return res;
  
}

