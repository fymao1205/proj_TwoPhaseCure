// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
//#include <omp.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include <Rmath.h>
#include <iterator>
#include <iostream>
#include <vector>
#include <algorithm>

using namespace Rcpp;
using namespace arma;


//[[Rcpp::export()]]
double hpwc_double(double x, NumericVector cuts, NumericVector levels, int logInd)
{
  
  NumericVector cut = sort_unique(cuts);
  int p = levels.size();
  double y;
  cut.push_front(0);
  cut.push_back(R_PosInf);
  
  if((cut[0] <= x) & (x <= cut[1])){
    y = levels[0];
  }
  if (p > 1.5) {
    for (int i=1; i<p; i++) {
      //y[(cut[i] <= x) & (x <= cut[i + 1])] = levels[i];
      if((cut[i] <= x) & (x <= cut[i+1])){
        y = levels[i];
      }
    }
  }
  if (logInd)
    y = log(y);
  return(y);
}

// [[Rcpp::export]]
double Hpwc_double(double x, NumericVector cuts, NumericVector levels, int logInd)
{
  NumericVector cut = sort_unique(cuts);
  int p = levels.size();
  //NumericVector y(x.size());
  double y;
  cut.push_front(0);
  cut.push_back(R_PosInf);
  NumericVector z;
  bool which = (cut[0] <= x) & (x <= cut[0 + 1]);
  if (which) {
    //y[which] = x[which];
    y = x*levels[0];
  }
  double wt = levels[0] * cut[1];
  if (p > 1.5) {
    for (int i = 1; i<p; i++) {
      which = (cut[i] <= x) & (x <= cut[i + 1]);
      if (which) {
        //NumericVector xwhich= x[which];
        double tmpx = wt + levels[i] * (x - cut[i]);
        y = tmpx;
      }
      wt = wt + levels[i] * (cut[i + 1] - cut[i]);
    }
  }
  
  //Rcout << "which is " << which << std::endl;
  //Rcout << "wt is " << wt << std::endl;
  //Rcout << "cut is " << cut << std::endl;
  //Rcout << "p is " << p << std::endl;
  //Rcout << "y is " << y << std::endl;
  
  if (logInd)
    y = log(y);
  return(y);
}

//[[Rcpp::export()]]
double ppwc_double(double q, NumericVector cuts, NumericVector levels, int lower, int logInd)
{
  double y;
  if (cuts[0]==0) {
    y = R::pexp(q, 1/levels[0], 0.0, 0.0);
  }else{
    //NumericVector qq(1);
    //qq[0] = q;
    y = Hpwc_double(q, cuts, levels, 0.0);
    if (logInd) {
      if (lower) {
        y = log(-(exp(-y)-1));
      }
      else {
        y = -y;
      }
    }
    else {
      if (lower) {
        y = -(exp(-y)-1);
      }
      else {
        y = exp(-y);
      }
    }
  }
  return(y);
}

//[[Rcpp::export()]]
NumericVector hpc(NumericVector x, NumericVector levels, NumericVector cuts, int logf)
{
  
  NumericVector cut = sort_unique(cuts);
  int p = levels.size();
  NumericVector y(x.size());
  cut.push_front(0);
  cut.push_back(R_PosInf);
  
  y[(cut[0] <= x) & (x < cut[1])] = levels[0];
  if (p > 1.5) {
    for (int i=1; i<p; i++) {
      y[(cut[i] <= x) & (x < cut[i + 1])] = levels[i];
    }
  }
  if (logf)
    y = log(y);
  return(y);
}

//[[Rcpp::export()]]
NumericVector Hpc(NumericVector x,  NumericVector levels, NumericVector cuts, int logf)
{
  NumericVector cut = sort_unique(cuts);
  int p = levels.size();
  NumericVector y(x.size());
  cut.push_front(0);
  cut.push_back(R_PosInf);
  NumericVector z;
  LogicalVector who = (cut[0] <= x) & (x < cut[0 + 1]);
  if (sum(who)) {
    y[who] = x[who];
    y = y*levels[0];
  }
  double su = levels[0] * cut[1];
  if (p > 1.5) {
    for (int i = 1; i<p; i++) {
      who = (cut[i] <= x) & (x < cut[i + 1]);
      if (sum(who)) {
        NumericVector xwho= x[who];
        NumericVector tmpx = su + levels[i] * (xwho - cut[i]);
        y[who] = tmpx;
      }
      su = su + levels[i] * (cut[i + 1] - cut[i]);
    }
  }
  if (logf)
    y = log(y);
  return(y);
}


//[[Rcpp::export()]]
NumericVector vppc(NumericVector q, NumericVector levels,  NumericVector cuts, int lower, int logf)
{
  NumericVector y(1);
  if (cuts[0]==0) {
    y = pexp(q, levels[0], lower, logf);
  }else{
    y = Hpc(q,  levels, cuts, 0.0);
    if (logf) {
      if (lower) {
        y = log(-(exp(-y)-1));
      }
      else {
        y = -y;
      }
    }
    else {
      if (lower) {
        y = -(exp(-y)-1);
      }
      else {
        y = exp(-y);
      }
    }
  }
  return(y);
}

//[[Rcpp::export()]]
NumericVector hpwc(NumericVector x, NumericVector cuts, NumericVector levels, int logInd)
{
  
  NumericVector cut = sort_unique(cuts);
  int p = levels.size();
  NumericVector y(x.size());
  cut.push_front(0);
  cut.push_back(R_PosInf);
  
  y[(cut[0] <= x) & (x <= cut[1])] = levels[0];
  if (p > 1.5) {
    for (int i=1; i<p; i++) {
      y[(cut[i] <= x) & (x <= cut[i + 1])] = levels[i];
    }
  }
  if (logInd)
    y = log(y);
  return(y);
}

//[[Rcpp::export()]]
NumericVector Hpwc(NumericVector x, NumericVector cuts, NumericVector levels, int logInd)
{
  NumericVector cut = sort_unique(cuts);
  int p = levels.size();
  NumericVector y(x.size());
  cut.push_front(0);
  cut.push_back(R_PosInf);
  NumericVector z;
  LogicalVector which = (cut[0] <= x) & (x <= cut[0 + 1]);
  if (sum(which)) {
    y[which] = x[which];
    y = y*levels[0];
  }
  double wt = levels[0] * cut[1];
  if (p > 1.5) {
    for (int i = 1; i<p; i++) {
      which = (cut[i] <= x) & (x <= cut[i + 1]);
      if (sum(which)) {
        NumericVector xwhich= x[which];
        NumericVector tmpx = wt + levels[i] * (xwhich - cut[i]);
        y[which] = tmpx;
      }
      wt = wt + levels[i] * (cut[i + 1] - cut[i]);
    }
  }
  
  //Rcout << "which is " << which << std::endl;
  //Rcout << "wt is " << wt << std::endl;
  //Rcout << "cut is " << cut << std::endl;
  //Rcout << "p is " << p << std::endl;
  //Rcout << "y is " << y << std::endl;
  
  if (logInd)
    y = log(y);
  return(y);
}


//[[Rcpp::export()]]
double ppwc(double q, NumericVector cuts, NumericVector levels, int lower, int logInd)
{
  double y;
  if (cuts[0]==0) {
    y = R::pexp(q, 1/levels[0], 0.0, 0.0);
  }else{
    NumericVector qq(1);
    qq[0] = q;
    y = Hpwc(qq, cuts, levels, 0.0)[0];
    if (logInd) {
      if (lower) {
        y = log(-(exp(-y)-1));
      }
      else {
        y = -y;
      }
    }
    else {
      if (lower) {
        y = -(exp(-y)-1);
      }
      else {
        y = exp(-y);
      }
    }
  }
  return(y);
}


//[[Rcpp::export()]]
NumericVector vppwc(NumericVector q, NumericVector cuts, NumericVector levels, int lower, int logInd)
{
  NumericVector y(1);
  if (cuts[0]==0) {
    y = pexp(q, levels[0], lower, logInd);
  }else{
    y = Hpwc(q, cuts, levels, 0.0);
    if (logInd) {
      if (lower) {
        y = log(-(exp(-y)-1));
      }
      else {
        y = -y;
      }
    }
    else {
      if (lower) {
        y = -(exp(-y)-1);
      }
      else {
        y = exp(-y);
      }
    }
  }
  return(y);
}

//[[Rcpp::export()]]
double dpwc(double x, NumericVector levels,  NumericVector cuts, int logInd)
{
  
  double y;
  if (cuts[0]==0) {
    y = R::dexp(x, 1/levels[0], 0.0);
  }else{
    NumericVector xx(1);
    xx[0] = x;
    y = hpwc(xx, cuts, levels, 0.0)[0] * ppwc(x, cuts, levels, 0.0, 0.0);
  }
  if (logInd)
    y = log(y);
  
  return(y);
}