#ifndef COMMONF_H
#define COMMONF_H

#include <RcppArmadillo.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include <Rmath.h>
#include <iostream>

using namespace Rcpp;

double hpwc_double(double x, NumericVector cuts, NumericVector levels, int logInd);
double Hpwc_double(double x, NumericVector cuts, NumericVector levels, int logInd);
double ppwc_double(double q, NumericVector cuts, NumericVector levels, int lower, int logInd);
double ppwc(double q, NumericVector cuts, NumericVector levels, int lower, int logInd);
#endif




