#include <Rcpp.h>
#include <unordered_map>

#include "common.h"

// [[Rcpp::plugins(cpp17)]] 

// [[Rcpp::export]]
double entropyE_(Rcpp::NumericVector& p) {
  Rcpp::NumericVector plogp = p * Rcpp::log(p);
  double ent = - Rcpp::sum(plogp);
  return ent;
}

// [[Rcpp::export]]
double entropy2_(Rcpp::NumericVector& p) {
  Rcpp::NumericVector plogp = p * log2(p);
  double ent = - Rcpp::sum(plogp);
  return ent;
}


// [[Rcpp::export]]
double entropy10_(Rcpp::NumericVector& p) {
  Rcpp::NumericVector plogp = p * Rcpp::log10(p);
  double ent = - Rcpp::sum(plogp);
  return ent;
}

////////////////////////////////////////////////////////

double mutual_information_worker_(Rcpp::NumericMatrix& ps, 
                                  std::function<double(double)> log_func) {
  if (ps.ncol() != 3) {
    Rcpp::stop("Unexpected");
  }
  
  Rcpp::NumericVector p_x = ps(Rcpp::_, 0);
  Rcpp::NumericVector p_y = ps(Rcpp::_, 1);
  Rcpp::NumericVector p_xy = ps(Rcpp::_, 2);
  
  size_t n = ps.nrow();
  double s = 0.0;
  
  for (size_t i = 0; i < n; ++i) {
    double p_xy_i = p_xy[i];
    s += p_xy_i * log_func( p_xy_i / (p_x[i] * p_y[i]));
  }
  
  return s;
}

// [[Rcpp::export]]
double mutual_informationE_(Rcpp::NumericMatrix& ps) {
  return mutual_information_worker_(ps, (double(*)(double))&std::log);
}


// [[Rcpp::export]]
double mutual_information2_(Rcpp::NumericMatrix& ps) {
  return mutual_information_worker_(ps, (double(*)(double))&std::log2);
}

// [[Rcpp::export]]
double mutual_information10_(Rcpp::NumericMatrix& ps) {
  return mutual_information_worker_(ps, (double(*)(double))&std::log10);
}
