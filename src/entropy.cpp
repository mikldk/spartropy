#include <Rcpp.h>
#include <unordered_map>

#include "common.h"

// [[Rcpp::plugins(cpp17)]] 

// [[Rcpp::export]]
double entropy(Rcpp::NumericVector p) {
  Rcpp::NumericVector plogp = p * Rcpp::log(p);
  double ent = - Rcpp::sum(plogp);
  return ent;
}

// [[Rcpp::export]]
double entropy2(Rcpp::NumericVector p) {
  Rcpp::NumericVector plogp = p * log2(p);
  double ent = - Rcpp::sum(plogp);
  return ent;
}


// [[Rcpp::export]]
double entropy10(Rcpp::NumericVector p) {
  Rcpp::NumericVector plogp = p * Rcpp::log10(p);
  double ent = - Rcpp::sum(plogp);
  return ent;
}

////////////////////////////////////////////////////////


// [[Rcpp::export]]
double mutual_information2(Rcpp::NumericMatrix ps) {
  if (ps.ncol() != 3) {
    Rcpp::stop("Unexpected");
  }
  
  Rcpp::NumericVector p_x = ps(Rcpp::_, 0);
  Rcpp::NumericVector p_y = ps(Rcpp::_, 1);
  Rcpp::NumericVector p_xy = ps(Rcpp::_, 2);
  
  
  return Rcpp::sum(p_xy * log2( p_xy / (p_x * p_y)));
}
