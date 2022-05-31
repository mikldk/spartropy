#include <Rcpp.h>
using namespace Rcpp;

#include "common.h"

// [[Rcpp::plugins(cpp17)]] 

// [[Rcpp::export]]
Rcpp::IntegerVector frequencies_(Rcpp::IntegerMatrix& x) {
  std::unordered_map<std::vector<int>, size_t, int_vector_hasher> counts;
  
  Rcpp::IntegerMatrix xt = Rcpp::transpose(x);
  size_t n = xt.ncol();
  
  for (size_t i = 0; i < n; ++i) {
    Rcpp::IntegerVector row = xt(Rcpp::_, i);
    std::vector<int> x = Rcpp::as< std::vector<int> >(row);
    
    counts[x]++;
  }
  
  Rcpp::IntegerVector y(counts.size());
  size_t i = 0;
  
  for (std::pair<std::vector<int>, int> element : counts) {
    y[i++] = element.second;
  }
  
  return y;
}

// [[Rcpp::export]]
Rcpp::NumericVector normalise_(Rcpp::IntegerVector& x) {
  double s = (double)Rcpp::sum(x);
  size_t n = x.size();
  
  Rcpp::NumericVector y(n);
  
  for (size_t i = 0; i < n; ++i) {
    y[i] = (double)x[i] / s;
  }
  
  return y;
}

//////////////////////////////

// https://stackoverflow.com/a/62118438
Rcpp::IntegerMatrix matrix_subset(Rcpp::IntegerMatrix& x, Rcpp::IntegerVector& y) { 
  
  // Determine the number of columns
  size_t n_cols_out = y.size();
  
  // Create an output matrix
  Rcpp::IntegerMatrix out = Rcpp::no_init(x.nrow(), n_cols_out);
  
  // Loop through each column and copy the data. 
  for (size_t z = 0; z < n_cols_out; ++z) {
    out(Rcpp::_, z) = x(Rcpp::_, y[z]);
  }
  
  return out;
}

void fill_maps(Rcpp::IntegerMatrix& x, 
               Rcpp::IntegerVector& is, 
               Rcpp::IntegerVector& js,
               
               std::unordered_map<std::vector<int>, size_t, int_vector_hasher>& n_is,
               std::unordered_map<std::vector<int>, size_t, int_vector_hasher>& n_js,
               std::unordered_map<std::pair< std::vector<int>, std::vector<int> >, size_t, int_vector_pair_hasher>& n_ijs) {
  
  size_t n = x.nrow();
  
  Rcpp::IntegerMatrix xt_is = Rcpp::transpose(matrix_subset(x, is));
  Rcpp::IntegerMatrix xt_js = Rcpp::transpose(matrix_subset(x, js));
  
  for (size_t k = 0; k < n; ++k) {
    
    Rcpp::IntegerVector row_is_k = xt_is(Rcpp::_, k);
    Rcpp::IntegerVector row_js_k = xt_js(Rcpp::_, k);
    
    std::vector<int> v_is_k = Rcpp::as< std::vector<int> >(row_is_k);
    std::vector<int> v_js_k = Rcpp::as< std::vector<int> >(row_js_k);
    
    n_is[v_is_k]++;
    n_js[v_js_k]++;
    
    std::pair< std::vector<int>, std::vector<int> > p = std::make_pair(v_is_k, v_js_k);
    n_ijs[p]++;
  }
}

// H(i | j)
// [[Rcpp::export]]
Rcpp::IntegerMatrix frequencies_2d_(
    Rcpp::IntegerMatrix& x, 
    Rcpp::IntegerVector& is, 
    Rcpp::IntegerVector& js) {
  
  std::unordered_map<std::vector<int>, size_t, int_vector_hasher> n_is;
  std::unordered_map<std::vector<int>, size_t, int_vector_hasher> n_js;
  std::unordered_map<std::pair< std::vector<int>, std::vector<int> >,
                     size_t, int_vector_pair_hasher> n_ijs;
  
  fill_maps(x, is, js,
            n_is,n_js, n_ijs);
  
  std::vector< std::tuple<int, int, int> > res_tuple;
  
  for (std::pair<std::vector<int>, int> e1 : n_is) {
    for (std::pair<std::vector<int>, int> e2 : n_js) {
      
      
      std::pair< std::vector<int>, std::vector<int> > p = std::make_pair(e1.first, e2.first);
      auto simultan = n_ijs.find(p);
      
      if (simultan == n_ijs.end()) {
        continue;
      }
      
      res_tuple.push_back(std::make_tuple(e1.second, e2.second, simultan->second));
    }
  }
  
  ///////////
  
  Rcpp::IntegerMatrix res(res_tuple.size(), 3);
  Rcpp::colnames(res) = Rcpp::CharacterVector::create("x", "y", "xy");
  
  for (size_t k = 0; k < res_tuple.size(); ++k) {
    std::tuple<int, int, int> v = res_tuple[k];
    res(k, 0) = std::get<0>(v);
    res(k, 1) = std::get<1>(v);
    res(k, 2) = std::get<2>(v);
  }
  
  return res;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix normalise_2d_(Rcpp::IntegerMatrix& x) {
  if (x.ncol() != 3) {
    Rcpp::stop("Unexpected");
  }
  
  double s = (double)Rcpp::sum(x(Rcpp::_, 2));
  size_t n = x.nrow();
  
  Rcpp::IntegerMatrix xt = Rcpp::transpose(x);
  
  Rcpp::NumericMatrix yt(3, n);
  
  for (size_t i = 0; i < n; ++i) {
    Rcpp::IntegerVector x_i = xt(Rcpp::_, i);
    
    yt(Rcpp::_, i) = NumericMatrix::create(
      (double)x_i[0] / s, 
      (double)x_i[1] / s, 
      (double)x_i[2] / s);
  }
  
  Rcpp::NumericMatrix y = Rcpp::transpose(yt);
  
  return y;
}



// H(i | j)
double mutual_information_implicit_worker(
    Rcpp::IntegerMatrix& x, 
    Rcpp::IntegerVector& is, 
    Rcpp::IntegerVector& js, 
    std::function<double(double)> log_func) {
  
  std::unordered_map<std::vector<int>, size_t, int_vector_hasher> n_is;
  std::unordered_map<std::vector<int>, size_t, int_vector_hasher> n_js;
  std::unordered_map<std::pair< std::vector<int>, std::vector<int> >,
                     size_t, int_vector_pair_hasher> n_ijs;
  
  fill_maps(x, is, js,
            n_is,n_js, n_ijs);
  
  double s = (double)x.nrow();
  double ent = 0.0;
  
  for (std::pair<std::vector<int>, int> e1 : n_is) {
    double p_x = (double)e1.second / s;
    
    for (std::pair<std::vector<int>, int> e2 : n_js) {
      
      
      std::pair< std::vector<int>, std::vector<int> > p = std::make_pair(e1.first, e2.first);
      auto simultan = n_ijs.find(p);
      
      if (simultan == n_ijs.end()) {
        continue;
      }
      
      double p_y = (double)e2.second  / s;
      double p_xy = (double)simultan->second / s;
      
      ent += p_xy * log_func( p_xy / (p_x * p_y));
    }
  }
  
  return ent;
}


// H(i | j)
// [[Rcpp::export]]
double mutual_informationE_implicit_(
    Rcpp::IntegerMatrix& x, 
    Rcpp::IntegerVector& is, 
    Rcpp::IntegerVector& js) {
  
  return mutual_information_implicit_worker(x, is, js, (double(*)(double))&std::log);
}

// H(i | j)
// [[Rcpp::export]]
double mutual_information2_implicit_(
    Rcpp::IntegerMatrix& x, 
    Rcpp::IntegerVector& is, 
    Rcpp::IntegerVector& js) {
  
  return mutual_information_implicit_worker(x, is, js, (double(*)(double))&std::log2);
}

// H(i | j)
// [[Rcpp::export]]
double mutual_information10_implicit_(
    Rcpp::IntegerMatrix& x, 
    Rcpp::IntegerVector& is, 
    Rcpp::IntegerVector& js) {
  
  return mutual_information_implicit_worker(x, is, js, (double(*)(double))&std::log10);
}
