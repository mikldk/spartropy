#include <Rcpp.h>

// [[Rcpp::plugins(cpp17)]] 

#ifndef SPARTROPY_COMMON_H
#define SPARTROPY_COMMON_H

  VECTORIZED_MATH_1(log2,::log2)

  
  // https://stackoverflow.com/a/72073933
  class int_vector_hasher {
  public:
    std::size_t operator()(std::vector<int> const& vec) const {
      std::size_t seed = vec.size();
      for (auto x : vec) {
        x = ((x >> 16) ^ x) * 0x45d9f3b;
        x = ((x >> 16) ^ x) * 0x45d9f3b;
        x = (x >> 16) ^ x;
        seed ^= x + 0x9e3779b9 + (seed << 6) + (seed >> 2);
      }
      return seed;
    }
  };

  class int_vector_pair_hasher {
  public:
    std::size_t operator()(std::pair< std::vector<int> , std::vector<int> > const& p) const {
      
      std::vector<int> v1 = p.first;
      std::vector<int> v2 = p.second;
      
      std::size_t seed = v1.size() + v2.size();
      
      for (auto x : v1) {
        x = ((x >> 16) ^ x) * 0x45d9f3b;
        x = ((x >> 16) ^ x) * 0x45d9f3b;
        x = (x >> 16) ^ x;
        seed ^= x + 0x9e3779b9 + (seed << 6) + (seed >> 2);
      }
      
      for (auto x : v2) {
        x = ((x >> 16) ^ x) * 0x45d9f3b;
        x = ((x >> 16) ^ x) * 0x45d9f3b;
        x = (x >> 16) ^ x;
        seed ^= x + 0x9e3779b9 + (seed << 6) + (seed >> 2);
      }
      
      return seed;
    }
  };
#endif
