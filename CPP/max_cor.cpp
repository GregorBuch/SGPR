#include <Rcpp.h>
#include <get_functions.cpp>
using namespace Rcpp;

// [[Rcpp::export]]
double max_cor(Rcpp::DoubleVector X, Rcpp::DoubleVector y, Rcpp::DoubleVector K,
               Rcpp::DoubleVector m, double alpha) {
  
  // Declarations
  int    n     = y.length();
  int    J     = K.length()-1;
  double max_z = 0;                      // highest correlation 
  
  for (int g=0; g<J; g++) {
    
    double z       = 0;                  // correlation of a variable
    Rcpp::DoubleVector Z(K[g+1]-K[g]);   // correlations of a variable group
    double max_g_z = 0;                  // highest correlation within a group
    double g_norm  = 0;                  // group norm 
    
    for (int v=K[g]; v<K[g+1]; v++) {
      Z[v-K[g]] = std::fabs(get_xty(X, y, n, v)/n);
      if (Z[v-K[g]] > max_g_z) max_g_z = Z[v-K[g]];
    }
    
    g_norm = get_norm(Z) / m[g];
    z      = max_g_z / (alpha+(1-alpha) * max_g_z/g_norm);
    if (z > max_z) max_z = z;
  }
  
  return(max_z);
}