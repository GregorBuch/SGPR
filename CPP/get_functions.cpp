#include <Rcpp.h>
// using namespace Rcpp;

// [[Rcpp::export]]
double get_loss(Rcpp::DoubleVector r, int n) {
  double l  = 0;
  for (int i=0;i<n;i++) l += std::pow(r[i],2);
  return(l);
}

// [[Rcpp::export]]
double get_xty(Rcpp::DoubleVector x, Rcpp::DoubleVector y, int n, int v) {
  double z   = 0;
  int    nv  = n*v;
  for (int i=0; i<n; i++) z += x[nv+i] * y[i];
  return(z);
}

// [[Rcpp::export]]
double get_norm(Rcpp::DoubleVector X) {
  int    p      = X.length();
  double X_norm = 0;
  for (int v=0; v<p; v++) X_norm += std::pow(X[v],2);
         X_norm = std::sqrt(X_norm);
  return(X_norm);
}
