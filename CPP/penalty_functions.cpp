#include <Rcpp.h>
// using namespace Rcpp;

// [[Rcpp::export]]
double mcp(double theta, double l, double gamma) {
  theta = std::fabs(theta);
  if (theta < gamma*l) return(l-theta/gamma);
  else return(0);
}

// [[Rcpp::export]]
double scad(double theta, double l, double gamma) {
  theta = std::fabs(theta);
  if (theta <= l) return(l);
  if (l<theta && theta<l*gamma) return((gamma*l-theta)/(gamma-1));
  else return(0);
}

// [[Rcpp::export]]
double ep(double theta, double l, double tau) {
  if (l != 0){
    theta = std::fabs(theta);
    return(l*exp(-tau/l*theta));
  }
  return(0);
}

// [[Rcpp::export]]
double lasso(double theta, double l) {
  if (theta >  l) return(theta-l);
  if (theta < -l) return(theta+l);
  return(0);
}

// [[Rcpp::export]]
double approx(double l_v, double l_g, double g_tau, double v_tau, double g_gamma, double v_gamma, 
              double old_b, double z_v, double g_norm, double g_norm_active,
              int method_v, int method_g) {
  
  double pen_v = 0;   // penalty on the variable-level
  double pen_g = 0;   // penalty on the group-level
  double pen_z = 0;   // penalty for correlation z
  
  switch(method_v) {
  case 1 :
    pen_v = l_v;
    break;
  case 2 :
    pen_v = mcp(old_b,l_v,v_gamma);
    break;
  case 3 :
    pen_v = ep(old_b,l_v,v_tau);
    break;
  default:
    pen_v = scad(old_b,l_v,v_gamma);
  }
  
  switch(method_g) {
  case 1 :
    pen_g = l_g;
    break;
  case 2 :
    pen_g = mcp(g_norm_active,l_g,g_gamma);
    break;
  case 3 :
    pen_g = ep(g_norm_active,l_g,g_tau);
    break;
  default:
    pen_g = scad(g_norm_active,l_g,g_gamma);
  }
  
  if (g_norm !=0) {pen_z = pen_v + pen_g*std::fabs(z_v)/g_norm;
  } else {
    pen_z = pen_v + pen_g;
  }
  
  return(pen_z);
}