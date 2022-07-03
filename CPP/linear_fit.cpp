#include <Rcpp.h>
// using namespace Rcpp;

// [[Rcpp::export]]
Rcpp::List lcdfit_linear(Rcpp::DoubleVector X, Rcpp::DoubleVector y,
                         Rcpp::DoubleVector G1, int G0, Rcpp::DoubleVector lam, double alpha, 
                         double prec,double ada_mult, double g_gamma, double v_gamma, double g_tau, double v_tau,
                         int max_iter, Rcpp::DoubleVector m, int own_l,
                         int method_v, int method_g) {
  
  // Initialize objects
  int l_init;                         // Start index of lambda sequence (0 or 1)
  int l_iter;                         // Iteration count for lambda
  int active_plus;                    // Indication whether the active set has changed
  int n           = y.length();       // No. observations
  int L           = lam.length();     // No. lambdas
  int J           = G1.length() - 1;  // No. variable groups to penalize
  int p           = X.length()/n;     // No. variables
  int max_l_iter  = max_iter/L;       // Maximum allowed iteration for a lambda
  int tot_iter    = 0;                // Total no. iterations
  
  double l_v;                         // Value of active lambda * alpha
  double l_g;                         // Value of active lambda * (1-alpha)
  double g_norm_tmp;                  // Temporary variable to build l2-norm of a group
  double g_norm_active;               // L2-norm of active variables of a group
  double z_v;                         // Correlation of variable v with r
  double pen_z;                       // Penalisation for z
  double delta;                       // Difference between actual and previous iteration
  double max_delta;                   // Largest difference between actual and previous iteration
  double rss  = get_loss(y,n);        // Residual sum of squares
  double sd_y = std::sqrt(rss/n);     // Standard deviation of y
  
  Rcpp::NumericVector dof (L);        // "degrees of freedom" for each lambda
  Rcpp::NumericVector loss (L);       // Loss / error for each lambda
  
  Rcpp::IntegerVector iter (L);       // No. iterations for each lmabda
  Rcpp::IntegerVector active (p);     // Indicator variable for variables for active set strategy
  Rcpp::IntegerVector active_g (J);   // Indicator variable for groups for active set strategy
  
  Rcpp::DoubleVector  b (L*p);        // Beta coefficient for each variable and lambda
  Rcpp::DoubleVector  old_b (p);      // Beta coefficient for each variable of the previous iteration
  Rcpp::DoubleVector  z (p);          // Correlation of each variable with r
  Rcpp::DoubleVector  r (n);          // Residuals 
  Rcpp::DoubleVector  g_norm (J);     // L2-norm of each variable group
  Rcpp::DoubleVector  ada_prec (L);   // Attained precision of each lambda
  
  // indicator variables:
  // i   for individuals [0:(n-1)]     
  // v   for variables   [0:(p-1)] OR [0:(G1-1)] OR [0:(G0-1)]    
  // g   for groups      [0:(J-1)]     
  // l   for lambdas     [l_init:(L-1)]     
  
  for (int v=0; v<p; v++) active[v]   = 0;
  for (int g=0; g<J; g++) active_g[g] = 0;
  for (int i=0; i<n; i++) r[i]        = y[i];
  
  // In case of a self-defined lambda sequence, start with the first value
  if (own_l) {
    l_init  = 0;
  } else {
    loss[0] = rss;
    l_init  = 1;
  }
  
  // Loop over lambda sequence
  for (int l=l_init; l<L; l++) {
    
    Rcpp::checkUserInterrupt();
    
    l_iter      = 0;
    ada_prec[l] = prec;
    
    // Warm start after first lambda
    if (l != 0) {
      for (int v=0; v<p; v++) old_b[v] = b[(l-1)*p+v];
    }
    
    // Initial approximation of l2-norm of each group 
    for (int g=0; g<J; g++) {
      g_norm_tmp = 0;
      for (int v=G1[g]; v<G1[g+1]; v++) {
        z[v]        = get_xty(X, r, n, v)/n + old_b[v];
        g_norm_tmp += std::pow(z[v],2);
      }
      g_norm[g]     = std::sqrt(g_norm_tmp);
    }
    
    // Activate variables until convergence or max_iter is reached
    while (tot_iter < max_iter) {
      
      // Update betas until convergence or max_iter is reached
      while (tot_iter < max_iter) {
        
        l_iter++;
        iter[l]++;
        tot_iter++;
        dof[l]    = 0;
        max_delta = 0;
        
        // Adaptive precision
        if (l_iter > max_l_iter && ada_mult !=1){
          l_iter     = 0;
          ada_prec[l] = ada_prec[l] * ada_mult;
        }
        
        // Update coefficients of unpenalized variable group
        for (int v=0; v<G0; v++) {
          delta    = get_xty(X, r, n, v)/n;
          b[l*p+v] = delta + old_b[v];
          
          if (std::fabs(delta) > max_delta) max_delta = std::fabs(delta);
          for (int i=0; i<n; i++) r[i] -= delta * X[n*v+i];
          dof[l]++;
        }
        
        // Update coefficients of penalized variable groups
        for (int g=0; g<J; g++) {
          
          if(active_g[g]){
            l_v  = lam[l] * alpha;
            l_g  = lam[l] * m[g] * (1-alpha);
            
            // Approximate l2-norm of the group 
            g_norm_active = 0;
            g_norm_tmp    = 0;
            for (int v=G1[g]; v<G1[g+1]; v++) {
              z[v]           = get_xty(X, r, n, v)/n + old_b[v];
              g_norm_tmp    += std::pow(z[v],2);
              g_norm_active += std::pow(old_b[v],2);
              
            }
            g_norm[g]     = std::sqrt(g_norm_tmp);
            g_norm_active = std::sqrt(g_norm_active);
            
            // Update coefficients of the variables within the group
            for (int v=G1[g]; v<G1[g+1]; v++) {
              if (active[v]) {
                
                z_v   = get_xty(X, r, n, v)/n + old_b[v];
                
                pen_z = approx(l_v, l_g, g_tau, v_tau, g_gamma, v_gamma,
                               old_b[v], z[v], g_norm[g], g_norm_active,
                               method_v, method_g);
                
                b[l*p+v] = lasso(z_v, pen_z);
                
                delta    = b[l*p+v] - old_b[v];
                if (delta != 0) {
                  if (std::fabs(delta) > max_delta) max_delta = std::fabs(delta);
                  for (int i=0; i<n; i++) r[i] -= delta*X[n*v+i];
                }
                
                dof[l] += std::fabs(b[l*p+v]) / std::fabs(z_v);
              }
            } 
          }        
        }
        
        // Check for convergence      
        for (int v=0; v<p; v++) old_b[v] = b[l*p+v];
        if (max_delta < ada_prec[l]*sd_y) break;
      }
      
      // Active set strategy
      active_plus = 0;
      for (int g=0; g<J; g++) {
        
        l_v = lam[l] * alpha;
        l_g = lam[l] * m[g] * (1-alpha);
        
        // Approximate l2-norm of the group 
        g_norm_active = 0;
        g_norm_tmp    = std::pow(g_norm[g],2);
        for (int v=G1[g]; v<G1[g+1]; v++) {
          if (active[v]==0) {
            z[v]         = get_xty(X, r, n, v)/n + old_b[v];
            g_norm_tmp  += std::pow(z[v],2);
          }
          g_norm_active += std::pow(old_b[v],2);
        }
        g_norm[g]     = std::sqrt(g_norm_tmp);
        g_norm_active = std::sqrt(g_norm_active);
        
        // Update coefficients of the variables within the group
        for (int v=G1[g]; v<G1[g+1]; v++) {
          if (active[v]==0) {
            
            z_v   = get_xty(X, r, n, v)/n + old_b[v];
            pen_z = approx(l_v, l_g, g_tau, v_tau, g_gamma, v_gamma,
                           old_b[v], z[v], g_norm[g], g_norm_active,
                           method_v, method_g);
            
            if (std::fabs(z_v) > pen_z) {
              active_plus++;
              active[v] = 1;
              b[l*p+v] = lasso(z_v, pen_z);
              for (int i=0; i<n; i++) r[i] -= b[l*p+v] * X[n*v+i];
              active_g[g] = 1;
              
            }
          }
        }
      }
      
      if (active_plus==0) {
        loss[l] = get_loss(r, n);
        break;
      }
      for (int v=0; v<p; v++) old_b[v] = b[l*p+v];
    }
  }
  
  return Rcpp::List::create(Rcpp::Named("beta")    = b,
                            Rcpp::Named("iter")    = iter,
                            Rcpp::Named("dof")     = dof,
                            Rcpp::Named("loss")    = loss, 
                            Rcpp::Named("adaprec") = ada_prec);
}

