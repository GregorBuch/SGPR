#include <Rcpp.h>
// using namespace Rcpp;

// [[Rcpp::export]]
double inv_log_link(double logit) {
  if (logit > 10) {
    return(1);
  } else if (logit < -10) {
    return(0);
  } else {
    return(exp(logit) / (1+exp(logit)));
  }
}

// [[Rcpp::export]]
Rcpp::List lcdfit_logistic(Rcpp::DoubleVector X, Rcpp::DoubleVector y,
                           Rcpp::DoubleVector G1, int G0, Rcpp::DoubleVector lam, double alpha, 
                           double prec,double ada_mult, double g_gamma, double v_gamma, double g_tau, double v_tau,
                           int max_iter, Rcpp::DoubleVector m, int own_l,
                           int method_v, int method_g) {
  
  // Initialize objects
  int l_init;                          // Start index of lambda sequence (0 or 1)
  int l_iter;                          // Iteration count for lambda
  int active_plus;                     // Indication whether the active set has changed
  int n           = y.length();        // No. observations
  int L           = lam.length();      // No. lambdas
  int J           = G1.length() - 1;   // No. variable groups to penalize
  int p           = X.length()/n;      // No. variables
  int max_l_iter  = max_iter/L;        // Maximum allowed iteration for a lambda
  int tot_iter    = 0;                 // Total no. iterations
  
  double l_v;                          // Value of active lambda * alpha
  double l_g;                          // Value of active lambda * (1-alpha)
  double g_norm_tmp;                   // Temporary variable to build l2-norm of a group
  double g_norm_active;                // L2-norm of active variables of a group
  double z_v;                          // Correlation of variable v with r
  double pen_z;                        // Penalisation for z
  double delta;                        // Difference between actual and previous iteration
  double delta_i;                      // Difference between actual and previous iteration of observation i
  double max_delta;                    // Largest difference between actual and previous iteration
  double old_b0;                       // Intercept of the previous iteration
  double loss0;                        // Null deviance of empty model
  double w    = 0.25;                  // IRLS-weight set to 0.25 = Hessian upper-bound (Krishnapuram & Hartemink, 2005)
  double w_v  = 1-alpha + w*alpha;     // IRLS-weight for the variable level
  double w_g  = alpha+w-w*alpha;       // IRLS-weight for the group level
  double av_y = 0;                     // mean(y)
  
  Rcpp::NumericVector dof (L);         // "degrees of freedom" for each lambda
  Rcpp::NumericVector loss (L);        // Loss / error for each lambda
  
  Rcpp::IntegerVector iter (L);        // No. iterations for each lmabda
  Rcpp::IntegerVector active (p);      // Indicator variable for variables for active set strategy
  Rcpp::IntegerVector active_g (J);    // Indicator variable for groups for active set strategy
  
  Rcpp::DoubleVector  b (L*p);         // Beta coefficient for each variable and lambda
  Rcpp::DoubleVector  b0 (L);          // Intercept for each lambda
  Rcpp::DoubleVector  old_b (p);       // Beta coefficient for each variable of the previous iteration
  Rcpp::DoubleVector  z (p);           // Correlation of each variable with r
  Rcpp::DoubleVector  r (n);           // Residuals 
  Rcpp::DoubleVector  g_norm (J);      // L2-norm of each variable group
  Rcpp::DoubleVector  ada_prec (L);    // Attained precision of each lambda
  Rcpp::DoubleVector  logit (n*L);     // Probability of y = 1 for each observatio and lambda
  
  // indicator variables:
  // i   for individuals [0:(n-1)]     
  // v   for variables   [0:(p-1)] OR [0:(G1-1)] OR [0:(G0-1)]    
  // g   for groups      [0:(J-1)]     
  // l   for lambdas     [l_init:(L-1)]     
  
  for (int v=0; v<p; v++) active[v]   = 0;
  for (int g=0; g<J; g++) active_g[g] = 0;
  for (int i=0; i<n; i++){
    r[i]  = y[i];
    av_y += y[i];
  }
  av_y = av_y/n;
  for (int i=0; i<n; i++) loss0 -= 2*y[i]*log(av_y) + 2*(1-y[i])*log(1-av_y);
  
  old_b0 = b0[0] = log(av_y/(1-av_y));
  
  // In case of a self-defined lambda sequence, start with the first value
  if (own_l) {
    l_init  = 0;
  } else {
    loss[0] = loss0;
    l_init  = 1;
  }
  
  // Adaptive rescaling of tuning parameters (Breheny & Huang, 2011)
  g_gamma = g_gamma / w;
  v_gamma = v_gamma / w;
  g_tau   = g_tau   * w;
  v_tau   = v_tau   * w;
  
  // Loop over lambda sequence
  for (int l=l_init; l<L; l++) {
    
    Rcpp::checkUserInterrupt();
    
    l_iter      = 0;
    ada_prec[l] = prec;
    
    // Warm start after first lambda
    if (l != 0) {
      old_b0 = b0[l-1];
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
          l_iter      = 0;
          ada_prec[l] = ada_prec[l] * ada_mult;
        }
        
        // Approximate loss
        loss[l] = 0;
        for (int i=0; i<n; i++) {
          double mu = inv_log_link(logit[i]);
          r[i] = (y[i] - mu) / w;
          if (y[i]==1) loss[l] -= 2*log(mu);
          else loss[l] -= 2*log(1-mu);
        }
        
        // Check for saturation
        if (loss[l]/loss0 < 0.05) {
          
          //warning("Model saturated; exiting...");
          Rcpp::warning("Warning: Lambda sequence reduced as saturated models were created.");
          for (int ll=l; ll<L; ll++) iter[ll] = NA_INTEGER;
          tot_iter = max_iter;
          break;
        }
        
        // Update intercept
        for (int i=0; i<n; i++) delta += r[i];
        delta = delta/n;
        b0[l] = delta + old_b0;
        for (int i=0; i<n; i++) {
          r[i]     -= delta;
          logit[i] += delta;
        }
        dof[l]    = 1;
        max_delta = std::fabs(delta);
        
        // Update coefficients of unpenalized variable group
        for (int v=0; v<G0; v++) {
          delta    = get_xty(X, r, n, v)/n;
          b[l*p+v] = delta + old_b[v];
          
          if (std::fabs(delta) > max_delta) max_delta = std::fabs(delta);
          for (int i=0; i<n; i++) {
            delta_i   = delta * X[n*v+i];
            r[i]     -= delta_i;
            logit[i] += delta_i;
          }
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
                               old_b[v], z[v], g_norm[g]*w_g, g_norm_active,
                               method_v, method_g);
                
                b[l*p+v] = lasso(w_v*z_v, pen_z)/w_v;
                
                delta    = b[l*p+v] - old_b[v];
                if (delta != 0) {
                  if (std::fabs(delta) > max_delta) max_delta = std::fabs(delta);
                  for (int i=0; i<n; i++){
                    delta_i   = delta * X[n*v+i];
                    r[i]     -= delta_i;
                    logit[i] += delta_i;
                  }
                }
                
                dof[l] += std::fabs(b[l*p+v]) / std::fabs(z_v);
              }
            } 
          }        
        }
        
        // Check for convergence      
        old_b0 = b0[l];
        for (int v=0; v<p; v++) old_b[v] = b[l*p+v];
        if (max_delta < ada_prec[l]) break;
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
                           old_b[v], z[v], g_norm[g]*w_g, g_norm_active,
                           method_v, method_g);
            
            if (w_v*std::fabs(z_v) > pen_z) {
              active_plus++;
              active[v] = 1;
              b[l*p+v] = lasso(w_v*z_v, pen_z)/w_v;
              for (int i=0; i<n; i++){
                delta_i   = b[l*p+v] * X[n*v+i];
                r[i]     -= delta_i;
                logit[i] += delta_i;
              }
              active_g[g] = 1;
              
            }
          }
        }
      }
      
      if (active_plus==0) {
        // update logit if you want to report them
        break;
      }
      old_b0 = b0[l];
      for (int v=0; v<p; v++) old_b[v] = b[l*p+v];
    }
  }
  
  return Rcpp::List::create(Rcpp::Named("intercept") = b0,
                            Rcpp::Named("beta")      = b,
                            Rcpp::Named("iter")      = iter,
                            Rcpp::Named("dof")       = dof,
                            Rcpp::Named("loss")      = loss, 
                            Rcpp::Named("adaprec")   = ada_prec);
}
