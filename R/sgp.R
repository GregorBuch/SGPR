sgp <- function(X, y, group=1:ncol(X), penalty=c("sgl","sgs","sgm","sge"), alpha=1/3, type=c("linear","logit"),Z=NULL,
                nlambda=100, lambda.min={if(nrow(X) > ncol(X)) 1e-4 else .05}, log.lambda=TRUE, lambdas,
                prec= 1e-4, ada_mult=2, max.iter=10000, standardize=TRUE, 
                vargamma=ifelse(pvar=="scad"|penalty=="sgs", 4, 3), grgamma=ifelse(pgr=="scad"|penalty=="sgs", 4, 3), vartau=1, grtau=1, 
                pvar=c("lasso","scad","mcp","exp"),pgr=c("lasso","scad","mcp","exp"),
                group.weight=rep(1,length(unique(group))), return=FALSE, ...) {
  
  type    <- match.arg(type) 
  if(!is.null(penalty)) penalty <- match.arg(penalty)
  if(!is.null(pvar))    pvar    <- match.arg(pvar)
  if(!is.null(pgr))     pgr     <- match.arg(pgr)
  
  # Validation and correction
  sgp        <- process.penalty(penalty,pvar,pgr,vargamma,grgamma,vartau,grtau,alpha)
  
  response   <- process.y(y, type)
  grouping   <- process.group(group,group.weight)
  predictors <- process.X(X, group)
  covariates <- process.Z(Z)
  
  n          <- length(response)
  if (nrow(predictors$X) != n) stop("Dimensions of X is not compatible with y", call.=FALSE)
  if (!is.null(Z) && (nrow(covariates$Z) != n)) stop("Dimensions of Z is not compatible with y", call.=FALSE)
  
  # Reorder processed dimensions
  to_order       <- F
  if (any(order(grouping) != 1:length(grouping))) {
    to_order          <- T  
    neworder          <- order(group)
    oldorder          <- match(1:length(group),neworder)
    
    grouping          <- grouping[neworder]
    
    predictors$X      <- predictors$X[, neworder]
    predictors$center <- predictors$center[neworder]
    predictors$scale  <- predictors$scale[neworder]
    
    group.weight      <- group.weight[order(unique(group))]
  }
  
  # Combine processed dimensions
  if (is.null(Z)) {
    dat      <- predictors
    Z_groups <- grouping
  } else {
    dat      <- list(X     = cbind(covariates$Z, predictors$X),
                     vars  = c(covariates$vars,predictors$vars),
                     center= c(covariates$center,predictors$center),
                     scale = c(covariates$scale,predictors$scale))
    Z_groups <- c(rep(0,ifelse(is.null(covariates),0,ncol(covariates))),grouping)
  }
  
  p          <- ncol(dat$X)
  
  # Setup lambdas
  if (missing(lambdas)) {
    if (is.null(covariates)) {
      lambdas <- process.lambda(dat$X,response,grouping,NULL,type,alpha,lambda.min,log.lambda,nlambda,group.weight,ada_mult)  
    } else{lambdas <- process.lambda(dat$X,response,grouping,covariates$Z,type,alpha,lambda.min,log.lambda,nlambda,group.weight,ada_mult)}
    
    #lam.max <- lambdas[1]
    own_l <- FALSE
  } else {
    #lam.max <- -1
    nlambda <- length(lambdas)
    own_l <- TRUE
  }
  
  # Indices for groupings
  J  <- as.integer(table(Z_groups))
  J0 <- as.integer(if (min(Z_groups)==0) J[1] else 0)
  JG <- as.integer(if (min(Z_groups)==0) cumsum(J) else c(0, cumsum(J)))
  if (J0) {
    lambdas[1] <- lambdas[1] + 1e-5
    own_l <- TRUE
  }
  
  # Call C++
  if (type=="linear") {
    fit <- lcdfit_linear(dat$X, response, JG, J0, lambdas, alpha, prec, ada_mult, grgamma, vargamma, grtau, vartau,
                         as.integer(max.iter), group.weight,as.integer(own_l),as.integer(sgp[[1]]),as.integer(sgp[[2]]))
    b        <- rbind(mean(y), matrix(fit[[1]], nrow=p))
    iter     <- fit[[2]]
    df       <- fit[[3]] + 1
    loss     <- fit[[4]]
    ada.prec <- fit[[5]]
  } else {
    fit <- lcdfit_logistic(dat$X, response, JG, J0, lambdas, alpha, prec, ada_mult, grgamma, vargamma, grtau, vartau,
                           as.integer(max.iter), group.weight,as.integer(own_l),as.integer(sgp[[1]]),as.integer(sgp[[2]]))
    b        <- rbind(fit[[1]], matrix(fit[[2]], nrow=p))
    iter     <- fit[[3]]
    df       <- fit[[4]]
    loss     <- fit[[5]]
    ada.prec <- fit[[6]]
  }
  
  # Remove failed lambdas
  conv    <- !is.na(iter)
  b       <- b[, conv, drop=FALSE]
  iter    <- iter[conv]
  lambdas <- lambdas[conv]
  df      <- df[conv]
  loss    <- loss[conv]
  
  if (iter[1] == max.iter) stop("The algorithm could not converge already at the first lambda. Try a higher value for alpha or another penalty", call.=FALSE)
  
  if(standardize){
    beta      <- matrix(0, nrow=1+p, ncol=ncol(b))
    beta[-1,] <- b[-1,] / dat$scale
    beta[1,]  <- b[1,] - dat$center %*% beta[-1, , drop=FALSE]
  }else{
    beta      <- b
  }
  
  # Labeling
  dimnames(beta) <- list(c("Intercept", dat$vars),
                         round(lambdas, digits=4))
  
  # Unstandardize
  if (to_order) beta[-1,] <- beta[1+oldorder,]
  
  erg <- structure(list(beta         = beta,
                        type         = type,
                        lambdas      = lambdas,
                        alpha        = alpha,
                        loss         = loss,
                        prec         = ada.prec,
                        n            = n,
                        penalty      = penalty,
                        df           = df,
                        iter         = iter,
                        group.weight = group.weight),
                        class        = "sgp")
  if (return){
    if (type == 'linear') {
      erg$y     <- response + attr(response, 'mean')
    } else {
      erg$y     <- response
    }
    erg$X     <- dat
    erg$group <- factor(grouping)
  } 
  
  return(erg)
}