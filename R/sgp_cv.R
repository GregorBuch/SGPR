sgp.cv <- function(X, y, group=1:ncol(X),Z=NULL, ..., nfolds=10, seed, fold, type,
                   returnY=FALSE, print.trace=FALSE) {
  
  # Fit a SGP
  tune_all        <- list(...)
  tune_all$X      <- X
  tune_all$y      <- y
  tune_all$Z      <- Z
  tune_all$group  <- group
  tune_all$return <- TRUE
  tune_all$type   <- type
  fit_all         <- do.call("sgp", tune_all)
  
  # Extract dimensions 
  X     <- fit_all$X$X
  y     <- fit_all$y
  n     <- fit_all$n
  
  # Generate folds
  if (!missing(seed)) set.seed(seed)
  if (missing(fold)) {
    if (fit_all$type=="linear") {
      fold <- sample((1:n %% nfolds)+1)
    } else {
      fold       <- integer(n)
      events     <- sum(y)
      fails      <- n-events
      fold[y==1] <- sample((1:events %% nfolds)+1)
      fold[y==0] <- sample(((events + 1:fails) %% nfolds)+1)
    }
  } else {
    nfolds <- max(fold)
  }
  
  # Performe cross-validation
  Loss <- Pred <- matrix(NA, n, length(fit_all$lambda))
  if (fit_all$type=="logit") class <- Pred
  
  tune_fold         <- list(...)
  tune_fold$group   <- fit_all$group
  tune_fold$lambdas <- fit_all$lambdas
  tune_fold$type    <- fit_all$type
  
  for (i in 1:nfolds) {
    if (print.trace) cat("Fitting SGP in fold: ", i, sep="","\n")
    
    X_out     <- X[fold==i, , drop=FALSE]
    y_out     <- y[fold==i]
    Z_out     <- Z[fold==i]
    
    tune_fold$X <- X[fold!=i, , drop=FALSE]
    tune_fold$y <- y[fold!=i]
    tune_fold$Z <- Z[fold!=i]
    
    fit_fold  <- suppressWarnings(do.call("sgp", tune_fold))
    
    intercept <- fit_fold$beta[1,]
    beta      <- fit_fold$beta[-1, ,drop=FALSE]
    eta       <- sweep(X_out %*% beta , 2, intercept , "+")
    
    if (fit_fold$type=="logit") {
      pred <- exp(eta)/(1+exp(eta))
    } else {
      pred <- eta
    }
    
    loss      <- get.loss(y_out, pred, fit_fold$type)
    pe        <- if (fit_all$type=="logit") {(pred < 0.5) == y_out} else NULL
    cv_l_l    <- length(fit_fold$lambdas)
    if(cv_l_l>ncol(Pred))cv_l_l <- ncol(Pred)
    
    Pred[fold==i, 1:cv_l_l]                             <- pred[,1:cv_l_l]
    Loss[fold==i, 1:cv_l_l]                             <- loss[,1:cv_l_l]
    if (fit_all$type=="logit") class[fold==i, 1:cv_l_l] <- pe[,1:cv_l_l]
  }
  
  # Remove failed lambdas
  conv     <- which(apply(is.finite(Pred), 2, all))
  Loss     <- Loss[, conv, drop=FALSE]
  Pred     <- Pred[, conv]
  lambdas  <- fit_all$lambdas[conv]
  
  # Fit null model
  if (is.null(fit_all$Z)) {
    null_model <- glm(y~1, family=ifelse(fit_all$type=="linear","gaussian","binomial"))
  } else {
    null_model <- glm(y~Z, family=ifelse(fit_all$type=="linear","gaussian","binomial"))
  }
  
  # Summarize
  cve      <- apply(Loss, 2, mean)
  cvse     <- apply(Loss, 2, sd) / sqrt(n)
  min      <- which.min(cve)
  null.dev <- mean(get.loss(y, predict(null_model, type="response"), fit_all$type))
  
  erg <- list(cve=cve, cvse=cvse, lambdas=lambdas, fit=fit_all, fold=fold,
              min=min, lambda.min=lambdas[min], null.dev=null.dev)
  if (fit_all$type=="logit") erg$pe <- apply(class, 2, mean)
  if (returnY) {
    if (fit_all$type=="linear") erg$Pred <- Pred + attr(y, "mean")
    else erg$Pred <- Pred
  }
  structure(erg, class="cv.sgp")
}