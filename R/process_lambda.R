process.lambda <- function(X, y,group, Z, type, alpha, lambda.min, log.lambda, nlambda, group.weight, ada_mult) {
  
  # Validation and correction of lambda related input
  if (nlambda < 2){
    warning("nlambda must be at least 2 and was set to its default value.")
    nlambda  <- 100
  } 
  if (ada_mult < 1){
    warning("ada_mult must be >= 1 and was set to its default value.")
    ada_mult <- 2
  } 
  
  # GLM with unpenalized variables
  n <- length(y)
  if (is.null(Z)) {
    fit <- glm(y~1, family=ifelse(type=="linear","gaussian","binomial"))
  } else {
    fit <- glm(y~Z, family=ifelse(type=="linear","gaussian","binomial"))
  }
  
  # Determine lambda.max
  if (type=="linear") {
    r <- fit$residuals
  } else {
    w <- fit$weights
    w <- w + 3*w*(1-alpha)*alpha
    if (max(w) < 1e-4) stop("Intercept or variables in Z result in a saturated model: no residuals left for selection", call.=FALSE)
    r <- residuals(fit, "working")*w
  }
  
  lambda.max <- max_cor(X, r, c(0, cumsum(table(group))), as.double(group.weight),alpha)
  
  # Generate lambda sequence
  if (log.lambda) { 
    if (lambda.min==0) {
      lambdas <- c(exp(seq(log(lambda.max), log(.001*lambda.max), length=nlambda-1)), 0)
    } else {
      lambdas <- exp(seq(log(lambda.max), log(lambda.min*lambda.max), length=nlambda))
    }
  } else { 
    if (lambda.min==0) {
      lambdas <- c(seq(lambda.max, 0.001*lambda.max, length = nlambda-1), 0)
    } else {
      lambdas <- seq(lambda.max, lambda.min*lambda.max, length = nlambda)
    }
  }
  lambdas
}