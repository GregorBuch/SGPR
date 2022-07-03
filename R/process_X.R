process.X <- function(X, group) {
  
  # Validation and correction of X
  if(anyNA(X)) stop("Missing data detected in X. Please remove or impute cases with NA's.", call.=FALSE)
  if(!is.matrix(X)) {
    tmp <- try(X <- model.matrix(~0+., X), silent=TRUE)
    if(inherits(tmp, "try-error")) stop("X must be a matrix or able to be coerced to a matrix", call.=FALSE)
  }
  if ( dim(X)[2] <= 1) stop("X must be a matrix with more than one column")
  if(mode(X)=="integer") mode(X) <- "double"
  
  if(length(group) != ncol(X)) stop ("Dimensions of group is not compatible with X", call.=FALSE)
  
  vars <- if(is.null(colnames(X))) paste0("Variable ", 1:ncol(X)) else colnames(X)
  
  Xsd    <- scale(X)
  center <- attributes(Xsd)$'scaled:center'
  scale  <- attributes(Xsd)$'scaled:scale'
  if (length(which(scale > 1e-6)) != ncol(X)) {
    stop ("Please remove constants (scale < 1e-6) from X.", call.=FALSE)
  }
  
  return(list(X=Xsd, vars=vars, center=center, scale=scale))
}
