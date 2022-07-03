process.Z <- function(Z) {
  
  # Validation and correction of Z
  if (is.null(Z)) return(NULL) 
  if(anyNA(Z)) stop("Missing data detected in Z. Please remove or impute cases with NA's.", call.=FALSE)
  if(!is.matrix(Z)) {
    tmp <- try(Z <- model.matrix(~0+., Z), silent=TRUE)
    if(inherits(tmp, "try-error")) stop("Z must be a matrix or able to be coerced to a matrix", call.=FALSE)
  }
  if(mode(Z)=="integer") mode(Z) <- "double"
  
  vars <- if(is.null(colnames(Z))) paste0("Covariate ", 1:ncol(Z)) else colnames(Z)
  
  Zsd     <- scale(Z)
  center <- attributes(Zsd)$'scaled:center'
  scale  <- attributes(Zsd)$'scaled:scale'
  if (length(which(scale > 1e-6)) != ncol(X)) {
    stop ("Please remove constants (scale < 1e-6) from Z.", call.=FALSE)
  }
  
  return(list(Z=Zsd, vars=vars, center=center, scale=scale))
}
