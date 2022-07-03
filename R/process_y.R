process.y <- function(y, type) {
  
  if (anyNA(y)) stop("Missing data detected in y. Please remove or impute cases with NA's.", call.=FALSE)
  
  if (is.data.frame(y)) y <- as.matrix(y)
  if (is.matrix(y)) {
    if (dim(y)[2]!=1) stop("y is multidinemsional, but only one-dimensional outcomes are supported!", call.=FALSE)
    y <- c(y)
  }
  
  # Logistic regression
  if (typeof(y) != "logical" && type=="logit") {
    tmp <- table(y)
    if (length(tmp) > 2) stop("If type='logit', y must be a binary variable, but this is not the case!", call.=FALSE)
    if (!identical(names(tmp), c("0", "1"))) {
      message(paste0("Logistic regression modeling Pr(y=", names(tmp)[2], ")"))
      y <- as.double(as.character(y) == names(tmp)[2])
    }
  }
  
  # Convert to double & linear regression
  if (typeof(y) != "double") {
    tryCatch(storage.mode(y) <- "double", warning=function(w) {
      stop("y must be numeric or be able to be converted to a numeric!", call.=FALSE)})
  }
  if (type=="linear") {
    ybar <- mean(y)
    y <- y - ybar
    attr(y, "mean") <- ybar
  }
  y
}
