get.loss <- function(y, pred, type) {
  n <- length(y)
  if (type=="linear") {
    loss <- (y-pred)^2
  } else {
      pred[pred < 0.00001] <- 0.00001
      pred[pred > 0.99999] <- 0.99999
      
      if (is.matrix(pred)) {
        loss <- matrix(NA, nrow=nrow(pred), ncol=ncol(pred))
        loss[y==1,] <- -2*log(pred[y==1, , drop=FALSE])
        loss[y==0,] <- -2*log(1-pred[y==0, , drop=FALSE])
      } else {
        loss <- double(length(y))
        loss[y==1] <- -2*log(pred[y==1])
        loss[y==0] <- -2*log(1-pred[y==0])
      }
  } 
  loss
}
