svaBatchCor <- function(dat, mmi, mm0,n.sv=NULL){
  dat <- as.matrix(dat)
  Y <- t(dat)
  if(is.null(n.sv))   n.sv <- (num.sv(dat,mmi,method="leek") - 1)
  o <- sva(dat,mmi,mm0,n.sv=n.sv)
  W <- o$sv
  alpha <- solve(t(W) %*% W) %*% t(W) %*% Y
  o$corrected <- t(Y - W %*% alpha)
  return(o)
}
