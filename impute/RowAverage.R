RowAverage <- function(x){
  mis <- is.na(x)
  rowM <- rowMeans(x, na.rm=T)
  rowM <- matrix(rowM, nc=ncol(x), nr=nrow(x), byrow=F)
  
  x[mis] <- rowM[mis]
  x
}
