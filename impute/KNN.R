KNN <- function(x, K=15, sim.method="EuDist"){
  miss.gene <- is.na(x)
  miss.row <- which(rowSums(miss.gene) != 0)
  xcomplete <- x[-miss.row, ]
  xincomplete <- x[miss.row, ]
  
  x[miss.row, ] <- t(apply(xincomplete, 1, function(row){
    row.miss <- is.na(row)
    row.exp <- which(row.miss)
    d <- similarityCal(row[-row.exp], xcomplete[, -row.exp, drop=F], sim.method)
    id.idx <- order(d, decreasing=T)[1:K]
    id.sel <- d[id.idx]
    const <- sum(id.sel)
    w <- 1/const*id.sel
    w <- matrix(w, nc=length(w), nr=1)
    row[row.exp] <- w%*%xcomplete[id.idx, row.exp, drop=F]
    return (row)
  }))
  
  return (x)
}

similarityCal<-function(vec, mat, method="EuDist"){
  switch(match.arg(method,methods),
         EuDist=1/sqrt(rowSums((mat-matrix(vec,nc=length(vec),nr=nrow(mat),byrow=TRUE))^2)),
         PCor=apply(mat,1,function(i) abs(cor(i,vec,use="everything",method="pearson"))),
         Angle=apply(mat,1,function(i) abs(sum(i * vec)/sqrt(sum(i^2)*sum(vec^2))))
}

