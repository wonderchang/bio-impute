SKNN <- function(xmiss, K=10, sim.method="EuDist"){
  miss.gene <- is.na(xmiss)
  miss.row <- which(rowSums(miss.gene) != 0)
  xincomplete <- xmiss[miss.row, ]
  miss.inc <- is.na(xincomplete)
  miss.origin <- order(rowSums(miss.inc))
  xincomplete <- xincomplete[order(rowSums(miss.inc)), ]
  xcomplete <- xmiss[-miss.row, ]
  xtmp <- matrix(nc=ncol(xincomplete))
  xmiss[miss.row[miss.origin], ] <- t(apply(xincomplete, 1, function(row){
    row.miss <- is.na(row)
    row.exp <- which(row.miss)
    d <- similarityCal(row[-row.exp], xcomplete[, -row.exp, drop=F], sim.method)
    id.idx <- order(d, decreasing=T)[1:K]
    id.sel <- d[id.idx]
    const <- sum(id.sel)
    w <- 1/const*id.sel
    row[row.exp] <- w%*%xcomplete[id.idx, row.exp, drop=F]
    xcomplete <<- rbind(xcomplete, row)
    return (row)
  }))
  
  return (xmiss)
}
similarityCal<-function(vec, mat, method="EuDist"){
  methods<-c("EuDist","cor","Angle")
  switch(match.arg(method,methods),
         EuDist=1/sqrt(rowSums((mat-matrix(vec,nc=length(vec),nr=nrow(mat),byrow=TRUE))^2)),
         cor=apply(mat,1,function(i) abs(cor(i,vec,use="everything",method="pearson"))),
         Angle=apply(mat,1,function(i) abs(sum(i * vec)/sqrt(sum(i^2)*sum(vec^2))))
  )
}
