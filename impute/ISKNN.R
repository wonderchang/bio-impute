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

#-----ISKNN impute-----

ISKNN <- function(xmiss, K=10, sim.method="EuDist", iter=2){
  
  xcomplete <- SKNN(xmiss, K, sim.method)
  
  miss.gene <- is.na(xmiss)
  miss.row <- which(rowSums(miss.gene)!=0)
  miss.exp <- lapply(miss.row, function(i) which(is.na(xmiss[i, ])))
  
  for(h in 1:iter) {
    xcomplete[miss.row, ] <- t(apply(xmiss[miss.row, ], 1, function(row){
      row.miss <- is.na(row)
      row.exp <- which(row.miss)
      d <- similarityCal(row[-row.exp], xcomplete[, -row.exp, drop=F], sim.method)
      id.idx <- order(d, decreasing=T)[2:(K+1)]
      id.sel <- d[id.idx]
      const <- 1/sum(id.sel)
      w <- const * id.sel
      row[row.exp] <- w %*% xcomplete[id.idx, row.exp, drop=F]
      return (row)
    }))
  }
  return(xcomplete) 
}
