#-----caculate distance 
similarityCal<-function(vec, mat, method="EuDist"){
  methods<-c("EuDist","cor","Angle")
  switch(match.arg(method,methods),
         EuDist=1/sqrt(rowSums((mat-matrix(vec,nc=length(vec),nr=nrow(mat),byrow=TRUE))^2)),
         cor=apply(mat,1,function(i) abs(cor(i,vec,use="everything",method="pearson"))),
         Angle=apply(mat,1,function(i) abs(sum(i * vec)/sqrt(sum(i^2)*sum(vec^2))))
  )
}

#-----WSLLS impute-----
WSLLS <- function(xmiss, K=10, sim.method="EuDist", order=2){

  miss.gene <- is.na(xmiss)
  miss.row <- which(rowSums(miss.gene)!=0)
  x.complete <- xmiss[-miss.row, ]
  x.incomplete <- xmiss[miss.row, ]
  miss.rate <- rowSums(is.na(x.incomplete))
  miss.rate.order <- order(miss.rate)
  therhold <- sum(miss.rate)/length(miss.rate)
  x.incomplete <- x.incomplete[miss.rate.order, ]
  x.imputed <- t(apply(x.incomplete, 1, function(row){
    row.miss <- which(is.na(row))
    sim <- similarityCal(row[-row.miss], x.complete[, -row.miss], sim.method)
    sim.id <- order(sim, decreasing=T)[1:K]
    x.diag <- diag(sim[sim.id]**order)
    row[row.miss] <- t(x.diag %*% x.complete[sim.id, row.miss, drop=FALSE]) %*% ginv(t(x.diag %*% x.complete[sim.id, -row.miss, drop=FALSE])) %*%row[-row.miss, drop=FALSE]
    if(length(row.miss) <= therhold){
      x.complete <<- rbind(x.complete, row)
    }
    return(row)
  }))
  
  xmiss[miss.row[miss.rate.order], ] <- x.imputed
  
  return(xmiss)
}

#-----IWSLLS impute-----
IWSLLS <- function(xmiss, K=10, sim.method="EuDist", iter=2, order=2){
  
  xcomplete <- WSLLS(xmiss, K, sim.method, order)

  miss.gene <- is.na(xmiss)
  miss.row <- which(rowSums(miss.gene)!=0)

  for(h in 1:iter) {
    xcomplete[miss.row, ] <- t(apply(xmiss[miss.row, ], 1, function(row){
      row.miss <- which(is.na(row))
      sim <- similarityCal(row[-row.miss], xcomplete[, -row.miss], sim.method)
      sim.id <- order(sim, decreasing=T)[2:(K+1)]
      x.diag <- diag(sim[sim.id]**order)
      row[row.miss] <- t(x.diag %*% xcomplete[sim.id, row.miss, drop=FALSE]) %*% ginv(t(x.diag %*% xcomplete[sim.id, -row.miss, drop=FALSE])) %*% row[-row.miss, drop=FALSE]
      return(row)
    }))
  }
  return(xcomplete) 
}
