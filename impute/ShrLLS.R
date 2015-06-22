#-----caculate distance 
similarityCal<-function(vec, mat, method="EuDist"){
  methods<-c("EuDist","cor","Angle")
  switch(match.arg(method,methods),
         EuDist=1/sqrt(rowSums((mat-matrix(vec,nc=length(vec),nr=nrow(mat),byrow=TRUE))^2)),
         cor=apply(mat,1,function(i) abs(cor(i,vec,use="everything",method="pearson"))),
         Angle=apply(mat,1,function(i) abs(sum(i * vec)/sqrt(sum(i^2)*sum(vec^2))))
  )
}

#-----ShrLLS impute-----
ShrLLS <- function(xmiss, K=10, sim.method="EuDist"){

  miss.gene <- is.na(xmiss)
  miss.row <- which(rowSums(miss.gene)!=0)
  x.complete <- xmiss[-miss.row, ]
  x.incomplete <- xmiss[miss.row, ]
  
  x.imputed <- t(apply(x.incomplete, 1, function(row){
    row.miss <- which(is.na(row))
    sim <- similarityCal(row[-row.miss], x.complete[, -row.miss], sim.method)
    sim.id <- order(sim, decreasing=T)[1:K]
    x.tmp <- ginv(t(x.complete[sim.id, -row.miss, drop=FALSE])) %*%row[-row.miss, drop=FALSE]
    x.Shrink <- as.vector((1 - (K-2)* var(x.tmp)/((length(row) - sum(is.na(row)))*sum(x.tmp**2)))) * x.tmp 
    row[row.miss] <- t(x.complete[sim.id, row.miss, drop=FALSE]) %*% x.Shrink
    return(row)
  }))
  
  xmiss[miss.row, ] <- x.imputed
  
  return(xmiss)
}
