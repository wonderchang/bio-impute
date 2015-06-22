
#-----row average impute-----

RowAverage <- function(x){
  mis <- is.na(x)
  rowM <- rowMeans(x, na.rm=T)
  rowM <- matrix(rowM, nc=ncol(x), nr=nrow(x), byrow=F)
  
  x[mis] <- rowM[mis]
  x
}


#-----ILLS impute-----

ILLS <- function(xmiss, K=10, sim.method="EuDist", iter=2){

  xcomplete <- RowAverage(xmiss);
  
  miss.gene <- is.na(xmiss)
  miss.row <- which(rowSums(miss.gene)!=0)
  
  for(h in 1:iter) {
    xcomplete[miss.row, ] <- t(apply(xmiss[miss.row, ], 1, function(row){
      row.miss <- which(is.na(row))
      sim <- similarityCal(row[-row.miss], xcomplete[, -row.miss, drop=F], sim.method)
      sim.id <- order(sim, decreasing=T)[2:(K+1)]
      row[row.miss] <- t(xcomplete[sim.id, row.miss, drop=FALSE]) %*% ginv(t(xcomplete[sim.id, -row.miss, drop=FALSE])) %*%row[-row.miss, drop=FALSE]
      return(row)
    }))
  }
  return(xcomplete) 
}

#-----caculate distance 

similarityCal<-function(vec, mat, method="EuDist"){
  methods<-c("EuDist","cor","Angle")
  switch(match.arg(method,methods),
         EuDist=1/sqrt(rowSums((mat-matrix(vec,nc=length(vec),nr=nrow(mat),byrow=TRUE))^2)),
         cor=apply(mat,1,function(i) abs(cor(i,vec,use="everything",method="pearson"))),
         Angle=apply(mat,1,function(i) abs(sum(i * vec)/sqrt(sum(i^2)*sum(vec^2))))
  )
}
