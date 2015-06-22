# Similarity method
similarityCal<-function(vec, mat, method="EuDist"){
  methods<-c("EuDist","cor","Angle")
  switch(match.arg(method,methods),
         EuDist=1/sqrt(rowSums((mat-matrix(vec,nc=length(vec),nr=nrow(mat),byrow=TRUE))^2)),
         cor=apply(mat,1,function(i) abs(cor(i,vec,use="everything",method="pearson"))),
         Angle=apply(mat,1,function(i) abs(sum(i * vec)/sqrt(sum(i^2)*sum(vec^2))))
  )
}

# 1. Zero
Zero <- function(x){
  x[is.na(x)] <- 0
  x
}

# 2. RowAverage
RowAverage <- function(x){
  mis <- is.na(x)
  rowM <- rowMeans(x, na.rm=T)
  rowM <- matrix(rowM, nc=ncol(x), nr=nrow(x), byrow=F)
  
  x[mis] <- rowM[mis]
  x
}

# 3. KNN
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

# 4. SKNN
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

# 5. IKNN
IKNN <- function(xmiss, K=10, sim.method="EuDist", iter=2){
  
  xcomplete <- RowAverage(xmiss)
  
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

# 6. ISKNN
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

# 7. LS
LS <- function(xmiss, K=10, sim.method="EuDist"){
  miss.gene <- is.na(xmiss)
  miss.row <- which(rowSums(miss.gene)!=0)
  
  xcomplete <- xmiss[-miss.row, ]
  xincomplete <- xmiss[miss.row, ]
  
  impute <- function(row) {
    
    row.exp <- which(is.na(row))
    gene <- row[-row.exp]
    cand_x <- xcomplete[, -row.exp, drop=F]
    
    sim <- similarityCal(gene,cand_x, method=sim.method)
    row.idx <- order(sim, decreasing=T)[1:K]
    row.r <- sim[row.idx]
    row.cand <- cand_x[row.idx, , drop=F]
    lg <- apply(row.cand, 1, function(x){lm(gene~x)$coefficients})
    row.impcand <- xcomplete[row.idx, row.exp, drop=F]
    
    y <- matrix(0, nc=ncol(row.impcand), nr=nrow(row.impcand))
    w <- (row.r**2/(1-row.r**2+0.000001))**2
    sw <- sum(w)
    w <- w/sw
    
    for(i in 1:nrow(row.impcand)) {
      y[i, ] <- lg[2, i] * row.impcand[i, ] + lg[1, i]
    }
    row[row.exp]<-apply(y, 2, function(x){sum(w*x)})
    return(row)
  }
  
  xmiss[miss.row, ] <- t(apply(xincomplete, 1, impute))
  
  return (xmiss)
}

# 8. LLS
LLS <- function(xmiss, K=10, sim.method="EuDist"){
  
  miss.gene <- is.na(xmiss)
  miss.row <- which(rowSums(miss.gene)!=0)
  x.complete <- xmiss[-miss.row, ]
  x.incomplete <- xmiss[miss.row, ]
  
  x.imputed <- t(apply(x.incomplete, 1, function(row){
    row.miss <- which(is.na(row))
    sim <- similarityCal(row[-row.miss], x.complete[, -row.miss], sim.method)
    sim.id <- order(sim, decreasing=T)[1:K]
    row[row.miss] <-t(x.complete[sim.id, row.miss, drop=FALSE]) %*% ginv(t(x.complete[sim.id, -row.miss, drop=FALSE])) %*%row[-row.miss, drop=FALSE]
    return(row)
  }))
  
  xmiss[miss.row, ] <- x.imputed
  
  return(xmiss)
}

# 9. SLLS
SLLS <- function(xmiss, K=10, sim.method="EuDist"){
  
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
    row[row.miss] <- t(x.complete[sim.id, row.miss, drop=FALSE]) %*% ginv(t(x.complete[sim.id, -row.miss, drop=FALSE])) %*%row[-row.miss, drop=FALSE]
    if(length(row.miss) < therhold){
      x.complete <<- rbind(x.complete, row)
    }
    return(row)
  }))
  
  xmiss[miss.row[miss.rate.order], ] <- x.imputed
  
  return(xmiss)
}

# 10. WLLS
WLLS <- function(xmiss, K=10, sim.method="EuDist", order=2){

  miss.gene <- is.na(xmiss)
  miss.row <- which(rowSums(miss.gene)!=0)
  x.complete <- xmiss[-miss.row, ]
  x.incomplete <- xmiss[miss.row, ]
  
  x.imputed <- t(apply(x.incomplete, 1, function(row){
    row.miss <- which(is.na(row))
    sim <- similarityCal(row[-row.miss], x.complete[, -row.miss], sim.method)
    sim.id <- order(sim, decreasing=T)[1:K]
    x.diag <- diag(sim[sim.id]**order)
    row[row.miss] <- t(x.diag %*% x.complete[sim.id, row.miss, drop=FALSE]) %*% ginv(t(x.diag %*% x.complete[sim.id, -row.miss, drop=FALSE])) %*% row[-row.miss, drop=FALSE]
    return(row)
  }))
  
  xmiss[miss.row, ] <- x.imputed
  
  return(xmiss)
}

# 11. ShrLLS
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

# 12. ISLLS
ISLLS <- function(xmiss, K=10, sim.method="EuDist", iter=2){
  
  xcomplete <- SLLS(xmiss, K, sim.method);
  
  miss.gene <- is.na(xmiss)
  miss.row <- which(rowSums(miss.gene)!=0)
  
  for(h in 1:iter) {
    xcomplete[miss.row, ] <- t(apply(xmiss[miss.row, ], 1, function(row){
      row.miss <- which(is.na(row))
      sim <- similarityCal(row[-row.miss], xcomplete[, -row.miss], sim.method)
      sim.id <- order(sim, decreasing=T)[2:(K+1)]
      row[row.miss] <- t(xcomplete[sim.id, row.miss, drop=FALSE]) %*% ginv(t(xcomplete[sim.id, -row.miss, drop=FALSE])) %*%row[-row.miss, drop=FALSE]
      return(row)
    }))
  }
  return(xcomplete) 
}

# 13. WSLLS
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

# 14. ShrSLLS
ShrSLLS <- function(xmiss, K=10, sim.method="EuDist"){

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
    x.tmp <- ginv(t(x.complete[sim.id, -row.miss, drop=FALSE])) %*%row[-row.miss, drop=FALSE]
    x.Shrink <- as.vector((1 - (K-2)* var(x.tmp)/((length(row) - sum(is.na(row)))*sum(x.tmp**2)))) * x.tmp 
    row[row.miss] <- t(x.complete[sim.id, row.miss, drop=FALSE]) %*% x.Shrink
    if(length(row.miss) <= therhold){
      x.complete <<- rbind(x.complete, row)
    }
    return(row)
  }))
  
  xmiss[miss.row[miss.rate.order], ] <- x.imputed
  
  return(xmiss)
}

# 15. IWLLS
IWLLS <- function(xmiss, K=10, sim.method="EuDist", iter=2, order=2){
  
  xcomplete <- WLLS(xmiss, K, sim.method, order)

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

# 16. ILLS
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

# 17. IShrLLS
IShrLLS <- function(xmiss, K=10, sim.method="EuDist", iter=2){
  
  xcomplete <- ShrLLS(xmiss, K, sim.method)

  miss.gene <- is.na(xmiss)
  miss.row <- which(rowSums(miss.gene)!=0)

  for(h in 1:iter) {
    xcomplete[miss.row, ] <- t(apply(xmiss[miss.row, ], 1, function(row){
      row.miss <- which(is.na(row))
      sim <- similarityCal(row[-row.miss], xcomplete[, -row.miss], sim.method)
      sim.id <- order(sim, decreasing=T)[2:(K+1)]
      x.tmp <- ginv(t(xcomplete[sim.id, -row.miss, drop=FALSE])) %*%row[-row.miss, drop=FALSE]
      x.Shrink <- as.vector((1 - (K-2)* var(x.tmp)/((length(row) - sum(is.na(row)))*sum(x.tmp**2)))) * x.tmp 
      row[row.miss] <- t(xcomplete[sim.id, row.miss, drop=FALSE]) %*% x.Shrink
      return(row)
    }))
  }
  return(xcomplete) 
}

# 18. ShrWLLS
ShrWLLS <- function(xmiss, K=10, sim.method="EuDist", order=2){
  
  miss.gene <- is.na(xmiss)
  miss.row <- which(rowSums(miss.gene)!=0)
  x.complete <- xmiss[-miss.row, ]
  x.incomplete <- xmiss[miss.row, ]
  
  x.imputed <- t(apply(x.incomplete, 1, function(row){
    row.miss <- which(is.na(row))
    sim <- similarityCal(row[-row.miss], x.complete[, -row.miss], sim.method)
    sim.id <- order(sim, decreasing=T)[1:K]
    x.diag <- diag(sim[sim.id]**order)
    x.tmp <- ginv(t(x.diag %*% x.complete[sim.id, -row.miss, drop=FALSE])) %*%row[-row.miss, drop=FALSE]
    x.Shrink <- as.vector((1 - (K-2)* var(x.tmp)/((length(row) - sum(is.na(row)))*sum(x.tmp**2)))) * x.tmp 
    row[row.miss] <- t(x.diag %*% x.complete[sim.id, row.miss, drop=FALSE]) %*% x.Shrink
    return(row)
  }))
  
  xmiss[miss.row, ] <- x.imputed
  
  return(xmiss)
}

# 19. IWSLLS
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

# 20. IShrSLLS
IShrSLLS <- function(xmiss, K=10, sim.method="EuDist", iter=2){

  xcomplete <- ShrSLLS(xmiss, K, sim.method)

  miss.gene <- is.na(xmiss)
  miss.row <- which(rowSums(miss.gene)!=0)

  for(h in 1:iter) {
    xcomplete[miss.row, ] <- t(apply(xmiss[miss.row, ], 1, function(row){
      row.miss <- which(is.na(row))
      sim <- similarityCal(row[-row.miss], xcomplete[, -row.miss], sim.method)
      sim.id <- order(sim, decreasing=T)[2:(K+1)]
      x.tmp <- ginv(t(xcomplete[sim.id, -row.miss, drop=FALSE])) %*%row[-row.miss, drop=FALSE]
      x.Shrink <- as.vector((1 - (K-2)* var(x.tmp)/((length(row) - sum(is.na(row)))*sum(x.tmp**2)))) * x.tmp 
      row[row.miss] <- t(xcomplete[sim.id, row.miss, drop=FALSE]) %*% x.Shrink
      return(row)
    }))
  }
  return(xcomplete) 
}

# 21. IShrWLLS
IShrWLLS <- function(xmiss, K=10, sim.method="EuDist", iter=2, order=2) {

  xcomplete <- ShrWLLS(xmiss, K, sim.method, order)

  miss.gene <- is.na(xmiss)
  miss.row <- which(rowSums(miss.gene)!=0)

  for(h in 1:iter) {
    xcomplete[miss.row, ] <- t(apply(xmiss[miss.row, ], 1, function(row){
      row.miss <- which(is.na(row))
      sim <- similarityCal(row[-row.miss], xcomplete[, -row.miss], sim.method)
      sim.id <- order(sim, decreasing=T)[2:(K+1)]
      x.diag <- diag(sim[sim.id]**order)
      x.tmp <- ginv(t(x.diag %*% xcomplete[sim.id, -row.miss, drop=FALSE])) %*%row[-row.miss, drop=FALSE]
      x.Shrink <- as.vector((1 - (K-2)* var(x.tmp)/((length(row) - sum(is.na(row)))*sum(x.tmp**2)))) * x.tmp 
      row[row.miss] <- t(x.diag %*% xcomplete[sim.id, row.miss, drop=FALSE]) %*% x.Shrink
      return(row)
    }))
  }
  return(xcomplete) 
}

# 22. ShrWSLLS
ShrWSLLS <- function(xmiss, K=10, sim.method="EuDist", order=2){
  
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
    x.tmp <- ginv(t(x.diag %*% x.complete[sim.id, -row.miss, drop=FALSE])) %*% row[-row.miss, drop=FALSE]
    x.Shrink <- as.vector((1 - (K-2)* var(x.tmp)/((length(row) - sum(is.na(row)))*sum(x.tmp**2)))) * x.tmp 
    row[row.miss] <- t(x.diag %*% x.complete[sim.id, row.miss, drop=FALSE]) %*% x.Shrink

    if(length(row.miss) <= therhold){
      x.complete <<- rbind(x.complete, row)
    }
    return(row)
  }))

  
  xmiss[miss.row[miss.rate.order], ] <- x.imputed
  
  return(xmiss)
}

# 23. IShrWSLLS
IShrWSLLS <- function(xmiss, K=10, sim.method="EuDist", iter=2, order=2) {

  xcomplete <- ShrWSLLS(xmiss, K, sim.method, order)

  miss.gene <- is.na(xmiss)
  miss.row <- which(rowSums(miss.gene)!=0)

  for(h in 1:iter) {
    xcomplete[miss.row, ] <- t(apply(xmiss[miss.row, ], 1, function(row){
      row.miss <- which(is.na(row))
      sim <- similarityCal(row[-row.miss], xcomplete[, -row.miss], sim.method)
      sim.id <- order(sim, decreasing=T)[2:(K+1)]
      x.diag <- diag(sim[sim.id]**order)
      x.tmp <- ginv(t(x.diag %*% xcomplete[sim.id, -row.miss, drop=FALSE])) %*%row[-row.miss, drop=FALSE]
      x.Shrink <- as.vector((1 - (K-2)* var(x.tmp)/((length(row) - sum(is.na(row)))*sum(x.tmp**2)))) * x.tmp 
      row[row.miss] <- t(x.diag %*% xcomplete[sim.id, row.miss, drop=FALSE]) %*% x.Shrink
      return(row)
    }))
  }
  return(xcomplete) 
}

