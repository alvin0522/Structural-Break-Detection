####  Version 2 ###############
## Modified at 08/09/2017
## Automatically select number of peaks/valleys

find_peaks <- function(x, span=0.25, lspan=0.05, tol=100){
  n <- length(x)
  y <- x
  mu.y.loc <- y
  mu.y.loc <- (x[1:(n-2)] + x[2:(n-1)] + x[3:n])/3
  mu.y.loc <- c(mu.y.loc[1], mu.y.loc, mu.y.loc[n-2])
  
  y.loess <- loess(x~I(1:n), span=span)
  y <- y.loess[[2]]
  sig.y <- var(y.loess$resid, na.rm=TRUE)^0.5
  DX.1 <- sign(diff(mu.y.loc, na.pad = FALSE))
  pks <- which(diff(DX.1, na.pad = FALSE) < 0 & DX.1[-(n-1)] > 0) + 1
  cpts <- NULL
  cpts.thresh <- NULL
  n.w <- floor(lspan*n/2)
  thresh <- 0.1
  out <- NULL
  for(pk in pks)
  {
    inner <- (pk-n.w):(pk+n.w)
    outer <- c((pk-2*n.w):(pk-n.w),(pk+2*n.w):(pk+n.w))
    if(sum(inner <0) == 0 & sum(outer<0) == 0){
      mu.y.outer <- mean(y[outer])
      if(!is.na(mu.y.outer)) 
        if (mean(y[inner])-mu.y.outer > thresh*sig.y) out <- c(out, pk)
    }
    
  }
  fit <- pamk(out)
  n.cluster <- fit$nc
  grp <- fit$pamobject$clustering
  vec <- split(out, grp)
  elm.n <- as.numeric(sapply(vec,length))
  cpts <- NULL
  for(i in 1:length(elm.n)){
    vec1 <- vec[which(elm.n==elm.n[i])]
    pk <- unlist(vec1)
    pk <- as.numeric(pk)
    thresh <- 0.1
    d.bks <- 1000
    while(d.bks > tol ){
      out <- NULL
      for(j in pk)
      {
        inner <- (j-n.w):(j+n.w)
        outer <- c((j-2*n.w):(j-n.w),(j+2*n.w):(j+n.w))
        if(sum(inner <0) == 0 & sum(outer<0) == 0){
          mu.y.outer <- mean(y[outer])
          if(!is.na(mu.y.outer)) 
            if (mean(y[inner])-mu.y.outer > thresh*sig.y) out <- c(out, j)
        }
        
      }
      
      if(is.null(out)){
        out <- out1
        break
      }
      brks <- pretty(out)
      count <- as.numeric(table(cut(out, pretty(out))))
      pos.set <-  order(count, decreasing = T)
      pos <- pos.set[length(pos.set)]
      lb <- brks[min(pos)]
      ub <- brks[max(pos)+1]
      pks <- pks[which(!(pks %in% seq(lb, ub, 1)))]
      #pks <- out
      thresh <- thresh + 0.1
      out1 <- out
    }
    resd <- y.loess$resid[out]
    cp <- out[which(abs(resd) == max(abs(resd)))]
    cpts <- c(cpts, cp)
}
  return(cpts)
  
}

####  Version 1 ###############

find_peaks2 <- function(x, ncp=1,  span=0.25, lspan=0.05, tol=100){
  n <- length(x)
  y <- x
  mu.y.loc <- y
  mu.y.loc <- (x[1:(n-2)] + x[2:(n-1)] + x[3:n])/3
  mu.y.loc <- c(mu.y.loc[1], mu.y.loc, mu.y.loc[n-2])
  
  y.loess <- loess(x~I(1:n), span=span)
  y <- y.loess[[2]]
  sig.y <- var(y.loess$resid, na.rm=TRUE)^0.5
  DX.1 <- sign(diff(mu.y.loc, na.pad = FALSE))
  pks <- which(diff(DX.1, na.pad = FALSE) < 0 & DX.1[-(n-1)] > 0) + 1
  cpts <- NULL
  cpts.thresh <- NULL
  n.w <- floor(lspan*n/2)
  for(i in 1:ncp){
    
    thresh <- 0.1
    # num <- 1
    d.bks <- 1000
    while(d.bks > tol ){
      out <- NULL
      for(pk in pks)
      {
        inner <- (pk-n.w):(pk+n.w)
        outer <- c((pk-2*n.w):(pk-n.w),(pk+2*n.w):(pk+n.w))
        if(sum(inner <0) == 0 & sum(outer<0) == 0){
          mu.y.outer <- mean(y[outer])
          if(!is.na(mu.y.outer)) 
            if (mean(y[inner])-mu.y.outer > thresh*sig.y) out <- c(out, pk)
        }
        
      }
      
      if(is.null(out)){
        out <- out1
        break
      }
      brks <- pretty(out)
      count <- as.numeric(table(cut(out, pretty(out))))
      pos.set <-  order(count, decreasing = T)
      pos <- pos.set[length(pos.set)]
      lb <- brks[min(pos)]
      ub <- brks[max(pos)+1]
      pks <- pks[which(!(pks %in% seq(lb, ub, 1)))]
      #pks <- out
      thresh <- thresh + 0.1
      out1 <- out
      d.bks <- min(diff(brks))
      count
      brks
    }
    resd <- y.loess$resid[out]
    cp <- out[which(abs(resd) == max(abs(resd)))]
    cpts <- c(cpts, cp)
    cpts.thresh <- c(cpts.thresh, thresh)
    lb <- min(out)
    ub <- max(out)
    if(ub < n/2){
      pks <- pks[which(!(pks %in% seq(1, ub, 1)))]
    }else if(lb > n/2){
      pks <- pks[which(!(pks %in% seq(lb, n, 1)))]
    }else{
      pks <- pks[which(!(pks %in% seq(lb, ub, 1)))]
    }
    
    
  }
  
  return(cpts)
  
}

