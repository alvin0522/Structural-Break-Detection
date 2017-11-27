if("moments" %in% rownames(installed.packages()) == FALSE)
  install.packages("moments")
if("MASS" %in% rownames(installed.packages()) == FALSE)
  install.packages("MASS")
library("moments")
library("MASS")

## Psi Calculation ##
psi_cal <- function(x,est,p,q){
  mxpq <- max(p,q)
  nx <- length(x)
  psih <- rep(1,nx)
  pdim <- p+q+1
  omega <- as.matrix(est[1],ncol=1,nrow=1)
  alpha <- as.matrix(est[2:(1+p)],ncol=1,nrow=p)
  if(q==0){
    beta <- 0
  }else
  {
    beta <- as.matrix(est[(p+2):pdim],ncol=1,nrow=q) 
  }
  
  for( t in 1:mxpq){
    psih[t] <- est[1]
  }
  for( t in (mxpq+1):nx){
    pos1 <- t-1
    pos2 <- t-p
    pos3 <- t-q
    psih[t]=psit(x, psih, omega, alpha, beta, pos1, pos2, pos3, q)
  }
  return(psih)
}

## First Four Moments Calculation ##
eps_dist2 <- function(eps){
  
  mts <- all.moments(eps, order.max=4, central = TRUE)
  mue <- mean(eps)
  vare <- var(eps)
  skewe <- mts[4]/vare^(3/2)
  kurte <- mts[5]/vare^2
  dist1 <- "empirical"
  
  return(list(mue = mue, vare=vare, skewe=skewe, kurte=kurte, dist1=dist1))
}

