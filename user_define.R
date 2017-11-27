## definiton of psih
psit <- function(x, psih, omega, alpha, beta, pos1, pos2, pos3, q){
  lagged.logx <- as.matrix(c(log(x[pos1:pos2])))
  lagged.psih <- as.matrix(c(psih[pos1:pos3]))
  if(q == 0) {
    psih=omega + t(alpha)%*%lagged.logx
  }else{
    psih=omega + t(alpha)%*%lagged.logx+t(beta)%*%lagged.psih
  }
  return(psih)  
}


## first derivative of psih
derpsit <- function(x, psih, pos1, pos2, pos3, pdim, q){
  if(q == 0){
    derpsi=matrix(c(1,log(x[pos1:pos2])),pdim,1) 
  }else{
    derpsi=matrix(c(1,log(x[pos1:pos2]),psih[pos1:pos3]),pdim,1) 
  }
  return(derpsi)
}

## second derivative of psih
der2psit <- function(pdim){
  der2psi=matrix(rep(0),pdim,pdim)  
}

##  Initial Value for K-inverse
InitMat <- function(x,p,q){
  y <- log(x)
  pdim <- p+q+1
  arma.model <- arima(y, order=c(max(p,q),0,q), include.mean=TRUE,method="CSS")
  Sigma <- arma.model$var.coef
  arma.phi <- as.numeric(arma.model$coef[1:p])
  arma.th <- as.numeric(arma.model$coef[(p+1):(pdim-1)])
  arma.mu <- as.numeric(arma.model$coef[pdim])
  
  #if(p==1 & q==0) {
  #  mat.vec <- c(-arma.mu, 1-arma.phi,
  #               1, 0)
  #}
  if(q == 0){
    R1 <- c(rep(-arma.mu, p), 1-sum(arma.phi))
    #R2 <- matrix(R1, nrow = 1, ncol= (p+length(arma.phi)))
    M1 <- diag(1, p)
    M1 <- cbind(M1, matrix(0, p, ncol = 1))
    mat.vec <- rbind(R1, M1)
    mat.vec <- as.matrix(mat.vec)
  }else{
    if(p==q){
      R1 <- c(rep(-arma.mu,p),rep(0,q), 1-sum(arma.phi))
      B1 <- diag(1,p,p)
      B2 <- cbind(B1, matrix(0,p,1))
      B3 <- diag(0,q,p)
      B4 <- cbind(diag(-1,q,q),matrix(0,q,1))
      mat.vec <- rbind(R1, cbind(B1,B2), cbind(B3,B4))
      mat.vec <- as.matrix(mat.vec)
    }else{
      R1 <- c(rep(-arma.mu,p),rep(0,q), 1-sum(arma.phi))
      B1 <- diag(1,p,p)
      B2 <- B1
      diag(B2)[p] <- 0
      B3 <- diag(0,q,p)
      B4 <- cbind(diag(-1,q,q),matrix(0,q,1))
      mat.vec <- rbind(R1,cbind(B1,B2), cbind(B3,B4))
      mat.vec <- as.matrix(mat.vec)
    }
  }
  D <- matrix(mat.vec, pdim, pdim, byrow=T)
  init.mat <- D%*%Sigma%*%t(D)
  init.mat <- solve(init.mat)
  return(init.mat)
}

## Initial Value Estimation
finitval <- function (x, p=1, q=1) {
    
    ## Take the log transform of the calibration data and define it as y
    y <- log(x)
    
    ## Fit an ARMA(m,n) model to the log transformed durations
    arma.model <- arima(y, order=c(max(p,q),0,q), include.mean=TRUE,method="CSS")
    
    ## OMEGA: omega vector is the intercept
    ## omega.hat <- ar.model$coef[length(ar.model$coef)]
    MEAN <- arma.model$coef[length(arma.model$coef)]
    intercept <- MEAN*(1-sum(arma.model$coef[1:p]))
    omega.hat <- intercept
    
    ## BETA: beta vector have the MA coefficients
    if(q==0) beta.hat <- 0
    if(q>0) beta.hat <- -1*arma.model$coef[(length(arma.model$coef)-q):(length(arma.model$coef)-1)]
    
    ## ALPHA: phi vector equals AR coefficients
    phi.hat <- arma.model$coef[1:max(p,q)]
    ## theta vector has length equal to that of the phi vector, and values equal to the 
    ## MA coefficients (betas) for the first q entries
    theta.hat <- rep(0, max(p,q))
    for(i in 1:length(beta.hat)) { theta.hat[i] <- beta.hat[i] }
    ## alpha vector equals phi-theta
    alpha.hat <- phi.hat-theta.hat		
    
    
    ### Var-cov matrix of estimates
    varcov <-vcov(arma.model)
    #vec.varcov <-c(varcov)
    
    ## Results
    results <- c(omega.hat, alpha.hat, beta.hat)
  }

## first derivative of SCAD penalty
p.lam.prime <- function(theta, initial, lambda, a,tau){
  n <- length(initial)
  penalty <- rep(0, n)
  for(i in 1:n){
    t1 <- abs(initial[i]) <= lambda
    t2 <- max((a*lambda- abs(initial[i])),0)/((a-1)*lambda)
    t3 <- abs(initial[i]) > lambda
    penalty[i] <- lambda*(t1+ t2*t3)/abs(initial[i]+tau)*theta[i]
  }
  return(penalty)
}

## second derivative of SCAD penalty
p.lam.2prime <- function(initial,lambda, a,tau){
  n <- length(initial)
  penalty <- rep(0, n)
  for(i in 1:n){
    t1 <- abs(initial[i]) <= lambda
    t2 <- max((a*lambda- abs(initial[i])),0)/((a-1)*lambda)
    t3 <- abs(initial[i]) > lambda
    penalty[i] <- lambda*(t1+ t2*t3)/(abs(initial[i])+tau)
    #penalty[i] <- lambda*(t1+ t2*t3)
  }
  return(penalty)
}
