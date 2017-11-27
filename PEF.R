## Method Average Vector Recursive Estimation (AVRE) ##
# Inputs are:
# x: durations (vector)
# initval: initial values for Log ACD model
# p: LogACD model, first parameter
# q: LogACD model, second parameter
# momente: first four moments (vector)

# Outputs are:
# estimation (vector)

recursive_mat <- function(x,initval,p,q,momente) {
  
  n <- length(x)
  pdim <- 1+p+q   # number of parameters to be estimated
  
  mue <- momente[1]
  vare <- momente[2]
  skewe <- momente[3]
  kurte <- momente[4]
  
  ########################################################################  	
  ## INITIALIZATION OF ARRAYS
  ########################################################################
  ## Identity matrix
  iden = diag(pdim)
  
  ## psi hat
  psih <- rep(1,n)
  
  ## k matrix (variance-covariance) and k inverse (observed information)
  kmat = array(NA, dim = c(pdim, pdim, n)) 
  kinv = array(NA, dim = c(pdim, pdim, n))
  termk1 = array(NA, dim=c(pdim,pdim,n))
  termk2 = array(NA, dim=c(pdim,pdim,n))
  
  ## Parameter estimates for each iteration
  thehat = array(NA, dim = c(pdim, 1, n))
  
  ## Derivative of psi and second derivative of psi
  derpsi<-matrix(rep(0),pdim,1)
  der2psi<-matrix(rep(0),pdim,pdim)
  
  ## Derivative of mu, sigsq, gamma, kappa; second derivates of mu, sigsq				
  dermu<-matrix(rep(0),pdim,1)				
  dersigsq<-matrix(rep(0),pdim,1)
  dergamma<-matrix(rep(0),pdim,1)
  derkappa<-matrix(rep(0),pdim,1)				
  der2mu<-matrix(rep(0),pdim,pdim)				
  der2sigsq<-matrix(rep(0),pdim,pdim)
  
  ## Derivative of m(t) and M(t) 
  derm<-matrix(rep(0),pdim,1)
  derqm<-matrix(rep(0),pdim,1)
  
  ## Derivative of quadratic variation, quadratic covariation, eta, rho
  dervm<-matrix(rep(0),pdim,1)
  dervqm<-matrix(rep(0),pdim,1)
  dervmqm<-matrix(rep(0),pdim,1)
  dereta<-matrix(rep(0),pdim,1)
  derrho<-matrix(rep(0),pdim,1)
  
  
  ## Optimal astr and bstr
  astr<-matrix(rep(0),pdim,1)
  bstr<-matrix(rep(0),pdim,1)
  
  ## Derivative of astr and bstr
  derastr<-matrix(rep(0),pdim,pdim)
  derbstr<-matrix(rep(0),pdim,pdim)
  
  
  ########################################################################
  ## INITIAL VALUES: for omega, alpha, beta, passed from main program
  ########################################################################
  
  #initial<- as.numeric(c(initval[1], initval[2], initval[3]))
  initial <- initval
  
  ########################################################################
  ## RECURSIVE FORMULAS FROM ESTIMATING EQUATIONS
  ########################################################################
  ## Put initial values into initial positions of arrays
  thehat[,,1]=initial
  
  
  
  init.mat <- InitMat(x,p,q)
  
  #init.mat <- diag(diag(init.mat))
  
  mxpq <- max(p,q)
  for (t in 1:mxpq){
    psih[t] = thehat[1,1,1]    # omega 
    thehat[,1,t] <- thehat[,1,1]
    kinv[,,t] <- init.mat     ###### KEY PART #############
    kmat[,,t] <- solve(kinv[,,t])
  }
  
  
  
  ## t=(mxpq+1):n					
  for (t in (mxpq+1):n)
  {			
    omega <- as.matrix(thehat[1,1,t-1],ncol=1,nrow=1)
    alpha <- as.matrix(thehat[2:(1+p),1,t-1],ncol=1,nrow=p)
    if(q != 0){
      beta <- as.matrix(thehat[(p+2):pdim,1,t-1],ncol=1,nrow=q)
    }
    pos1 <- t-1
    pos2 <- t-p
    pos3 <- t-q
    
    psih[t]=psit(x,psih,omega, alpha, beta, pos1, pos2, pos3, q)
    if(abs(psih[t]) > log(.Machine$double.xmax/max(momente))/4) {
      psih[t] <- sign(psih[t])*log(.Machine$double.xmax/max(momente))/8
    }
    
    ## Define derivatives of psih(t) wrt theta: pdim*1 vector							
    ## First derivative									
    derpsi=derpsit(x, psih, pos1, pos2, pos3, pdim, q)
    
    
    
    ## Second derivative: pdim*pdim matrix									
    der2psi=der2psit(pdim)   
    
    ## mu(t), sigsq(t), gamma(t), kappa(t)							
    mu=mue*exp(psih[t]) 							
    sigsq=vare*exp(2*psih[t])				
    gamma=skewe*exp(3*psih[t]) # recall this is third central moment							
    kappa=kurte*exp(4*psih[t]) # recall this is fourth central moment
    
    ## First Derivatives of mu(t),sigsq(t), gamma(t) and kappa(t)												
    dermu=mue*exp(psih[t])*derpsi							
    dersigsq=2*vare*exp(2*psih[t])*derpsi
    dergamma=3*skewe*exp(3*psih[t])*derpsi
    derkappa=4*kurte*exp(4*psih[t])*derpsi
    
    ## Second Derivatives of mu(t) and sigsq(t)							
    der2mu <- mu*((derpsi)%*%t(derpsi) + der2psi)
    der2sigsq <- 2*sigsq*(der2psi + 2*(derpsi)%*%t(derpsi))
    
    ## Compute m(t) and M(t)							
    m = x[t]-mu       							
    qm = m**2-sigsq
    
    ## Compute Quadratic variations of m(t) and M(t) and 
    ## covariance of (m(t), M(t))							
    vm = sigsq
    vqm = kappa-vm**2  
    vmqm = gamma
    
    ## Define rho^2(t) 															
    termr = 1-((vmqm**2)/(vm*vqm))								
    rho = 1/termr
    if(is.na(rho) | is.infinite(rho))
      rho <- 1
    
    ## Define eta 								
    eta = vmqm/(vm*vqm)	
    if(is.na(eta) | is.infinite((eta)))
      eta = 0
    
    ## Define vectors astr and bstr						
    astr = rho*(-dermu/vm + dersigsq*eta)							
    bstr = rho*(dermu*eta - dersigsq/vqm)
    
    ## Define Derivatives of m(t) and M(t) 							
    derm = -dermu
    derqm = 2*m*derm - dersigsq	
    
    ## Derivatives of variations <m>(t), <M>(t) and <m,M>(t)							
    dervm = dersigsq						
    dervqm = derkappa - 2*sigsq*dersigsq
    dervmqm = dergamma
    
    #Note: derrho=0 - and it is!
    ru=vm*vqm
    rv=vm*vqm-vmqm**2
    rdu=vm*dervqm+vqm*dervm
    rdv=vm*dervqm+vqm*dervm -2*vmqm*dervmqm
    derrho=(rv*rdu-ru*rdv)/(rv**2)
    if(sum(is.na(derrho)) >0)
      derrho <- matrix(rep(0),pdim,1)
    ## Derivative of eta(t)							
    num1=vm*vqm*dervmqm
    num2=vmqm*(vm*dervqm + vqm*dervm)
    den=(vm**2)*(vqm**2)
    dereta=(num1-num2)/den
    if(sum(is.infinite(dereta)) >0 | sum(is.na(dereta))>0)
      dereta <- matrix(rep(0),pdim,1)
    ## Derivatives of astr and bstr							
    ## For astr						
    terma1 = der2mu/vm - dermu%*%t(dervm)/(vm**2)
    terma2 = der2sigsq*eta + dersigsq%*%t(dereta)
    terma3 = -dermu/vm + dersigsq*eta							
    derastr = -rho*terma1 + rho*terma2 + derrho%*%t(terma3)
    #derastr = -rho*terma1 + rho*terma2 
    
    ## For bstr						
    termb1 = der2mu*eta + dermu%*%t(dereta)													
    termb21 = der2sigsq/vqm 
    termb22 = (dersigsq%*%t(dervqm))/(vqm**2)
    if(sum(is.infinite(termb22)>0 | sum(is.na(termb22))>0))
      termb22 <- matrix(0, ncol=pdim, nrow = pdim)
    termb2 = termb21 - termb22
    termb3 = dermu*eta - dersigsq/vqm
    derbstr = rho*termb1 - rho*termb2 + derrho%*%t(termb3)
    #derbstr = rho*termb1 - rho*termb2
    
    ### Recursive Formulas						
    ## Compute Kinv(t): Information 							
    termk1[,,t] = astr%*%t(derm) + m * derastr 							
    termk2[,,t] = bstr%*%t(derqm) + qm * derbstr
    
    
    ## Need Numerical FIX when termk1 and termk2 are NA
    if ( sum(as.numeric(is.na(termk1[,,t]))) > 0){
      foundFlag1 = 0
      for (w in 1:(t-1)) {
        if ( (foundFlag1 == 0) & sum(as.numeric(is.na(termk1[,,t-w]))) == 0 ) {
          termk1[,,t] = termk1[,,(t-w)]
          foundFlag1 = 1
        }
        else
          termk1[,,t]= matrix(rep(.Machine$double.xmax/3,pdim*pdim),pdim,pdim)
      }
    }
    
    if ( sum(as.numeric(is.na(termk2[,,t]))) > 0){
      foundFlag2 = 0
      for (w in 1:(t-1)) {
        if ( (foundFlag2 == 0) & sum(as.numeric(is.na(termk2[,,t-w]))) == 0 ) {
          termk2[,,t] = termk2[,,(t-w)]
          foundFlag2 = 1
        }
        else
          termk2[,,t]= matrix(rep(.Machine$double.xmax/3,pdim*pdim),pdim,pdim)
      }
    }
    kinv[,,t] = kinv[,,t-1] - (termk1[,,t]+termk2[,,t])
    
    s <- svd(kinv[,,t])
    inf.index <- which(is.infinite(s$d))
    D <- diag(1/s$d)
    for(i in 1:length(s$d)){
      if(is.infinite(D[i,i])) 
        D[i,i] <- .Machine$double.xmin*100
    }
    kmat[,,t] <- s$v%*%D%*%t(s$u)
    
    
    
    
    ## compute thehat[t]						
    termt = astr*m + bstr*qm
    if(sum(as.numeric(is.na(termt))) >0 ) {
      #cat("NA occurs at time:",t,"\n")
      break
    }
    
    thehat[,,t] = thehat[,,t-1]+kmat[,,t]%*%termt
  }
  
  
  
  finalest <- return(est = thehat[,1,])
  
}



## Method Average Vector Recursive Estimation (AVRE - Penalized version) ##
# Inputs are:
# x: durations (vector)
# initval: initial values for Log ACD model
# p: LogACD model, first parameter
# q: LogACD model, second parameter
# momente: first four moments (vector)

# Outputs are:
# estimation (vector)

recursive_mat2 <- function(x,initval,p,q,momente,lambda,a,tau=0) {
  
  n <- length(x)
  pdim <- 1+p+q   # number of parameters to be estimated
  #estimate <- matrix(0, ncol = pdim, nrow = iteration)
  err <- rep(1,pdim)
  mue <- momente[1]
  vare <- momente[2]
  skewe <- momente[3]
  kurte <- momente[4]
  
  ########################################################################  	
  ## INITIALIZATION OF ARRAYS
  ########################################################################
  ## Identity matrix
  iden = diag(pdim)
  
  ## psi hat
  psih <- rep(1,n)
  
  ## k matrix (variance-covariance) and k inverse (observed information)
  kmat = array(NA, dim = c(pdim, pdim, n)) 
  kinv = array(NA, dim = c(pdim, pdim, n))
  kmat2 = array(NA, dim = c(pdim, pdim, n)) 
  kinv2 = array(NA, dim = c(pdim, pdim, n))
  termk1 = array(NA, dim=c(pdim,pdim,n))
  termk2 = array(NA, dim=c(pdim,pdim,n))
  
  ## Parameter estimates for each iteration
  thehat = array(NA, dim = c(pdim, 1, n))
  
  ## Derivative of psi and second derivative of psi
  derpsi<-matrix(rep(0),pdim,1)
  der2psi<-matrix(rep(0),pdim,pdim)
  
  ## Derivative of mu, sigsq, gamma, kappa; second derivates of mu, sigsq				
  dermu<-matrix(rep(0),pdim,1)				
  dersigsq<-matrix(rep(0),pdim,1)
  dergamma<-matrix(rep(0),pdim,1)
  derkappa<-matrix(rep(0),pdim,1)				
  der2mu<-matrix(rep(0),pdim,pdim)				
  der2sigsq<-matrix(rep(0),pdim,pdim)
  
  ## Derivative of m(t) and M(t) 
  derm<-matrix(rep(0),pdim,1)
  derqm<-matrix(rep(0),pdim,1)
  
  ## Derivative of quadratic variation, quadratic covariation, eta, rho
  dervm<-matrix(rep(0),pdim,1)
  dervqm<-matrix(rep(0),pdim,1)
  dervmqm<-matrix(rep(0),pdim,1)
  dereta<-matrix(rep(0),pdim,1)
  derrho<-matrix(rep(0),pdim,1)
  
  
  ## Optimal astr and bstr
  astr<-matrix(rep(0),pdim,1)
  bstr<-matrix(rep(0),pdim,1)
  
  ## Derivative of astr and bstr
  derastr<-matrix(rep(0),pdim,pdim)
  derbstr<-matrix(rep(0),pdim,pdim)
  
  
  ########################################################################
  ## INITIAL VALUES: for omega, alpha, beta, passed from main program
  ########################################################################
  
  #initial<- as.numeric(c(initval[1], initval[2], initval[3]))
  initial <- initval
  
  ########################################################################
  ## RECURSIVE FORMULAS FROM ESTIMATING EQUATIONS
  ########################################################################
  ## Put initial values into initial positions of arrays
  
  thehat[,,1] <- initial
  init.mat <- InitMat(x,p,q)
  
  sum_termk <- 0
  mxpq <- max(p,q)
  for (t in 1:mxpq){
    psih[t] = thehat[1,1,1]    # omega 
    thehat[,1,t] <- thehat[,1,1]
    kinv[,,t] <- init.mat     ###### KEY PART #############
    #kinv2[,,t] <- init.mat
    kmat[,,t] <- solve(kinv[,,t])
    #kmat2[,,t] <- solve(kinv2[,,t])
  }
  
  pp <- NULL
  
  ## t=(mxpq+1):n					
  for (t in (mxpq+1):n)
  {			
    omega <- as.matrix(thehat[1,1,t-1],ncol=1,nrow=1)
    alpha <- as.matrix(thehat[2:(1+p),1,t-1],ncol=1,nrow=p)
    if(q != 0){
      beta <- as.matrix(thehat[(p+2):pdim,1,t-1],ncol=1,nrow=q)
    }
    pos1 <- t-1
    pos2 <- t-p
    pos3 <- t-q
    
    psih[t]=psit(x,psih,omega, alpha, beta, pos1, pos2, pos3, q)
    if(abs(psih[t]) > log(.Machine$double.xmax/max(momente))/4) {
      psih[t] <- sign(psih[t])*log(.Machine$double.xmax/max(momente))/8
    }
    
    ## Define derivatives of psih(t) wrt theta: pdim*1 vector							
    ## First derivative									
    derpsi=derpsit(x, psih, pos1, pos2, pos3, pdim, q)
    
    
    
    ## Second derivative: pdim*pdim matrix									
    der2psi=der2psit(pdim)   
    
    ## mu(t), sigsq(t), gamma(t), kappa(t)							
    mu=exp(psih[t])*mue 							
    sigsq=vare*exp(2*psih[t])			
    gamma=skewe*exp(3*psih[t])	 # recall this is third central moment							
    kappa=kurte*exp(4*psih[t])	 # recall this is fourth central moment
    
    ## First Derivatives of mu(t),sigsq(t), gamma(t) and kappa(t)												
    dermu=exp(psih[t])*derpsi*mue							
    dersigsq=2*vare*exp(2*psih[t])*derpsi
    dergamma=3*skewe*exp(3*psih[t])*derpsi	
    derkappa=4*kurte*exp(4*psih[t])*derpsi
    
    ## Second Derivatives of mu(t) and sigsq(t)							
    der2mu <- mu*((derpsi)%*%t(derpsi) + der2psi)
    der2sigsq <- 2*sigsq*(der2psi + 2*(derpsi)%*%t(derpsi))
    
    ## Compute m(t) and M(t)							
    m = x[t]-mu       							
    qm = m**2-sigsq
    
    ## Compute Quadratic variations of m(t) and M(t) and 
    ## covariance of (m(t), M(t))							
    vm = sigsq
    vqm = kappa-vm**2  
    vmqm = gamma
    
    ## Define rho^2(t) 															
    termr = 1-((vmqm**2)/(vm*vqm))								
    rho = 1/termr
    if(is.na(rho) | is.infinite(rho))
      rho <- 1
    
    ## Define eta 								
    eta = vmqm/(vm*vqm)	
    if(is.na(eta) | is.infinite((eta)))
      eta = 0
    
    ## Define vectors astr and bstr						
    astr = rho*(-dermu/vm + dersigsq*eta)							
    bstr = rho*(dermu*eta - dersigsq/vqm)
    
    ## Define Derivatives of m(t) and M(t) 							
    derm = -dermu
    derqm = 2*m*derm - dersigsq	
    
    ## Derivatives of variations <m>(t), <M>(t) and <m,M>(t)							
    dervm = dersigsq						
    dervqm = derkappa - 2*sigsq*dersigsq
    dervmqm = dergamma
    
    #Note: derrho=0 - and it is!
    ru=vm*vqm
    rv=vm*vqm-vmqm**2
    rdu=vm*dervqm+vqm*dervm
    rdv=vm*dervqm+vqm*dervm -2*vmqm*dervmqm
    derrho=(rv*rdu-ru*rdv)/(rv**2)
    if(sum(is.na(derrho)) >0)
      derrho <- matrix(rep(0),pdim,1)
    ## Derivative of eta(t)							
    num1=vm*vqm*dervmqm
    num2=vmqm*(vm*dervqm + vqm*dervm)
    den=(vm**2)*(vqm**2)
    dereta=(num1-num2)/den
    if(sum(is.infinite(dereta)) >0 | sum(is.na(dereta))>0)
      dereta <- matrix(rep(0),pdim,1)
    ## Derivatives of astr and bstr							
    ## For astr						
    terma1 = der2mu/vm - dermu%*%t(dervm)/(vm**2)
    terma2 = der2sigsq*eta + dersigsq%*%t(dereta)
    terma3 = -dermu/vm + dersigsq*eta							
    derastr = -rho*terma1 + rho*terma2 + derrho%*%t(terma3)
    #derastr = -rho*terma1 + rho*terma2 
    
    ## For bstr						
    termb1 = der2mu*eta + dermu%*%t(dereta)													
    termb21 = der2sigsq/vqm 
    termb22 = (dersigsq%*%t(dervqm))/(vqm**2)
    if(sum(is.infinite(termb22)>0 | sum(is.na(termb22))>0))
      termb22 <- matrix(0, ncol=pdim, nrow = pdim)
    termb2 = termb21 - termb22
    termb3 = dermu*eta - dersigsq/vqm
    derbstr = rho*termb1 - rho*termb2 + derrho%*%t(termb3)
    
    
    ### Recursive Formulas						
    ## Compute Kinv(t): Information 							
    termk1[,,t] = astr%*%t(derm) + m * derastr 							
    termk2[,,t] = bstr%*%t(derqm) + qm * derbstr 
    
    ## add derivative of penalty here
    penalty <- p.lam.prime(thehat[,,t-1], thehat[,,1] ,lambda, a, tau)
    pp <- rbind(pp, penalty)
    #penalty <- p.lam.prime(thehat[,,t-1],lambda, a)
    penalty <- matrix(penalty, ncol = 1, nrow = pdim)
    
    ## Need Numerical FIX when termk1 and termk2 are NA
    if ( sum(as.numeric(is.na(termk1[,,t]))) > 0){
      foundFlag1 = 0
      for (w in 1:(t-1)) {
        if ( (foundFlag1 == 0) & sum(as.numeric(is.na(termk1[,,t-w]))) == 0 ) {
          termk1[,,t] = termk1[,,(t-w)]
          foundFlag1 = 1
        }
        else
          termk1[,,t]= matrix(rep(.Machine$double.xmax/3,pdim*pdim),pdim,pdim)
      }
    }
    
    if ( sum(as.numeric(is.na(termk2[,,t]))) > 0){
      foundFlag2 = 0
      for (w in 1:(t-1)) {
        if ( (foundFlag2 == 0) & sum(as.numeric(is.na(termk2[,,t-w]))) == 0 ) {
          termk2[,,t] = termk2[,,(t-w)]
          foundFlag2 = 1
        }
        else
          termk2[,,t]= matrix(rep(.Machine$double.xmax/3,pdim*pdim),pdim,pdim)
      }
    }
    
    
    p_lambda <- p.lam.2prime(thehat[,,1], lambda, a, tau)
    p_lambda <- diag(p_lambda)
    if(t >= 0*n){
      kinv[,,t] = kinv[,,t-1] - (termk1[,,t]+termk2[,,t] - p_lambda)
      
    }else{
      kinv[,,t] = kinv[,,t-1] - (termk1[,,t]+termk2[,,t])
      
    }
    
    
    s <- svd(kinv[,,t])
    
    inf.index <- which(is.infinite(s$d))
    D <- diag(1/s$d)
    for(i in 1:length(s$d)){
      if(is.infinite(D[i,i])) 
        D[i,i] <- .Machine$double.xmin*100
    }
    kmat[,,t]<- s$v%*%D%*%t(s$u)
    
    termt = astr*m + bstr*qm
    if(sum(as.numeric(is.na(termt))) >0 ) {
      break
    }
    thehat[,,t] = thehat[,,t-1] + kmat[,,t]%*%(termt - penalty)
    
  }
  
  
  
  
  finalest <- return( thehat[,1,])
  
  
}



