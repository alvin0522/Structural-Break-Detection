######
## Simulation of Piecewise LogACD(p,q) process
######
# Allows different order of p and q
# Allows different error distributions
# Allows (p,0) model. 
# Modified and tested at 08/05/2017

mix.logacd <- function (n,nb,para,cpts, ps, qs, feps,par1, par2) 
{
  
  ########################################################################
  ## Check for stationarity constraint
  ########################################################################
  n.seg=length(para)
  
  #len.para=rep(0,n.seg)
  #for (i in 1:n.seg){
  #  len.para[i]=length(para[[i]])-1
  #}
  
  # nn is a vector containing the change-location of each segment
  nn=rep(0,n.seg)
  for (j in 1:(n.seg-1)){nn[j]=n*cpts[j]}
  nn[n.seg]=n
  
  #if (sum(alpha)+sum(beta)>=1) 
  #{ cat("Error: The simulated process is not weakly stationary")
  #} else {
  
  ## Start simulating Durations
  nt <- n + nb  # Total number of simulated datapoints
  omega <- NULL
  for(i in 1:n.seg){
    omega <- c(omega, para[[i]][1])
  }
  
  alpha <- matrix(0,nrow = 1)
  for(i in 1:n.seg){
    mat <- matrix(para[[i]][2:(1+ps[i])], nrow = 1)
    alpha <- rbind.fill.matrix(alpha, mat)
  }
  
  beta <- matrix(0,nrow = 1)
  for(i in 1:n.seg){
    if(qs[i] == 0){
      mat <- matrix(NA, nrow = 1)
    }else{
      mat <- matrix(para[[i]][(2+ps[i]):(1+ps[i]+qs[i])], nrow = 1)
    }
    beta <- rbind.fill.matrix(beta, mat)
  }
  
  mxpq <- max(ps[1],qs[1])		
  
  
  ########################################################################
  ## Initialize psis, xs, x, and psi vectors
  ########################################################################
  psis <- rep(1,nt)
  logxs <- rep(1,nt)
  xs <- rep(1,nt)
  
  lagged.logxs <- rep(1,ps[1])
  lagged.psis <- rep(1,qs[1])
  
  x <-rep(NA,n)	# we save last n entries of xs into x 
  psi <-rep(NA,n)	# we save last n entries of psis into psi
  eps <- rep(1, nt)
  
  ########################################################################
  ## Generate errors eps from a positive-valued distribution
  ## Choices are: exponential, Weibull, or Gamma
  ########################################################################
  nn.seg <- nb + nn
  if (feps[1]=="exponential") 
  {
    ## Generate errors from exponential distribution
    eps[1:nn.seg[1]] <- rexp(nn.seg[1],par1[1])
  } else if (feps[1]=="weibull") {		
    ## Generate errors from weibull distribution			
    eps[1:nn.seg[1]] <- rweibull(nn.seg[1], shape=par1[1], scale=par2[1]) 	 
  } else if (feps[1]=="gamma") {
    ## Generate errors from gamma distribution			
    # randomly generate errors from gamma distribution			
    eps[1:nn.seg[1]] <-  rgamma(nn.seg[1],shape=par1[1],scale=par2[1])
  }
  for(j in 2:n.seg){
    indx <- c((nn.seg[j-1]+1):nn.seg[j])
    if (feps[j]=="exponential") 
    {
      ## Generate errors from exponential distribution
      eps[(nn.seg[j-1]+1):nn.seg[j]] <-  rexp(length(indx),par1[j])
    } else if (feps[j]=="weibull") {		
      ## Generate errors from weibull distribution			
      eps[(nn.seg[j-1]+1):nn.seg[j]] <-  rweibull(length(indx), shape=par1[j], scale=par2[j]) 	 
    } else if (feps[j]=="gamma") {
      ## Generate errors from gamma distribution			
      # randomly generate errors from gamma distribution			
      eps[(nn.seg[j-1]+1):nn.seg[j]] <- rgamma(length(indx),shape=par1[j],scale=par2[j])
    }
  }
  
  
  ########################################################################
  ## Compute psis and xs
  ########################################################################
  
  for (t in 1:mxpq){
    psis[t] = omega[1]
    logxs[t]= psis[t] + log(eps[t]) -log(mean(eps))	
    xs[t]=exp(logxs[t])
  }
  
  for (t in (mxpq+1):nt) 
  {
    lagged.logxs <- c(logxs[(t-1):(t-ps[1])])
    lagged.psis <- c(psis[(t-1):(t-qs[1])])
    if(is.na(beta[2,])){
      psis[t] = omega[1] + alpha[2,(1:ps[1])]%*%lagged.logxs
    }else{
      psis[t] = omega[1] + alpha[2,(1:ps[1])]%*%lagged.logxs + beta[2,(1:qs[1])]%*%lagged.psis 
    }
    logxs[t] = psis[t] + log(eps[t]) 	-log(mean(eps))			
    xs[t] = exp(logxs[t] )
    for(j in 2:n.seg){
      lagged.logxs <- c(logxs[(t-1):(t-ps[j])])
      lagged.psis <- c(psis[(t-1):(t-qs[j])])
      if(is.na(beta[j+1,])){
        psis_star = omega[j] + alpha[j+1,(1:ps[j])]%*%lagged.logxs
      }else{
        psis_star = omega[j] + alpha[j+1,(1:ps[j])]%*%lagged.logxs + beta[j+1,(1:qs[1])]%*%lagged.psis 
      }
      logxs[t] = logxs[t] + (psis_star - psis[t])*((t - nb)>nn[j-1])
      if((t - nb)>nn[j-1]){
        psis[t] <- psis_star
      }
      xs[t] = exp(logxs[t])
    }
  }
  
  #}
  ## End Durations Simulation of psis and xs
  
  ## Save last n values as time series into x and psi (discard first nb as burn-in)		
  x <- ts(xs[(nb+1):nt])	 	
  psi <- ts(psis[(nb+1):nt])	
  
  ## output as dataframe, Col 1 has x and Col 2 has psi		
  return(list(y=x, para=para, cp=nn, 
              error.info=list(dist=feps, par1=par1, par2=par2)))
}