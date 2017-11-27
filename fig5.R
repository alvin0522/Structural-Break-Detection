rm(list=ls())
library("ggplot2")
library("plyr")
setwd("path where all codes are stored")
source("./multiplot.R")
source("./core_funs.R")
source("./user_define.R")
source("./PEF.R")
source("./mix_logacd_sim.R")
source("./multiplot.R")
source("./find_peaks.R")

##### Simulation Set Up ####
set.seed(12345)
para <- list(seq1=c(0.1,0.3),
             seq2=c(0.1,-0.2),
             seq3=c(0.1,0.3))
cpts <- c(0.4, 0.7)    #0.5, 0.7; 0.3, 0.8
ps <- c(1,1,1)
qs <- c(0,0,0)
feps <- rep("weibull", 3)
par1 <- rep(3, 3)
par2 <- rep(4, 3)

sim <- mix.logacd(7500, 500, para, cpts, ps, qs, feps, par1, par2)


dat <- sim$y
dat1 <- data.frame(y=dat, x=c(1:length(dat)))


##### trace-plot  ####
dat <- sim$y
model.order = c(5,0)
pen.par =  c(2,3.7)
p <- model.order[1]
q <- model.order[2]
pdim <- p+q+1
lambda <- pen.par[1]
a <- pen.par[2]
est <- finitval(dat,p,q)
if (q==0) est <- est[-(pdim+1)]
psih <- psi_cal(dat,est,p,q)
eps <- dat/exp(psih)
ed <- eps_dist2(eps)
momente <- c(ed$mue, ed$vare, ed$skewe, ed$kurte)
difference <- log(momente[1]) - mean(log(eps))
res2 <- recursive_mat2(x = dat, initval = est, p = p, q = q,
                       momente = momente, lambda = lambda, a = a)
res2[1,] <- res2[1,] + difference


##### Apply FindPeaks Procedure ####
span <- 0.25
x <- res2[2, -c(1:1500)]
n <- length(x)
y <- x
mu.y.loc <- y
mu.y.loc <- (x[1:(n-2)] + x[2:(n-1)] + x[3:n])/3
mu.y.loc <- c(mu.y.loc[1], mu.y.loc, mu.y.loc[n-2])

y.loess <- loess(x~I(1:n), span=span)
yy <- y.loess[[2]]
y <- y.loess[[2]]
sig.y <- var(y.loess$resid, na.rm=TRUE)^0.5
DX.1 <- sign(diff(mu.y.loc, na.pad = FALSE))
pks <- which(diff(DX.1, na.pad = FALSE) < 0 & DX.1[-(n-1)] > 0) + 1
cpts <- NULL
cpts.thresh <- NULL
lspan=0.05
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
points(x=out, y=yy[out], col="blue")

dat1 <- data.frame(y= x, x=1:length(x))
dat2 <- data.frame(y= yy, x=1:length(x))

## step 1: c=0.1 ##

dat3 <- dat2
dat3$y[-out] <- NA


fig1 <- ggplot(dat1, aes(y=y,x=x))+
  geom_line()+
  geom_line(data=dat2, aes(y=y,x=x),linetype="longdash", size=1)+
  geom_point(data=dat3,aes(y=y,x=x),colour="#999999")+
  theme_bw()+
  ylab("Estimated Values")+xlab("(a)")+
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=15,face="bold"),
        strip.text.x = element_text(size = 18),
        strip.text.y = element_text(size = 18),
        legend.text = element_text(size = 18),
        legend.position = "right",
        panel.grid.major.x=element_blank())

## step 2: c=0.5 ##
### find var out ###
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
thresh <- 0.5
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
dat3 <- dat2
dat3$y[-out] <- NA

fig2 <- ggplot(dat1, aes(y=y,x=x))+
  geom_line()+
  geom_line(data=dat2, aes(y=y,x=x),linetype="longdash", size=1)+
  geom_point(data=dat3,aes(y=y,x=x),colour="#999999")+
  theme_bw()+
  ylab("Estimated Values")+xlab("(b)")+
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=15,face="bold"),
        strip.text.x = element_text(size = 18),
        strip.text.y = element_text(size = 18),
        legend.text = element_text(size = 18),
        legend.position = "right",
        panel.grid.major.x=element_blank())

## step 3: c=1 ##
### find var out ###
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
thresh <- 1
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
dat3 <- dat2
dat3$y[-out] <- NA

fig3 <- ggplot(dat1, aes(y=y,x=x))+
  geom_line()+
  geom_line(data=dat2, aes(y=y,x=x),linetype="longdash", size=1)+
  geom_point(data=dat3,aes(y=y,x=x),colour="#999999")+
  theme_bw()+
  ylab("Estimated Values")+xlab("(c)")+
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=15,face="bold"),
        strip.text.x = element_text(size = 18),
        strip.text.y = element_text(size = 18),
        legend.text = element_text(size = 18),
        legend.position = "right",
        panel.grid.major.x=element_blank())

## last step : c=1.5 ##
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
thresh <- 1.8
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
dat3 <- dat2
dat3$y[-out] <- NA

fig4 <- ggplot(dat1, aes(y=y,x=x))+
  geom_line()+
  geom_line(data=dat2, aes(y=y,x=x),linetype="longdash", size=1)+
  geom_point(data=dat3,aes(y=y,x=x),colour="#999999")+
  theme_bw()+
  ylab("Estimated Values")+xlab("(d)")+
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=15,face="bold"),
        strip.text.x = element_text(size = 18),
        strip.text.y = element_text(size = 18),
        legend.text = element_text(size = 18),
        legend.position = "right",
        panel.grid.major.x=element_blank())

#### Multiplot ####
multiplot(fig1, fig3, fig2,  fig4, cols=2)