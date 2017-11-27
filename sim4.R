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
para <- list(seq1=c(0.5,-0.3),
             seq2=c(0.3,0.3),
             seq3=c(0.5,-0.2),
             seq4=c(0.2,0.4))
cpts <- c(0.4, 0.6, 0.8)   
ps <- c(1,1,1,1)
qs <- c(0,0,0,0)
feps <- rep("gamma", 4)
par1 <- c(3,2,3,2)
par2 <- c(4,5,4,5)


#### Data Generation & Parameter Estiamtion using PEF ######
set.seed(12345)
PEF <- matrix(0, nrow=1)
B <- 1
for(j in 1:B){
  print(j)
  sim <- mix.logacd(8000, 500, para, cpts, ps, qs, feps, par1, par2)
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
  x <- res2[1,-c(1:1500)]
  p <- find_peaks(x) +1500
  v <- find_peaks(-x) + 1500
  
  mat <- matrix(sort(c(p,v)), nrow = 1)
  PEF <- rbind.fill.matrix(PEF, mat)
}


res2[1,] <- res2[1,]+difference
aa <- res2
aa[,c(1:1500)] <- NA
dat3 <- data.frame(y=as.numeric(t(aa)),
                   pars = rep(1:6, each=8000),
                   x= rep(1:8000, 6),
                   method = rep("PEF", 8000*6))
dat3$pars <- factor(dat3$pars)
labels <- list(expression(omega), expression(alpha[1]),
               expression(alpha[2]), expression(alpha[3]),
               expression(alpha[4]), expression(alpha[5]))


##### Simulated Data Plot ####
fig3 <- ggplot(data = dat3, aes(y=y, x=x, linetype=pars))+
  geom_line( )+
  geom_vline(xintercept=c(3200,4800,6400), colour="black",
             linetype="longdash", size=1)+
  theme_bw()+
  ylab("Estimated Values")+xlab("(b)")+
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=15,face="bold"),
        strip.text.x = element_text(size = 18),
        strip.text.y = element_text(size = 18),
        legend.text = element_text(size = 18),
        legend.position = "right",
        panel.grid.major.x=element_blank())+
  guides(linetype = guide_legend( title = " "))+
  scale_linetype_manual(breaks = levels(dat3$pars),
                        values=1:6,labels =labels)

##### Trace Plot  ####
dat4 <- data.frame(y=sim$y,
                   x= c(1:8000),
                   method = c("Duration", 8000))

fig4 <- ggplot(data=dat4, aes(y=y,x=x))+
  geom_line(colour="#999999")+
  geom_vline(xintercept=c(3200,4800,6400), colour="black",
             linetype="longdash", size=1)+
  theme_bw()+
  ylab("Durations")+xlab("(a)")+
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=15,face="bold"),
        strip.text.x = element_text(size = 18),
        strip.text.y = element_text(size = 18),
        legend.text = element_text(size = 18),
        legend.position = "right",
        panel.grid.major.x=element_blank())

##### Multiple Plot #####
multiplot(fig4,fig3,cols=2)
