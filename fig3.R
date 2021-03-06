rm(list=ls())
library("ggplot2")
library("plyr")
setwd("path where all codes are stored")
source("./multiplot.R")
source("./core_funs.R")
source("./user_define.R")
source("./PEF.R")
source("./mix_logacd_sim.R")

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
res1 <- recursive_mat(x = dat, initval = est, p = p, q = q,
                      momente = momente)
res1[1,] <- res1[1,] + difference

##### Black-White Plot  ####
aa <- res2
aa[,c(1:1500)] <- NA
dat3 <- data.frame(y=as.numeric(t(aa)),
                   pars = rep(1:6, each=7500),
                   x= rep(1:7500, 6),
                   method = rep("PEF", 7500*6))

labels <- list(expression(omega), expression(alpha[1]),
               expression(alpha[2]), expression(alpha[3]),
               expression(alpha[4]), expression(alpha[5]))


dat3$pars <- factor(dat3$pars)
fig <- ggplot(data = dat3, aes(y=y, x=x, linetype=pars))+
  geom_line()+
  geom_vline(xintercept=c(3000,5250), colour="black",
             linetype="dotted", size=1.5)+
  theme_bw()+
  ylab("Estimated Values")+xlab("Index")+
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
