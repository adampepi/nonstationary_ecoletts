##Code for  state space models and simulations for Pepi, Holyoak & Karban 2021, Ecology Letters
rm(list = ls())
library(rjags)
library(R2jags)
library(AICcmodavg)
library(lattice)
library(MCMCvis)
library(tidyverse)
library(tidybayes)
library(ggplot2)
library(ggstance)
library(bayestestR)
setwd("~/Documents/Research/dissertation/time series analyses/nonstationary_modelling")
cats1<-read.csv('bodega.cats2.csv')
str(cats1)
cats<-cats1
cats$Precip<-as.numeric(scale(cats$precip))
cats$cat.count<-as.integer(cats$lupine.count)
#First part
cats2<-cats1[1:19,]
cats2$Precip<-as.numeric(scale(cats2$precip))
cats2$cat.count<-as.integer(cats2$lupine.count)
str(cats2)
#Second part
cats3<-cats1[19:34,]
cats3$Precip<-as.numeric(scale(cats3$precip))
cats3$cat.count<-as.integer(cats3$lupine.count)
str(cats)

### Delayed DD Gompertz precip

#Full series

sink("tigermodelprecipgompdelayed.jags")  # name of the txt file
cat("
    model{
  #Priors
    mean.r[1] ~ dnorm(0, 0.1)  # Prior for mean growth rate
    mean.r[2] ~ dnorm(0, 0.1)
    sigma.proc ~ dunif(0, 10)             # Prior for sd of state process
    sigma2.proc <- pow(sigma.proc, 2)
    tau.proc <- pow(sigma.proc, -2)
          beta_precip ~ dnorm(0,0.01)
    beta_dd1 ~ dnorm(0,0.01)
    beta_dd2 ~ dnorm(0,0.01)
    beta_0 ~ dnorm(0, 0.01)
      #State process
    for (t in 3:T){
    mean.r[t] <- beta_0 + beta_dd1 * logN.est[t-1] + beta_dd2 * logN.est[t-2] + beta_precip * precip[t-1] 
    }
      #Error model
    for(t in 1:T){    
    logN.est[t] ~ dnorm(mean.r[t], tau.proc)   
    N.est[t] <- exp(logN.est[t])
  # Observation process
    y[t] ~  dpois(N.est[t]*logarea[t])
    }
             }
    
    ", fill=T)
sink()

jags.data15<- list(y = as.integer(cats$cat.count), T = length(cats$Year),precip=cats$Precip,logarea=cats$lupine.area)
year<-cats$Year
inits15 <- function(){list(sigma.proc =0.9, beta_0= 3,beta_dd1= -0.1, beta_dd2= -0.35, beta_precip=0.03)}
parameters15<- c("mean.r", "sigma2.proc", "N.est","logN.est","beta_0","beta_dd1",
                 'beta_dd2',"beta_precip")
ni <- 20000
nt <- 1
nb <- 1000
nc <- 3
pgompdelayed <- jags(jags.data15, inits15, parameters15, "tigermodelprecipgompdelayed.jags", n.chains = nc,
                     n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())
pgd<-jags.model("tigermodelprecipgompdelayed.jags",data=jags.data15, inits=inits15,n.chains = nc)

summary(pgompdelayed)
print(pgompdelayed, digits = 3)

theta.samples1 <- coda.samples(pgd, n.iter=10000, thin=10, c("beta_precip","beta_dd1","beta_dd2","beta_0",'sigma2.proc'))

par(mar=rep(2,4)); plot(theta.samples1)

#First half


sink("tigermodelprecipgompdelayed2.jags")  # name of the txt file
cat("
    model{
  #Priors
    mean.r[1] ~ dnorm(0, 0.1)  # Prior for mean growth rate
    mean.r[2] ~ dnorm(0, 0.1)
    sigma.proc ~ dunif(0, 10)             # Prior for sd of state process
    sigma2.proc <- pow(sigma.proc, 2)
    tau.proc <- pow(sigma.proc, -2)
          beta_precip ~ dnorm(0,0.01)
    beta_dd1 ~ dnorm(0,0.01)
    beta_dd2 ~ dnorm(0,0.01)
    beta_0 ~ dnorm(0, 0.01)
      #State process
    for (t in 3:T){
    mean.r[t] <- beta_0 + beta_dd1 * logN.est[t-1] + beta_dd2 * logN.est[t-2] + beta_precip * precip[t-1] 
    }
      #Error model
    for(t in 1:T){    
    logN.est[t] ~ dnorm(mean.r[t], tau.proc)   
    N.est[t] <- exp(logN.est[t])
  # Observation process
    y[t] ~  dpois(N.est[t]*logarea[t])
    }
             }
    
    ", fill=T)
sink()

jags.data15<- list(y = as.integer(cats2$cat.count), T = length(cats2$Year),precip=cats2$Precip,logarea=cats2$lupine.area)
year<-cats$Year
inits15 <- function(){list(sigma.proc =0.9, beta_0= 3,beta_dd1= -0.1, beta_dd2= -0.35, beta_precip=0.03)}
parameters15<- c("mean.r", "sigma2.proc", "N.est","logN.est","beta_0","beta_dd1",
                 'beta_dd2',"beta_precip")
ni <- 20000
nt <- 1
nb <- 1000
nc <- 3
pgompdelayed2 <- jags(jags.data15, inits15, parameters15, "tigermodelprecipgompdelayed2.jags", n.chains = nc,
                     n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())
pgd2<-jags.model("tigermodelprecipgompdelayed2.jags",data=jags.data15, inits=inits15,n.chains = nc)

summary(pgompdelayed)
print(pgompdelayed, digits = 3)


theta.samples2 <- coda.samples(pgd2, n.iter=10000, thin=10, c("beta_precip","beta_dd1","beta_dd2","beta_0","sigma2.proc"))

par(mar=rep(2,4)); plot(theta.samples2)

#Second half

sink("tigermodelprecipgompdelayed3.jags")  # name of the txt file
cat("
    model{
  #Priors
    mean.r[1] ~ dnorm(0, 0.1)  # Prior for mean growth rate
    mean.r[2] ~ dnorm(0, 0.1)
    sigma.proc ~ dunif(0, 10)             # Prior for sd of state process
    sigma2.proc <- pow(sigma.proc, 2)
    tau.proc <- pow(sigma.proc, -2)
          beta_precip ~ dnorm(0,0.01)
    beta_dd1 ~ dnorm(0,0.01)
    beta_dd2 ~ dnorm(0,0.01)
    beta_0 ~ dnorm(0, 0.01)
      #State process
    for (t in 3:T){
    mean.r[t] <- beta_0 + beta_dd1 * logN.est[t-1] + beta_dd2 * logN.est[t-2] + beta_precip * precip[t-1] 
    }
      #Error model
    for(t in 1:T){    
    logN.est[t] ~ dnorm(mean.r[t], tau.proc)   
    N.est[t] <- exp(logN.est[t])
  # Observation process
    y[t] ~  dpois(N.est[t]*logarea[t])
    }
             }
    
    ", fill=T)
sink()

jags.data15<- list(y = as.integer(cats3$cat.count), T = length(cats3$Year),precip=cats3$Precip,logarea=cats3$lupine.area)
year<-cats$Year
inits15 <- function(){list(sigma.proc =0.9, beta_0= 3,beta_dd1= -0.1, beta_dd2= -0.35, beta_precip=0.03)}
parameters15<- c("mean.r", "sigma2.proc", "N.est","logN.est","beta_0","beta_dd1",
                 'beta_dd2',"beta_precip")
ni <- 20000
nt <- 1
nb <- 1000
nc <- 3
pgompdelayed3 <- jags(jags.data15, inits15, parameters15, "tigermodelprecipgompdelayed3.jags", n.chains = nc,
                     n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())
pgd3<-jags.model("tigermodelprecipgompdelayed3.jags",data=jags.data15, inits=inits15,n.chains = nc)

summary(pgompdelayed)
print(pgompdelayed, digits = 3)
plotModelOutput(pgompdelayed,log((cats$lupine.count+1)/cats$lupine.area))

theta.samples3 <- coda.samples(pgd3, n.iter=10000, thin=10, c("beta_precip","beta_dd1","beta_dd2","beta_0",'sigma2.proc'))

par(mar=rep(2,4)); plot(theta.samples3)

coda::traceplot(theta.samples3)

png("pgomp12.png", width=600, height=400)
MCMCplot(object = pgompdelayed, 
         object2 = pgompdelayed2, 
         col='black',
         col2='blue',
         params = c("beta_0",'beta_dd1','beta_dd2',"beta_precip"),
         labels = c('Intercept', 'Direct density-dependence', 'Delayed density-dependence', 'Precipitation'),
         xlim=c(-3,3))
legend('topright',inset=.05,col = c('black',"blue","red"),lty=1,lwd=3,legend=c('Whole series',"Before threshold","After threshold"))
dev.off()
png("pgomp3.png", width=600, height=400, bg = "transparent")
MCMCplot(object = pgompdelayed3, 
         col='red',
         params = c("beta_0",'beta_dd1','beta_dd2',"beta_precip"),
         labels = c('Intercept', 'Direct density-dependence', 'Delayed density-dependence', 'Precipitation'),
         xlim=c(-3,3))


dev.off()


m1<-tidy_draws(theta.samples1)
m1
spread_draws()
point_interval(m1)

III<-theta.samples1 %>%
  gather_draws(beta_0,beta_dd1,beta_dd2,beta_precip)%>%
  mutate(Model = "Whole series")
I<-theta.samples2 %>%
  gather_draws(beta_0,beta_dd1,beta_dd2,beta_precip)%>%
  mutate(Model = "Part I")
II<-theta.samples3 %>%
  gather_draws(beta_0,beta_dd1,beta_dd2,beta_precip)%>%
 mutate(Model = "Part II")

overall<-bind_rows(III,II,I)
overall$Model<-as.factor(overall$Model)

levels(overall$Model)
overall$Model = factor(overall$Model,levels(overall$Model)[c(2,1,3)])
levels(overall$Model)

overall%>%
ggplot(aes(y = .variable , x = .value,color=Model)) +
  stat_halfeyeh(position = position_dodgev(height = .9),.width = c(.90, .95))+geom_vline(xintercept=0,lty=2)+xlab(label="Value")+
  ylab(label="Parameter")+
  scale_y_discrete(labels = c(expression(alpha["0"]),expression(alpha["1"]),expression(alpha["2"]),expression(beta["Precip"])))+
  theme_classic()+theme(text=element_text(size=20))+guides(color = guide_legend(reverse = TRUE))+ scale_color_manual(values=c("red", "blue", "black"))

ci(overall.1$beta_dd1[overall.1$Model=="Whole series"],method="HDI",ci=0.9)
ci(overall.1$beta_dd1[overall.1$Model=="Part I"],method="HDI",ci=0.9)
ci(overall.1$beta_dd1[overall.1$Model=="Part II"],method="HDI",ci=0.9)
ci(overall.1$beta_dd2[overall.1$Model=="Whole series"],method="HDI",ci=0.9)
ci(overall.1$beta_dd2[overall.1$Model=="Part I"],method="HDI",ci=0.9)
ci(overall.1$beta_dd2[overall.1$Model=="Part II"],method="HDI",ci=0.9)
ci(overall.1$beta_precip[overall.1$Model=="Whole series"],method="HDI",ci=0.9)
ci(overall.1$beta_precip[overall.1$Model=="Part I"],method="HDI",ci=0.9)
ci(overall.1$beta_precip[overall.1$Model=="Part II"],method="HDI",ci=0.9)




III.1<-theta.samples1 %>%
  spread_draws(beta_0,beta_dd1,beta_dd2,beta_precip)%>%
  mutate(Model = "Whole series")
I.1<-theta.samples2 %>%
  spread_draws(beta_0,beta_dd1,beta_dd2,beta_precip)%>%
  mutate(Model = "Part I")
II.1<-theta.samples3 %>%
  spread_draws(beta_0,beta_dd1,beta_dd2,beta_precip)%>%
  mutate(Model = "Part II")

overall.1<-bind_rows(III.1,II.1,I.1)
overall.1$Model<-as.factor(overall.1$Model)

overall.1

overall.1$beta1<-overall.1$beta_dd1
overall.1$beta2<-overall.1$beta_dd2

sample(overall.1$beta_0[overall.1$Model=="Whole series"],size=1)

xes<-seq(-2,2,by=.1)
croyama<-function(x)-0.25*x^2

curve1<-croyama(xes)
curve<-data.frame(beta1=xes,beta2=curve1)

xx<-c(-2,0,2)
yy<-c(-1,1,-1)
beta1<-c(0.085,-0.619,0.45)
beta2<-c(-0.325,-0.258,-0.239)
model<-c("Whole series","Part I","Part II")

royama<-data.frame(beta1=beta1,beta2=beta2,Model=model)

triangle<-data.frame(beta1=xx,beta2=yy)

vline<-data.frame(x1=0,x2=0,y1=-1,y2=1)

labels<-data.frame(quadrant=c("I","II","III","IV","I'","II'","III'","IV'"),y=c(0.35,0.35,-.6,-.6,0.35,0.35,-1.15,-1.15),x=c(0.3,-0.3,-.7,.7,-1.5,1.5,-.7,.7))
str(overall.1)
detach("package:biwavelet", unload=TRUE)
ggplot()+geom_point(data=overall.1,mapping=aes(x=beta1,y=beta2, color=Model),size=1,alpha = 0.25)+xlim(-2,2)+ylim(-1.2,1)+theme_classic()+geom_point(data=royama,mapping=aes(x=beta1,y=beta2, fill=Model),size=4,shape=21)+xlab(expression(alpha["1"]))+ylab(expression(alpha["2"]))+geom_line(data=curve,mapping=aes(x=beta1,y=beta2))+
  geom_polygon(data=triangle,mapping=aes(x=beta1,y=beta2),colour="black",fill=NA)+geom_segment(data=vline,mapping=aes(x=x1,xend=x2,y=y1,yend=y2))+
  geom_text(data=labels,mapping=aes(x=x,y=y,label=quadrant),size=7)+
geom_segment(aes(x =royama$beta1[2],y = royama$beta2[2],xend = royama$beta1[3],yend = royama$beta2[3]),arrow=arrow(),data=royama)+theme(text=element_text(size=20))


##### Simlations -- using posterior
library(biwavelet)
cats1$logmean<-log(cats1$cat.mean)
cats1$logntm1<-log(cats1$ntm1)
cats1$logntm2<-log(cats1$ntm1)
cats1$precipscaled<-scale(cats1$precip)

cats2<-cats1[2:34,]


cats3<-cats1[3:34,]
oseries<-cbind(cats1$Year,cats1$logmean)

nrands<-1000
wtc.0 = wt(oseries)

##Generalised simulation function


##values drawn from posterior

simulation<-function(a0samp=overall.1$beta_0[overall.1$Model=="Whole series"],
                     a1samp=overall.1$beta_dd1[overall.1$Model=="Whole series"],
                     a2samp=overall.1$beta_dd2[overall.1$Model=="Whole series"],
                     bprecipsamp=overall.1$beta_precip[overall.1$Model=="Whole series"],
                     reps=1000
  
)
{
##fixed values
log1<-function (N,Ntm1,a0,a1,a2,b1,p)  a0 + a1*N+a2*Ntm1+b1*p
tf<-34 #run time
n0<-0.293  #pop init size
n1<--2.957
n2<--1.822
samples<-vector(length=reps)
for(i in 1:reps){

a0<-sample(a0samp,size=1)
a1<-sample(a1samp,size=1)
a2<-sample(a2samp,size=1)
b1<-sample(bprecipsamp,size=1)

n<-rep(NA,tf)  #make vector
n[1] = n0 #put init pop
n[2] = n1
n[3]=n2
precip<-cats1$precipscaled

for(t in 3:(tf-1)){  #t-1 to match lengths
  n[t+1]<-log1(N=n[t],Ntm1=n[t-1],a0=a0,a1=a1,a2=a2,b1=b1,p=precip[t])
}


sim1<-cbind(cats1$Year,n)
wtc.1 = wt(sim1)
samples[i]<-wdist(wtc.0$wave,wtc.1$wave)
}
print(samples)
}


s1<-simulation(a0samp=overall.1$beta_0[overall.1$Model=="Whole series"],
               a1samp=overall.1$beta_dd1[overall.1$Model=="Whole series"],
               a2samp=overall.1$beta_dd2[overall.1$Model=="Whole series"],
               bprecipsamp=overall.1$beta_precip[overall.1$Model=="Whole series"],
               reps=10000)

hist(s1)
s1.1<-as.data.frame(s1)
s1.1$sim<-"Whole Series"
str(s1.1)
ggplot(s1.1,aes(s1))+geom_density()

###First Part
s2<-simulation(a0samp=overall.1$beta_0[overall.1$Model=="Part I"],
               a1samp=overall.1$beta_dd1[overall.1$Model=="Part I"],
               a2samp=overall.1$beta_dd2[overall.1$Model=="Part I"],
               bprecipsamp=overall.1$beta_precip[overall.1$Model=="Part I"],
               reps=10000)

hist(s2)
s2.1<-as.data.frame(s2)
s2.1$sim<-"Part I"
s2.1$s1<-s2.1$s2
ggplot(s2.1,aes(s2))+geom_density()


###Second Part
s3<-simulation(a0samp=overall.1$beta_0[overall.1$Model=="Part II"],
               a1samp=overall.1$beta_dd1[overall.1$Model=="Part II"],
               a2samp=overall.1$beta_dd2[overall.1$Model=="Part II"],
               bprecipsamp=overall.1$beta_precip[overall.1$Model=="Part II"],
               reps=10000)

hist(s3)
s3.1<-as.data.frame(s3)
s3.1$sim<-"Part II"
s3.1$s1<-s3.1$s3
ggplot(s3.1,aes(s3))+geom_density()

###Whole series --- precip from part II

s4<-simulation(a0samp=overall.1$beta_0[overall.1$Model=="Whole series"],
               a1samp=overall.1$beta_dd1[overall.1$Model=="Whole series"],
               a2samp=overall.1$beta_dd2[overall.1$Model=="Whole series"],
               bprecipsamp=overall.1$beta_precip[overall.1$Model=="Part II"],
               reps=10000)

hist(s4)
s4.1<-as.data.frame(s4)
s4.1$sim<-"Whole Series Precip Same"
s4.1$s1<-s4.1$s4
ggplot(s1.1,aes(s4))+geom_density()

###First Part
s5<-simulation(a0samp=overall.1$beta_0[overall.1$Model=="Part I"],
               a1samp=overall.1$beta_dd1[overall.1$Model=="Part I"],
               a2samp=overall.1$beta_dd2[overall.1$Model=="Part I"],
               bprecipsamp=overall.1$beta_precip[overall.1$Model=="Part II"],
               reps=10000)

hist(s5)
s5.1<-as.data.frame(s5)
s5.1$sim<-"Part I Precip Same"
s5.1$s1<-s5.1$s5
ggplot(s5.1,aes(s5))+geom_density()

##Part II no DD

s6<-simulation(a0samp=overall.1$beta_0[overall.1$Model=="Part II"],
               a1samp=0,
               a2samp=0,
               bprecipsamp=overall.1$beta_precip[overall.1$Model=="Part II"],
               reps=10000)

hist(s6)
s6.1<-as.data.frame(s6)
s6.1$sim<-"Part II No DD"
s6.1$s1<-s6.1$s6
ggplot(s6.1,aes(s6))+geom_density()

s1.1

str(s1.1)
str(s2.1)
str(s3.1)
str(s4.1)
str(s5.1)
str(s6.1)

values<-c(s1.1$s1,s2.1$s1,s3.1$s1,s4.1$s1,s5.1$s1,s6.1$s1)
length(values)
simulations<-c(s1.1$sim,s2.1$sim,s3.1$sim,s4.1$sim,s5.1$sim,s6.1$sim)
length(simulations)
allsims<-data.frame(values=values,simulation=simulations)
str(allsims)
allsims$simulation<-as.factor(allsims$simulation)
str(allsims)

ggplot(allsims,aes(values,lty=simulation))+geom_density()



hdi1<-ci(allsims$values[allsims$simulation=='Whole Series'], method="HDI",ci=0.95)
hdi2<-ci(allsims$values[allsims$simulation=='Whole Series Precip Same'], method="HDI",ci=0.95)
hdi3<-ci(allsims$values[allsims$simulation=='Part II No DD'], method="HDI",ci=0.95)
hdi4<-ci(allsims$values[allsims$simulation=='Part II'], method="HDI",ci=0.95)
hdi5<-ci(allsims$values[allsims$simulation=='Part I Precip Same'], method="HDI",ci=0.95)
hdi6<-ci(allsims$values[allsims$simulation=='Part I'], method="HDI",ci=0.95)
hdi1
hdi2
hdi3
hdi4
hdi5
hdi6
hdi1<-ci(allsims$values[allsims$simulation=='Whole Series'], method="HDI",ci=0.9)
hdi2<-ci(allsims$values[allsims$simulation=='Whole Series Precip Same'], method="HDI",ci=0.9)
hdi3<-ci(allsims$values[allsims$simulation=='Part II No DD'], method="HDI",ci=0.9)
hdi4<-ci(allsims$values[allsims$simulation=='Part II'], method="HDI",ci=0.9)
hdi5<-ci(allsims$values[allsims$simulation=='Part I Precip Same'], method="HDI",ci=0.9)
hdi6<-ci(allsims$values[allsims$simulation=='Part I'], method="HDI",ci=0.9)
hdi1
hdi2
hdi3
hdi4
hdi5
hdi6
str(allsims)
ggplot(data=allsims,aes(y = simulation , x = values)) +
  stat_halfeyeh(position = position_dodgev(height = .9),.width = c(.90, .95))+xlab(label="Value")+
  ylab(label="Scenario")+
  theme_classic()+theme(text=element_text(size=20))

#####Posterior predictive checks
cats1<-read.csv('bodega.cats2.csv')
str(cats1)
cats<-cats1
cats$Precip<-as.numeric(scale(cats$precip))
cats$cat.count<-as.integer(cats$lupine.count)

cats2<-cats1[1:19,]
cats2$Precip<-as.numeric(scale(cats2$precip))
cats2$cat.count<-as.integer(cats2$lupine.count)
str(cats2)

cats3<-cats1[19:34,]
cats3$Precip<-as.numeric(scale(cats3$precip))
cats3$cat.count<-as.integer(cats3$lupine.count)
str(cats)

sink("tigermodelprecipgompdelayedppcheck.jags")  # name of the txt file
cat("
    model{
  #Priors
    mean.r[1] ~ dnorm(0, 0.1)  # Prior for mean growth rate
    mean.r[2] ~ dnorm(0, 0.1)
    sigma.proc ~ dunif(0, 10)             # Prior for sd of state process
    sigma2.proc <- pow(sigma.proc, 2)
    tau.proc <- pow(sigma.proc, -2)
          beta_precip ~ dnorm(0,0.01)
    beta_dd1 ~ dnorm(0,0.01)
    beta_dd2 ~ dnorm(0,0.01)
    beta_0 ~ dnorm(0, 0.01)
      #State process
    for (t in 3:T){
    mean.r[t] <- beta_0 + beta_dd1 * logN.est[t-1] + beta_dd2 * logN.est[t-2] + beta_precip * precip[t-1] 
    }
      #Error model
    for(t in 1:T){    
    logN.est[t] ~ dnorm(mean.r[t], tau.proc)   
    N.est[t] <- exp(logN.est[t])
    
  # Observation process
    y[t] ~  dpois(N.est[t]*logarea[t])
    y.new[t] ~  dpois(N.est[t]*logarea[t])
    res[t]<-y[t]-N.est[t]
    res.new[t]<-y.new[t]-N.est[t]
    }
    #Derived parameters
  fit <- sum(res[])
  fit.new <- sum(res.new[])
             }
    
    ", fill=T)
sink()

jags.data15<- list(y = as.integer(cats$cat.count), T = length(cats$Year),precip=cats$Precip,logarea=cats$lupine.area)
year<-cats$Year
inits15 <- function(){list(sigma.proc =0.9, beta_0= 3,beta_dd1= -0.1, beta_dd2= -0.35, beta_precip=0.03)}
parameters15<- c("mean.r", "sigma2.proc", "N.est","logN.est","beta_0","beta_dd1",
                 'beta_dd2',"beta_precip","fit",'fit.new','y.new')
ni <- 20000
nt <- 1
nb <- 1000
nc <- 3
pgompdelayed <- jags(jags.data15, inits15, parameters15, "tigermodelprecipgompdelayedppcheck.jags", n.chains = nc,
                     n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())
pgd<-jags.model("tigermodelprecipgompdelayed.jags",data=jags.data15, inits=inits15,n.chains = nc)

summary(pgompdelayed)
print(pgompdelayed, digits = 3)

theta.samples1 <- coda.samples(pgd, n.iter=10000, thin=10, c("beta_precip","beta_dd1","beta_dd2","beta_0"))

par(mar=rep(1,4)); plot(theta.samples1,oma=c(2,2,2,2))

library(jagsUI)
pgompdelayed <- jags(jags.data15, inits15, parameters15, "tigermodelprecipgompdelayedppcheck.jags", n.chains = nc,
                     n.thin = nt, n.iter = ni, n.burnin = nb)
pp.check(x=pgompdelayed, observed = 'fit', simulated = 'fit.new')

crosscorr.plot(theta.samples1)
crosscorr(theta.samples1)
#pp.check(x=pgompdelayed, observed = 'y', simulated = 'y.new')

str(pgompdelayed)
ynew<-pgompdelayed$sims.list$y.new
year<-cats$Year
y = as.integer(cats$cat.count)
y
plot(x=year,y=y+1,type='l',log='y',xlab="Year",ylab="Count")
lines(x=year,y=ynew[sample(nrow(ynew),size=1,replace=FALSE),]+1, col="blue")
lines(x=year,y=ynew[sample(nrow(ynew),size=1,replace=FALSE),]+1, col="blue")
lines(x=year,y=ynew[sample(nrow(ynew),size=1,replace=FALSE),]+1, col="blue")
lines(x=year,y=ynew[sample(nrow(ynew),size=1,replace=FALSE),]+1, col="blue")
lines(x=year,y=ynew[sample(nrow(ynew),size=1,replace=FALSE),]+1, col="blue")
lines(x=year,y=ynew[sample(nrow(ynew),size=1,replace=FALSE),]+1, col="blue")
lines(x=year,y=ynew[sample(nrow(ynew),size=1,replace=FALSE),]+1, col="blue")
lines(x=year,y=ynew[sample(nrow(ynew),size=1,replace=FALSE),]+1, col="blue")
lines(x=year,y=ynew[sample(nrow(ynew),size=1,replace=FALSE),]+1, col="blue")
lines(x=year,y=ynew[sample(nrow(ynew),size=1,replace=FALSE),]+1, col="blue")
lines(x=year,y=ynew[sample(nrow(ynew),size=1,replace=FALSE),]+1, col="blue")
lines(x=year,y=ynew[sample(nrow(ynew),size=1,replace=FALSE),]+1, col="blue")
lines(x=year,y=ynew[sample(nrow(ynew),size=1,replace=FALSE),]+1, col="blue")
lines(x=year,y=ynew[sample(nrow(ynew),size=1,replace=FALSE),]+1, col="blue")
lines(x=year,y=ynew[sample(nrow(ynew),size=1,replace=FALSE),]+1, col="blue")
lines(x=year,y=y+1, lwd=3)


#First half


sink("tigermodelprecipgompdelayed2ppcheck.jags")  # name of the txt file
cat("
    model{
  #Priors
    mean.r[1] ~ dnorm(0, 0.1)  # Prior for mean growth rate
    mean.r[2] ~ dnorm(0, 0.1)
    sigma.proc ~ dunif(0, 10)             # Prior for sd of state process
    sigma2.proc <- pow(sigma.proc, 2)
    tau.proc <- pow(sigma.proc, -2)
          beta_precip ~ dnorm(0,0.01)
    beta_dd1 ~ dnorm(0,0.01)
    beta_dd2 ~ dnorm(0,0.01)
    beta_0 ~ dnorm(0, 0.01)
      #State process
    for (t in 3:T){
    mean.r[t] <- beta_0 + beta_dd1 * logN.est[t-1] + beta_dd2 * logN.est[t-2] + beta_precip * precip[t-1] 
    }
      #Error model
    for(t in 1:T){    
    logN.est[t] ~ dnorm(mean.r[t], tau.proc)   
    N.est[t] <- exp(logN.est[t])
  # Observation process
    y[t] ~  dpois(N.est[t]*logarea[t])
    y.new[t] ~  dpois(N.est[t]*logarea[t])
    res[t]<-y[t]-N.est[t]
    res.new[t]<-y.new[t]-N.est[t]
    }
    #Derived parameters
  fit <- sum(res[])
  fit.new <- sum(res.new[])
  }
    
    ", fill=T)
sink()

jags.data15<- list(y = as.integer(cats2$cat.count), T = length(cats2$Year),precip=cats2$Precip,logarea=cats2$lupine.area)
year<-cats$Year
inits15 <- function(){list(sigma.proc =0.9, beta_0= 3,beta_dd1= -0.1, beta_dd2= -0.35, beta_precip=0.03)}
parameters15<- c("mean.r", "sigma2.proc", "N.est","logN.est","beta_0","beta_dd1",
                 'beta_dd2',"beta_precip","fit",'fit.new','y.new')
ni <- 20000
nt <- 1
nb <- 1000
nc <- 3
pgompdelayed2 <- jags(jags.data15, inits15, parameters15, "tigermodelprecipgompdelayed2ppcheck.jags", n.chains = nc,
                      n.thin = nt, n.iter = ni, n.burnin = nb)
pgd2<-jags.model("tigermodelprecipgompdelayed2.jags",data=jags.data15, inits=inits15,n.chains = nc)

summary(pgompdelayed)
print(pgompdelayed, digits = 3)


theta.samples2 <- coda.samples(pgd2, n.iter=10000, thin=10, c("beta_precip","beta_dd1","beta_dd2","beta_0"))

par(mar=rep(1,4)); plot(theta.samples2)


pp.check(x=pgompdelayed2, observed = 'fit', simulated = 'fit.new')
crosscorr(theta.samples2)

ynew<-pgompdelayed2$sims.list$y.new
ynew

y = as.integer(cats2$cat.count)
y
year<-cats2$Year
plot(x=cats2$Year,y=y+1,type='l',log='y',xlab="Year",ylab="Count")
lines(x=year,y=ynew[sample(nrow(ynew),size=1,replace=FALSE),]+1, col="blue")
lines(x=year,y=ynew[sample(nrow(ynew),size=1,replace=FALSE),]+1, col="blue")
lines(x=year,y=ynew[sample(nrow(ynew),size=1,replace=FALSE),]+1, col="blue")
lines(x=year,y=ynew[sample(nrow(ynew),size=1,replace=FALSE),]+1, col="blue")
lines(x=year,y=ynew[sample(nrow(ynew),size=1,replace=FALSE),]+1, col="blue")
lines(x=year,y=ynew[sample(nrow(ynew),size=1,replace=FALSE),]+1, col="blue")
lines(x=year,y=ynew[sample(nrow(ynew),size=1,replace=FALSE),]+1, col="blue")
lines(x=year,y=ynew[sample(nrow(ynew),size=1,replace=FALSE),]+1, col="blue")
lines(x=year,y=ynew[sample(nrow(ynew),size=1,replace=FALSE),]+1, col="blue")
lines(x=year,y=ynew[sample(nrow(ynew),size=1,replace=FALSE),]+1, col="blue")
lines(x=year,y=ynew[sample(nrow(ynew),size=1,replace=FALSE),]+1, col="blue")
lines(x=year,y=ynew[sample(nrow(ynew),size=1,replace=FALSE),]+1, col="blue")
lines(x=year,y=ynew[sample(nrow(ynew),size=1,replace=FALSE),]+1, col="blue")
lines(x=year,y=ynew[sample(nrow(ynew),size=1,replace=FALSE),]+1, col="blue")
lines(x=year,y=ynew[sample(nrow(ynew),size=1,replace=FALSE),]+1, col="blue")
lines(x=year,y=y+1, lwd=3)


#Second half

sink("tigermodelprecipgompdelayed3ppcheck.jags")  # name of the txt file
cat("
    model{
  #Priors
    mean.r[1] ~ dnorm(0, 0.1)  # Prior for mean growth rate
    mean.r[2] ~ dnorm(0, 0.1)
    sigma.proc ~ dunif(0, 10)             # Prior for sd of state process
    sigma2.proc <- pow(sigma.proc, 2)
    tau.proc <- pow(sigma.proc, -2)
          beta_precip ~ dnorm(0,0.01)
    beta_dd1 ~ dnorm(0,0.01)
    beta_dd2 ~ dnorm(0,0.01)
    beta_0 ~ dnorm(0, 0.01)
      #State process
    for (t in 3:T){
    mean.r[t] <- beta_0 + beta_dd1 * logN.est[t-1] + beta_dd2 * logN.est[t-2] + beta_precip * precip[t-1] 
    }
      #Error model
    for(t in 1:T){    
    logN.est[t] ~ dnorm(mean.r[t], tau.proc)   
    N.est[t] <- exp(logN.est[t])
  # Observation process
    y[t] ~  dpois(N.est[t]*logarea[t])
    y.new[t] ~  dpois(N.est[t]*logarea[t])
    res[t]<-y[t]-N.est[t]
    res.new[t]<-y.new[t]-N.est[t]
    }
    #Derived parameters
  fit <- sum(res[])
  fit.new <- sum(res.new[])
    
             }
    
    ", fill=T)
sink()

jags.data15<- list(y = as.integer(cats3$cat.count), T = length(cats3$Year),precip=cats3$Precip,logarea=cats3$lupine.area)
year<-cats$Year
inits15 <- function(){list(sigma.proc =0.9, beta_0= 3,beta_dd1= -0.1, beta_dd2= -0.35, beta_precip=0.03)}
parameters15<- c("mean.r", "sigma2.proc", "N.est","logN.est","beta_0","beta_dd1",
                 'beta_dd2',"beta_precip","fit",'fit.new','y.new')
ni <- 20000
nt <- 1
nb <- 1000
nc <- 3
pgompdelayed3 <- jags(jags.data15, inits15, parameters15, "tigermodelprecipgompdelayed3ppcheck.jags", n.chains = nc,
                      n.thin = nt, n.iter = ni, n.burnin = nb )
pgd3<-jags.model("tigermodelprecipgompdelayed3.jags",data=jags.data15, inits=inits15,n.chains = nc)

summary(pgompdelayed3)
print(pgompdelayed3, digits = 3)

theta.samples3<- coda.samples(pgd3, n.iter=10000, thin=10, c("beta_precip","beta_dd1","beta_dd2","beta_0"))

par(mar=rep(1,4)); plot(theta.samples3)

crosscorr(theta.samples3)

pp.check(x=pgompdelayed3, observed = 'fit', simulated = 'fit.new')


ynew<-pgompdelayed3$sims.list$y.new
ynew

y = as.integer(cats3$cat.count)
y
year<-cats3$Year
dev.off()
plot(x=year,y=y+1,type='l',log='y',xlab="Year",ylab="Count")
lines(x=year,y=ynew[sample(nrow(ynew),size=1,replace=FALSE),]+1, col="blue")
lines(x=year,y=ynew[sample(nrow(ynew),size=1,replace=FALSE),]+1, col="blue")
lines(x=year,y=ynew[sample(nrow(ynew),size=1,replace=FALSE),]+1, col="blue")
lines(x=year,y=ynew[sample(nrow(ynew),size=1,replace=FALSE),]+1, col="blue")
lines(x=year,y=ynew[sample(nrow(ynew),size=1,replace=FALSE),]+1, col="blue")
lines(x=year,y=ynew[sample(nrow(ynew),size=1,replace=FALSE),]+1, col="blue")
lines(x=year,y=ynew[sample(nrow(ynew),size=1,replace=FALSE),]+1, col="blue")
lines(x=year,y=ynew[sample(nrow(ynew),size=1,replace=FALSE),]+1, col="blue")
lines(x=year,y=ynew[sample(nrow(ynew),size=1,replace=FALSE),]+1, col="blue")
lines(x=year,y=ynew[sample(nrow(ynew),size=1,replace=FALSE),]+1, col="blue")
lines(x=year,y=ynew[sample(nrow(ynew),size=1,replace=FALSE),]+1, col="blue")
lines(x=year,y=ynew[sample(nrow(ynew),size=1,replace=FALSE),]+1, col="blue")
lines(x=year,y=ynew[sample(nrow(ynew),size=1,replace=FALSE),]+1, col="blue")
lines(x=year,y=ynew[sample(nrow(ynew),size=1,replace=FALSE),]+1, col="blue")
lines(x=year,y=ynew[sample(nrow(ynew),size=1,replace=FALSE),]+1, col="blue")
lines(x=year,y=y+1, lwd=3)

##Simulation plots
samples<-matrix(0,ncol=34-3,nrow=100)

samples

simulation2<-function(a0samp=overall.1$beta_0[overall.1$Model=="Whole series"],
                     a1samp=overall.1$beta_dd1[overall.1$Model=="Whole series"],
                     a2samp=overall.1$beta_dd2[overall.1$Model=="Whole series"],
                     bprecipsamp=overall.1$beta_precip[overall.1$Model=="Whole series"],
                     reps=100
                     
)
{
  ##fixed values
  log1<-function (N,Ntm1,a0,a1,a2,b1,p)  a0 + a1*N+a2*Ntm1+b1*p
  tf<-34 #run time
  n0<-0.293  #pop init size
  n1<--2.957
  n2<--1.822
  samples<-matrix(0,ncol=tf,nrow=reps)
  for(i in 1:reps){
    
    a0<-sample(a0samp,size=1)
    a1<-sample(a1samp,size=1)
    a2<-sample(a2samp,size=1)
    b1<-sample(bprecipsamp,size=1)
    
    n<-rep(NA,tf)  #make vector
    n[1] = n0 #put init pop
    n[2] = n1
    n[3]=n2
    precip<-cats1$precipscaled
    
    for(t in 3:(tf-1)){  #t-1 to match lengths
      n[t+1]<-log1(N=n[t],Ntm1=n[t-1],a0=a0,a1=a1,a2=a2,b1=b1,p=precip[t])
    }
    
    
    sim1<-cbind(cats1$Year,n)
    samples[i,]<-n
  }
  print(samples)
}



##Have to have JagsUI unloaded?
par(mfrow=c(2,3))


s1<-simulation2(a0samp=overall.1$beta_0[overall.1$Model=="Whole series"],
               a1samp=overall.1$beta_dd1[overall.1$Model=="Whole series"],
               a2samp=overall.1$beta_dd2[overall.1$Model=="Whole series"],
               bprecipsamp=overall.1$beta_precip[overall.1$Model=="Whole series"],
               reps=100)

s1.1<-exp(s1)

s1.1[1,]
y = cats$lupine.count/cats$lupine.area
year<-cats$Year
plot(x=year,y=y+1,type='l',log='y',xlab="Year",ylab="Density",ylim=c(1,20),main="Whole series")
for(i in 1:50){
  lines(x=year,y=s1.1[sample(nrow(s1.1),size=1,replace=FALSE),]+1, col=rgb(0,0,1,0.1))
}

lines(x=year,y=y+1, lwd=3)



###First Part
s2<-simulation2(a0samp=overall.1$beta_0[overall.1$Model=="Part I"],
               a1samp=overall.1$beta_dd1[overall.1$Model=="Part I"],
               a2samp=overall.1$beta_dd2[overall.1$Model=="Part I"],
               bprecipsamp=overall.1$beta_precip[overall.1$Model=="Part I"],
               reps=100)

s1.1<-exp(s2)
y = cats$lupine.count/cats$lupine.area
year<-cats$Year
plot(x=year,y=y+1,type='l',log='y',xlab="Year",ylab="Density",ylim=c(1,20),main='Part I')
for(i in 1:50){
  lines(x=year,y=s1.1[sample(nrow(s1.1),size=1,replace=FALSE),]+1, col=rgb(0,0,1,0.1))
}
lines(x=year,y=y+1, lwd=3)

###Second Part
s3<-simulation2(a0samp=overall.1$beta_0[overall.1$Model=="Part II"],
               a1samp=overall.1$beta_dd1[overall.1$Model=="Part II"],
               a2samp=overall.1$beta_dd2[overall.1$Model=="Part II"],
               bprecipsamp=overall.1$beta_precip[overall.1$Model=="Part II"],
               reps=100)

s1.1<-exp(s3)
y = cats$lupine.count/cats$lupine.area
year<-cats$Year
plot(x=year,y=y+1,type='l',log='y',xlab="Year",ylab="Density",ylim=c(1,20),main='Part II')
for(i in 1:50){
  lines(x=year,y=s1.1[sample(nrow(s1.1),size=1,replace=FALSE),]+1, col=rgb(0,0,1,0.1))
}
lines(x=year,y=y+1, lwd=3)
###Whole series --- precip from part II

s4<-simulation2(a0samp=overall.1$beta_0[overall.1$Model=="Whole series"],
               a1samp=overall.1$beta_dd1[overall.1$Model=="Whole series"],
               a2samp=overall.1$beta_dd2[overall.1$Model=="Whole series"],
               bprecipsamp=overall.1$beta_precip[overall.1$Model=="Part II"],
               reps=100)

s1.1<-exp(s4)
y = cats$lupine.count/cats$lupine.area
year<-cats$Year
plot(x=year,y=y+1,type='l',log='y',xlab="Year",ylab="Density",ylim=c(1,20),main='Whole Series Precip Same')
for(i in 1:50){
  lines(x=year,y=s1.1[sample(nrow(s1.1),size=1,replace=FALSE),]+1, col=rgb(0,0,1,0.1))
}
lines(x=year,y=y+1, lwd=3)

###First Part
s5<-simulation2(a0samp=overall.1$beta_0[overall.1$Model=="Part I"],
               a1samp=overall.1$beta_dd1[overall.1$Model=="Part I"],
               a2samp=overall.1$beta_dd2[overall.1$Model=="Part I"],
               bprecipsamp=overall.1$beta_precip[overall.1$Model=="Part II"],
               reps=100)

s1.1<-exp(s5)
y = cats$lupine.count/cats$lupine.area
year<-cats$Year
plot(x=year,y=y+1,type='l',log='y',xlab="Year",ylab="Density",ylim=c(1,20),main="Part I Precip Same")
for(i in 1:30){
  lines(x=year,y=s1.1[sample(nrow(s1.1),size=1,replace=FALSE),]+1, col=rgb(0,0,1,0.1))
}
lines(x=year,y=y+1, lwd=3)

##Part II no DD

s6<-simulation2(a0samp=overall.1$beta_0[overall.1$Model=="Part II"],
               a1samp=0,
               a2samp=0,
               bprecipsamp=overall.1$beta_precip[overall.1$Model=="Part II"],
               reps=100)

s1.1<-exp(s6)
y = cats$lupine.count/cats$lupine.area
year<-cats$Year
plot(x=year,y=y+1,type='l',log='y',xlab="Year",ylab="Density",ylim=c(1,20),main='Part II No DD')
for(i in 1:30){
  lines(x=year,y=s1.1[sample(nrow(s1.1),size=1,replace=FALSE),]+1, col=rgb(0,0,1,0.1))
}
lines(x=year,y=y+1, lwd=3)

