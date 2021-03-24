
rm(list = ls())
library(rjags)
library(R2jags)
library(AICcmodavg)
library(lattice)

cats1<-read.csv('bodega.cats3.csv')
str(cats1)
cats<-cats1
cats$Precip<-as.numeric(scale(cats$precip))
cats$cat.count<-as.integer(cats$lupine.count)

###Model 1: Gompertz DD
sink("tigermodelgomp.jags")  # name of the txt file
cat("
    model{
      # Priors and constraints
       mean.r[1] ~ dnorm(0, 0.1)  # Prior for mean growth rate
       sigma.proc ~ dunif(0, 10)             # Prior for sd of state process
       sigma2.proc <- pow(sigma.proc, 2)
       tau.proc <- pow(sigma.proc, -2)
          beta_dd1 ~ dnorm(0,0.01)
       beta_0 ~ dnorm(0, 0.01)
      #State process
    for (t in 2:T){
    mean.r[t] <- beta_0 + beta_dd1 * logN.est[t-1] 
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

jags.data1<- list(y = as.integer(cats$cat.count), T = length(cats$Year),logarea=cats$lupine.area)
year<-cats$Year
inits1 <- function(){list(sigma.proc = runif(1, 0, 1),
                          beta_0= runif(1, 0, 1),beta_dd1= runif(1, 0, 1)
)}
parameters1 <- c("mean.r", "sigma2.proc", "N.est","logN.est","beta_0","beta_dd1")
ni <- 2000
nt <- 1
nb <- 100
nc <- 3
gomp<- jags(jags.data1, inits1, parameters1, "tigermodelgomp.jags", n.chains = nc,
            n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd()) 
g<-jags.model("tigermodelgomp.jags",data=jags.data1, inits=inits1,n.chains = nc)

summary(gomp)
print(gomp, digits = 3)
Y<-log((cats$lupine.count+1)/cats$lupine.area)
plotModelOutput <- function(jagsmodel, Y) {
  attach.jags(jagsmodel)
  x <- seq(1,length(Y))
  XPred <- cbind(apply(logN.est,2,quantile,0.025), apply(logN.est,2,mean), apply(logN.est,2,quantile,0.975))
  ylims <- c(min(c(Y,XPred), na.rm=TRUE), max(c(Y,XPred), na.rm=TRUE))
  plot(Y, col="white",ylim=ylims, xlab="t",ylab="ln(Caterpillars on lupine)")
  polygon(c(x,rev(x)), c(XPred[,1], rev(XPred[,3])), col="grey70",border=NA)
  lines(XPred[,2])
  lines(Y,lty=2)
}

plotModelOutput(gomp,log((cats$lupine.count+1)/cats$lupine.area))

theta.samples1 <- coda.samples(g, n.iter=10000, thin=10, c("beta_0","beta_dd1"))

par(mar=rep(1,4)); plot(theta.samples1)

### Model 2 Gompertz and precip 
sink("tigermodelgompprecip.jags")  # name of the txt file
cat("
    model{
       # Priors and constraints
        mean.r[1] ~ dnorm(0, 0.1)  # Prior for mean growth rate
       sigma.proc ~ dunif(0, 10)             # Prior for sd of state process
    sigma2.proc <- pow(sigma.proc, 2)
    tau.proc <- pow(sigma.proc, -2)
        beta_dd1 ~ dnorm(0,0.01)
    beta_precip ~ dnorm(0,0.01)
    beta_0 ~ dnorm(0, 0.01)
      #State process
    for (t in 2:T){
    mean.r[t] <- beta_0 + beta_dd1 * logN.est[t-1] + beta_precip * precip[t-1] 
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

jags.data2<- list(y = as.integer(cats$cat.count), T = length(cats$Year),precip=cats$Precip,logarea=cats$lupine.area)
year<-cats$Year
inits2<- function(){list(sigma.proc = runif(1, 0, 1),beta_0= runif(1, 0, 1),
                         beta_dd1= runif(1, 0, 1),beta_precip= runif(1, 0, 1))}
parameters2<- c("mean.r", "sigma2.proc", "N.est","logN.est","beta_0","beta_dd1",
                "beta_precip")
ni <- 20000
nt <- 1
nb <- 1000
nc <- 3
gompprecip <- jags(jags.data2, inits2, parameters2, "tigermodelgompprecip.jags", n.chains = nc,
                   n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())
gp<-jags.model("tigermodelgompprecip.jags",data=jags.data2, inits=inits2,n.chains = nc)
summary(gompprecip)
print(gompprecip, digits = 3)

plotModelOutput(gompprecip,log((cats$lupine.count+1)/cats$lupine.area))

theta.samples2 <- coda.samples(gp, n.iter=10000, thin=10, c("beta_0","beta_dd1","beta_precip"))

par(mar=rep(1,4)); plot(theta.samples2)

###Model 3: Precip
sink("tigermodelprecip.jags")  # name of the txt file
cat("
    model{
       # Priors and constraints
        mean.r[1] ~ dnorm(0, 0.1)  # Prior for mean growth rate
       sigma.proc ~ dunif(0, 10)             # Prior for sd of state process
    sigma2.proc <- pow(sigma.proc, 2)
    tau.proc <- pow(sigma.proc, -2)
        beta_precip ~ dnorm(0,0.01)
    beta_0 ~ dnorm(0, 0.01)
    
  #State process
    for (t in 2:T){
    mean.r[t] <- beta_0 + logN.est[t-1] + beta_precip * precip[t-1] 
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

jags.data3<- list(y = as.integer(cats$cat.count), T = length(cats$Year), 
                  precip=cats$Precip,logarea=cats$lupine.area)
year<-cats$Year
inits3 <- function(){list(sigma.proc = runif(1, 0, 1),
                          beta_0= 5,
                          beta_precip= -0.014)}
parameters3 <- c("mean.r", "sigma2.proc", "N.est","logN.est","beta_0",
                 'beta_precip')
ni <- 20000
nt <- 1
nb <- 1000
nc <- 3
precip <- jags(jags.data3, inits3, parameters3, "tigermodelprecip.jags", n.chains = nc,
               n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())
p<-jags.model("tigermodelprecip.jags",data=jags.data3, inits=inits3,n.chains = nc)

summary(precip)
print(precip, digits = 3)

plotModelOutput(precip,log((cats$lupine.count+1)/cats$lupine.area))

theta.samples3 <- coda.samples(p, n.iter=10000, thin=10, c("beta_0","beta_precip"))

par(mar=rep(1,4)); plot(theta.samples3)

###Model 4: Delayed DD Gompertz
sink("tigermodelgompdelayed.jags")  # name of the txt file
cat("
    model{
    # Priors and constraints
      mean.r[1] ~ dnorm(0, 0.1)  # Prior for mean growth rate
    mean.r[2] ~ dnorm(0, 0.1)
    sigma.proc ~ dunif(0, 10)             # Prior for sd of state process
    sigma2.proc <- pow(sigma.proc, 2)
    tau.proc <- pow(sigma.proc, -2)
    beta_dd1 ~ dnorm(0,0.01)
    beta_dd2 ~ dnorm(0,0.01)
    beta_0 ~ dnorm(0, 0.01)
  #State process
    for (t in 3:T){
    mean.r[t] <- beta_0 + beta_dd1 * logN.est[t-1] + beta_dd2 * logN.est[t-2] 
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
jags.data4<- list(y = as.integer(cats$cat.count), T = length(cats$Year),logarea=cats$lupine.area)
year<-cats$Year
inits4 <- function(){list(sigma.proc = runif(1, 0, 1),
                          beta_0= runif(1, 0, 1),beta_dd1= runif(1, 0, 1),beta_dd2= runif(1, 0, 1))}
parameters4 <- c("mean.r", "sigma2.proc", "N.est","logN.est","beta_0","beta_dd1",
                 'beta_dd2')
ni <- 20000
nt <- 1
nb <- 1000
nc <- 3
gompdelayed <- jags(jags.data4, inits4, parameters4, "tigermodelgompdelayed.jags", n.chains = nc,
                    n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())
gd<-jags.model("tigermodelgompdelayed.jags",data=jags.data4, inits=inits4,n.chains = nc)
summary(gompdelayed)
print(gompdelayed, digits = 3)

plotModelOutput(gompdelayed,log((cats$lupine.count+1)/cats$lupine.area))

theta.samples4 <- coda.samples(gd, n.iter=10000, thin=10, c("beta_0","beta_dd1","beta_dd2"))

par(mar=rep(1,4)); plot(theta.samples5)

###Model 5: Intercept
sink("tigermodelint.jags")  # name of the txt file
cat("
    model{
        # Priors and constraints
    mean.r[1] ~ dnorm(0, 0.1)  # Prior for mean growth rate
       sigma.proc ~ dunif(0, 10)             # Prior for sd of state process
    sigma2.proc <- pow(sigma.proc, 2)
    tau.proc <- pow(sigma.proc, -2)
      alpha ~ dgamma(1, 1)
    beta_0 ~ dnorm(0, 0.01)
    #State process
    for (t in 2:T){
    mean.r[t] <- beta_0 +  logN.est[t-1] 
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
jags.data5<- list(y = as.integer(cats$cat.count), T = length(cats$Year),logarea=cats$lupine.area)
year<-cats$Year
inits5 <- function(){list(sigma.proc = runif(1, 0, 1),
                          beta_0= runif(1, 0, 1))}
parameters5 <- c("mean.r", "sigma2.proc", "N.est","logN.est","beta_0")
ni <- 20000
nt <- 1
nb <- 1000
nc <- 3
intercept <- jags(jags.data5, inits5, parameters5, "tigermodelint.jags", n.chains = nc,
                  n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())
i<-jags.model("tigermodelint.jags",data=jags.data5, inits=inits5,n.chains = nc)

summary(intercept)
print(intercept, digits = 3)
plotModelOutput(intercept,log((cats$lupine.count+1)/cats$lupine.area))

theta.samples5 <- coda.samples(i, n.iter=10000, thin=10, c("beta_0"))

par(mar=rep(1,4)); plot(theta.samples5)

###Model 6: Delayed DD Gompertz precip

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

jags.data6<- list(y = as.integer(cats$cat.count), T = length(cats$Year),precip=cats$Precip,logarea=cats$lupine.area)
year<-cats$Year
inits6 <- function(){list(sigma.proc =0.9, beta_0= 3,beta_dd1= -0.1, beta_dd2= -0.35, beta_precip=0.03)}
parameters6<- c("mean.r", "sigma2.proc", "N.est","logN.est","beta_0","beta_dd1",
                'beta_dd2',"beta_precip")
ni <- 20000
nt <- 1
nb <- 1000
nc <- 3
pgompdelayed <- jags(jags.data6, inits6, parameters6, "tigermodelprecipgompdelayed.jags", n.chains = nc,
                     n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())
pgd<-jags.model("tigermodelprecipgompdelayed.jags",data=jags.data6, inits=inits6,n.chains = nc)

summary(pgompdelayed)
print(pgompdelayed, digits = 3)
plotModelOutput(pgompdelayed,log((cats$lupine.count+1)/cats$lupine.area))

theta.samples6 <- coda.samples(pgd, n.iter=10000, thin=10, c("beta_precip","beta_dd1","beta_dd2","beta_0"))

par(mar=rep(1,4)); plot(theta.samples6)
######Model selection



s <- lapply(jsamples1, unclass)
sapply(s, sum)

sum(jsamples1$deviance)
gpdic<-jags.samples(gp,c("deviance", "WAIC"), type="mean",n.iter=100000)
pdic<-jags.samples(p,c("deviance", "WAIC"), type="mean",n.iter=100000)
gdic<-jags.samples(g,c("deviance", "WAIC"), type="mean",n.iter=100000)
idic<-jags.samples(i,c("deviance", "WAIC"), type="mean",n.iter=100000)
gddic<-jags.samples(gd,c("deviance", "WAIC"), type="mean",n.iter=100000)
pgddic<-jags.samples(pgd,c("deviance", "WAIC"), type="mean",n.iter=100000)



gpdic1<-sum(gpdic$deviance +  gpdic$WAIC)
gdic1<-sum(gdic$deviance +  gdic$WAIC)
pdic1<-sum(pdic$deviance +  pdic$WAIC)
idic1<-sum(idic$deviance +  idic$WAIC)
gddic1<-sum(gddic$deviance +  gddic$WAIC)
pgddic1<-sum(pgddic$deviance +  pgddic$WAIC)


DICs<-c(gpdic1,
        gdic1,
        pdic1,
        idic1,
        gddic1,
        pgddic1)
names(DICs)<-c('gpdic1',
               'gdic1',
               'pdic1',
               'idic1',
               'gddic1',
               'pgddic1')
DICs

DICsorted<-sort(DICs)
DICsorted
dDIC<-DICsorted-DICsorted[1]
dDIC
