rm(list=ls())
library(biwavelet)

setwd("~/aapepi@ucdavis.edu/dissertation/time series analyses")
cats1<-read.csv('bodega.cats2.csv')
cats2<-read.csv('threeseries2.csv')
str(cats1)

str(cats2)
lupinecats<-log((cats1$lupine.count+1)/cats1$lupine.area)
meancats<-log(cats1$cat.mean)
hemcats<-log(cats1$hemlock.count)
hist(lupinecats)
hist(cats1$precip)
t1<-cbind(cats1$Year,scale(lupinecats))
t2<-cbind(cats1$Year,scale(cats1$precip))
hist(t1)
hist(t2)

nrands = 10000
par(mfrow=c(1,1))
wtc.AB = wtc(t1, t2, nrands = nrands)
plot(wtc.AB, plot.phase = T, lty.coi = 1, col.coi = "grey", lwd.coi = 2, 
     lwd.sig = 2, arrow.lwd = 0.03, arrow.len = 0.12, ylab = "Period", xlab = "Year", 
     plot.cb = F, main = "Wavelet Coherence: Cats vs Precip")


par(mfrow=c(4,1),mar = c(2,4,1.5,4))

plot(x=cats1$Year,y=as.integer(cats1$cat.mean*100),type='b',xlab=NA,ylab="Caterpillars",log = "y",ylim=c(1,500),main="Caterpillars and Precipitation")

par(new = T)
with(cats1, plot(x=cats1$Year, y=cats1$precip, axes=F, xlab=NA, ylab=NA,type='b',col='blue',lty=2))
axis(side = 4)
mtext(side = 4, line = 3, 'Precipitation, cm',cex=0.7)


wt1<-wt(t1)
plot(wt1, type = "power.corr.norm", main = "Caterpillars",ylab="Period",xlab=NA)
wt2<-wt(t2)
plot(wt2, type = "power.corr.norm", main = "Precipitation",ylab="Period",xlab=NA)


plot(wtc.AB, plot.phase = T, lty.coi = 1, col.coi = "grey", lwd.coi = 2, 
     lwd.sig = 2, arrow.lwd = 0.03, arrow.len = 0.03, ylab = "Period", xlab = "Year", 
     plot.cb = F, main = "Wavelet Coherence")




#Global wavelets

wt1$gws <- apply(wt1$power, 1, mean)

wt1$gws.sig <- wt.sig(d=cbind(cats1$Year,scale(lupinecats)), dt=wt1$dt, scale=wt1$scale, sig.test=1, 
                      sig.level=0.90, dof=length(cats1$Year)-wt1$scale, 
                      mother='morlet', lag1=0.0)
                      
                      

wt2$gws <- apply(wt2$power, 1, mean)

wt2$gws.sig <- wt.sig(d=cbind(cats1$Year,scale(cats1$precip)), dt=wt2$dt, scale=wt2$scale, sig.test=1, 
                      sig.level=0.90, dof=length(cats1$Year)-wt2$scale, 
                      mother='morlet', lag1=0.0)

wtc.AB$gws <- apply(wtc.AB$power, 1, mean)

wtc.AB

wtc.AB$gws.sig <- wt.sig(d=cbind(cats1$Year,scale(cats1$precip)), dt=wtc.AB$dt, scale=wtc.AB$scale, sig.test=1, 
                      sig.level=0.90, dof=length(cats1$Year)-wtc.AB$scale, 
                      mother='morlet', lag1=0.0)

par(mfrow=c(4,1),mar = c(4,4,1.5,4))

plot(wt1$gws, wt1$period, type="l",
     xlab="Global Wavelet Spectrum", ylab=NA,
     log="y", ylim=rev(range(wt1$period)), xlim=range(c(wt1$gws, wt1$gws.sig$signif)))
lines(wt2$gws, wt1$period, type="l",
     xlab="Global Wavelet Spectrum", ylab=NA,
     log="y", ylim=rev(range(wt1$period)), xlim=range(c(wt1$gws, wt1$gws.sig$signif)),lty=2)

plot(wt1$gws, wt1$period, type="l",
     xlab="Global Wavelet Spectrum", ylab=NA,
     log="y", ylim=rev(range(wt1$period)), xlim=range(c(wt1$gws, wt1$gws.sig$signif)))
lines(wt1$gws.sig$signif, wt1$period, lty=2, col='red')  

plot(wt2$gws, wt2$period, type="l",
     xlab="Global Wavelet Spectrum", ylab=NA,
     log="y", ylim=rev(range(wt2$period)), xlim=range(c(wt2$gws, wt2$gws.sig$signif)))
lines(wt2$gws.sig$signif, wt2$period, lty=2, col='red')  


plot(wtc.AB$gws, wtc.AB$period, type="l",
     xlab="Global Coherence", ylab=NA,
     log="y", ylim=rev(range(wtc.AB$period)), xlim=range(c(wtc.AB$gws, wtc.AB$gws.sig$signif)))


