##Code for change-point models for Pepi, Holyoak & Karban 2021, Ecology Letters


library(AICcmodavg)
#library(remotes)
#install_github("MatthieuStigler/seglm")
library(seglm)
library(tsDyn)

cats1<-read.csv('bodega.cats3.csv')
str(cats1)

cats1$logmean<-log(cats1$cat.mean)
cats1$logmeanlag1<-NA
cats1$logmeanlag2<-NA
cats1$logmeanlag1[2:34]<-cats1$logmean[1:33]
cats1$logmeanlag2[3:34]<-cats1$logmean[1:32]
cats1$preciplag1<-NA
cats1$preciplag1[2:34]<-scale(cats1$precip[1:33])
bmr<-log(cats1$cat.mean)
year<-cats1$Year

setar3<-setar(bmr,thDelay=0,m=2,thVar=year,nthresh=1)
summary(setar3)

ar2<-linear(bmr,m=2)
summary(ar2)

AIC(ar2)
AIC(setar3)
AIC(ar2)-AIC(setar3)


#SEGLM

cats1lag2<-cats1[3:34,]
str(cats1lag2)
seglm1<-seglm_lm(logmean~logmeanlag1+logmeanlag2+preciplag1,data=cats1lag2, th_var = "Year", nthresh =1)
seglm1
summary(seglm1)

lm1<-lm(logmean~logmeanlag1+logmeanlag2+preciplag1,data=cats1lag2)
summary(lm1)
AIC(seglm1)
AIC(lm1)



