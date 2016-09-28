

setwd("C:/Users/Jon/Dropbox/GIT/CLEANING_MOBA_FILES")

PlotOverlapingFlagsLength=function(a,weights,lengths) {
  
  cat("
  color
  GREEN:  time-based  (100)
  BLUE: univar in respect to GA (010)
  MAGENTA:  residuals (w-l) in respect to GA  (001)
  RED: combinations  (111,110,011,101)
  
  pch
  FULL BUBBLE:  c(100,010,001)
  +:  time-based & univar in respect to GA  (110)
  X: univar & residuals (w-l) in respect to GA  (011)
  RHOMB: time-based & residuals (w-l) in respect to GA (101)
  FULL: all possible flags (111)
  ")
  
  
par(mfrow=c(3,4))
for (time in 1:12) {
  w=a[,weights[time]]
  l=a[,lengths[time]]
  lt=paste(lmsk2[,time],lmsk3[,time],pmsk3[,time],sep="") # length test, pmask = residuals
  
  COL=rep(8,length(w))
  COL[lt=="100"]=3 # GREEN:  time-based
  COL[lt=="010"]=4 # BLUE: standard devs in respect to GA
  COL[lt=="001"]=6 # MAGENTA:  residuals (w-l) in respect to GA
  COL[lt %in% c("111","110","011","101")]=2 # RED: combinations
  
  PCH=rep(1,length(w))
  PCH[lt %in% c("100","010","001")]=19 # FULL
  PCH[lt=="110"]=3 # +:  time-based & standard devs in respect to GA
  PCH[lt=="011"]=4 # X: standard devs & residuals (w-l) in respect to GA
  PCH[lt=="101"]=5 # RHOMB: time-based & residuals (w-l) in respect to GA
  PCH[lt=="111"]=19 # FULL: all possible flags
  
  ORD=rep(2,length(w))
  ORD[which(lt=="000")]=1
  
  dat=data.frame(w=w,l=l,COL=COL,PCH=PCH,ORD=ORD)
  dat=dat[order(dat$ORD),] # to prevent hidden colored values
  plot(dat$w~dat$l,col=dat$COL,pch=dat$PCH,main=time,xlab=lengths[time],ylab=weights[time])
  }


} # end of function
