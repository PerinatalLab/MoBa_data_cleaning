

#  manually correct wrong weight/height entries
#  (lost digits, zeros, decimal points)

#setwd("C:/Users/Jon/Dropbox/GIT/CLEANING_MOBA_FILES")

correctOCRandTypos=function(a) { 

  
#  timepoint 1
time=1
w1=a[,weights[time]]
l1=a[,lengths[time]]
slc=which(a$GestationalAge>257) # slc= selection. only TERMS
MD=median(w1,na.rm=T); SD=sd(w1,na.rm=T)
#plot(w1~l1)
#plot(w1[slc]~l1[slc])

# very low (single-digit) length
#plot(w1~l1,xlim=c(0,10))
# --->>> nonsense values 

# weight single-digit values
ths=which((w1>(MD-4*SD)/1000)&(w1<((MD+4*SD)/1000))) # ths=these
use=ths[ths %in% slc]
COL=rep(1,dim(a)[1]); COL[use]=2
#plot(w1~l1,col=COL,ylim=c(0,30)); abline(h=c((MD-4*SD)/1000,(MD+4*SD)/1000))
# --->>> none

# weight double-digit values
ths=which((w1>((MD-4*SD)/100))&(w1<((MD+4*SD)/100))) # ths=these
use=ths[ths %in% slc]
COL=rep(1,dim(a)[1]); COL[use]=2
#plot(w1~l1,col=COL,ylim=c(0,100)); abline(h=c((MD-4*SD)/100,(MD+4*SD)/100))
# --->>> none

# weight tripple-digit values
ths=which((w1>((MD-4*SD)/10))&(w1<((MD+4*SD)/10))) # ths=these
use=ths[ths %in% slc]
COL=rep(1,dim(a)[1]); COL[use]=2
#plot(w1~l1,col=COL,ylim=c(0,1000)); abline(h=c((MD-4*SD)/10,(MD+4*SD)/10))
# SOLUTION
a[use,weights]
a[use,weights[time]]=NA # individual decision

##########################################
#  timepoint 2
time=2
w1=a[,weights[time]]
l1=a[,lengths[time]]
slc=which(a$GestationalAge>257) # slc= selection
MD=median(w1,na.rm=T); SD=sd(w1,na.rm=T)
#plot(w1~l1)

# very low (single-digit) length
# plot(w1~l1,xlim=c(0,10))
# --->>> nonsense values 

# weight three-digit values
ths=which((w1>((MD-4*SD)/10))&(w1<((MD+4*SD)/10))) # ths=these
use=ths[ths %in% slc]
COL=rep(1,dim(a)[1]); COL[use]=2
plot(w1~l1,col=COL,ylim=c(0,8000))
abline(h=c((MD-4*SD)/10,(MD+4*SD)/10))
abline(h=c((MD-4*SD),(MD+4*SD)),lty=2)
# SOLUTION
a[use,weights] # note that there are four datapoints, but two have no length values
a[use,weights[time]]=a[use,weights[time]]*10+5

##########################################
#  timepoint 3
time=3
w1=a[,weights[time]]
l1=a[,lengths[time]]
slc=which(a$GestationalAge>240) # slc= selection
MD=median(w1,na.rm=T); SD=sd(w1,na.rm=T)
LMD=median(l1,na.rm=T); LSD=sd(l1,na.rm=T) #(for length)
plot(w1~l1)

# weight single-digit values
ths=which((w1>(MD-4*SD)/1000)&(w1<((MD+4*SD)/1000))) # ths=these
use=ths[ths %in% slc]
COL=rep(1,dim(a)[1]); COL[use]=2
plot(w1~l1,col=COL,ylim=c(0,500)); abline(h=c((MD-4*SD)/1000,(MD+4*SD)/1000))
# SOLUTION
a[use,weights]
a[use,weights[time]]=NA # individual decision

# weight two-digit values
ths=which((w1>(MD-4*SD)/100)&(w1<((MD+4*SD)/100))) # ths=these
use=ths[ths %in% slc]
COL=rep(1,dim(a)[1]); COL[use]=2
plot(w1~l1,col=COL,ylim=c(0,1000)); abline(h=c((MD-4*SD)/100,(MD+4*SD)/100))
# --->>> none

# weight three-digit values
ths=which((w1>((MD-4*SD)/10))&(w1<((MD+4*SD)/10))) # ths=these
use=ths[ths %in% slc]
COL=rep(1,dim(a)[1]); COL[use]=2
plot(w1~l1,col=COL,ylim=c(0,1500)); abline(h=c((MD-4*SD)/10,(MD+4*SD)/10))
# SOLUTION
a[use,weights]
a[use,weights[time]]=a[use,weights[time]]*10+5

# problems with too low length
ths=which((l1>((LMD-4*LSD)/10))&(l1<((LMD+4*LSD)/10))) # ths=these
use=ths[ths %in% slc]
COL=rep(1,dim(a)[1]); COL[use]=2
plot(w1~l1,col=COL,xlim=c(0,100))
abline(v=c((LMD-4*LSD)/10,(LMD+4*LSD)/10))
abline(v=c((LMD-4*LSD),(LMD+4*LSD)),lty=2)
# SOLUTION
a[use,lengths]
a[use,lengths[time]]=a[use,lengths[time]]*10

##################################
#  timepoint 4
time=4
w1=a[,weights[time]]
l1=a[,lengths[time]]
MD=median(w1[(w1<20e3)&(w1>2e3)],na.rm=T); SD=sd(w1[(w1<20e3)&(w1>2e3)],na.rm=T)
LMD=median(l1[(l1<80)&(l1>45)],na.rm=T); LSD=sd(l1[(l1<80)&(l1>45)],na.rm=T)
plot(w1~l1)

# problems with too low length (ten times lower)
ths=which((l1>((LMD-4*LSD)/10))&(l1<((LMD+4*LSD)/10))) # ths=these
use=ths[ths %in% slc]
COL=rep(1,dim(a)[1]); COL[use]=2
plot(w1~l1,col=COL,xlim=c(0,100))
abline(v=c((LMD-4*LSD)/10,(LMD+4*LSD)/10))
abline(v=c((LMD-4*LSD),(LMD+4*LSD)),lty=2)
# SOLUTION
a[use,lengths]
a[use,lengths[time]]=a[use,lengths[time]]*10

# problems with too high weight (extra zero)
ths=which((w1>((MD-4*SD)*10))&(w1<((MD+4*SD)*10))) # ths=these
COL=rep(1,dim(a)[1]); COL[ths]=2
plot(w1~l1,col=COL)
abline(h=c((MD-4*SD)*10,(MD+4*SD)*10))
abline(h=c((MD-4*SD),(MD+4*SD)),lty=2)
# SOLUTION
a[ths,weights]
a[ths,weights[time]]=floor(a[ths,weights[time]]/10)

# problems with too low weight -  lost three zeros
ths=which((w1>(MD-4*SD)/1000)&(w1<((MD+4*SD)/1000))) # ths=these
COL=rep(1,dim(a)[1]); COL[ths]=2
plot(w1~l1,col=COL,ylim=c(0,1000)); abline(h=c((MD-4*SD)/1000,(MD+4*SD)/1000))
# SOLUTION
a[ths,weights]
a[ths,weights[time]]=NA #  delete that value

# problems with too low weight - lost two zeros
ths=which((w1>((MD-4*SD)/100))&(w1<((MD+4*SD)/100))) # ths=these
COL=rep(1,dim(a)[1]); COL[ths]=2
plot(w1~l1,col=COL,ylim=c(0,10000)); abline(h=c((MD-4*SD)/100,(MD+4*SD)/100))
# SOLUTION
a[ths,weights]
a[ths,weights[time]]=a[ths,weights[time]]*100+50

# problems with too low weight - lost single digit
ths=which((w1>((MD-4*SD)/10))&(w1<((MD+4*SD)/10))) # ths=these
COL=rep(1,dim(a)[1]); COL[ths]=2
plot(w1~l1,col=COL,ylim=c(0,20000))
abline(h=c((MD-4*SD)/10,(MD+4*SD)/10))
abline(h=c((MD-4*SD),(MD+4*SD)),lty=2)
# SOLUTION
a[ths,weights]
a[ths,weights[time]]=a[ths,weights[time]]*10+5


##############################
#  timepoint 5
time=5
w1=a[,weights[time]]
l1=a[,lengths[time]]
WMD=median(w1[(w1<2e4)&(w1>5e3)],na.rm=T); WSD=sd(w1[(w1<2e4)&(w1>5e3)],na.rm=T)
LMD=median(l1[(w1>50)],na.rm=T); LSD=sd(l1[(w1>50)],na.rm=T)
plot(w1~l1)

# problems with very low (single-digit) length
ths=which((l1>((LMD-4*LSD)/10))&(l1<((LMD+4*LSD)/10))) # ths=these
COL=rep(1,dim(a)[1]); COL[ths]=2
plot(w1~l1,col=COL)
abline(v=c((LMD-4*LSD)/10,(LMD+4*LSD)/10))
abline(v=c((LMD-4*LSD),(LMD+4*LSD)),lty=2)
# SOLUTION
a[ths,lengths]
a[ths,lengths[time]]=a[ths,lengths[time]]*10

# problems with very high (one zero too many) weight
ths=which((w1>((WMD-4*WSD)*10))&(w1<((WMD+4*WSD)*10))) # ths=these
COL=rep(1,dim(a)[1]); COL[ths]=2
plot(w1~l1,col=COL)
abline(h=c((WMD-4*WSD)*10,(WMD+4*WSD)*10))
abline(h=c((WMD-4*WSD),(WMD+4*WSD)),lty=2)
# SOLUTION
a[ths,weights]
a[ths,weights[time]]=floor(a[ths,weights[time]]/10)

# problems with too low weight -  lost three zeros
ths=which((w1>(WMD-4*WSD)/1000)&(w1<((WMD+4*WSD)/1000))) # ths=these
COL=rep(1,dim(a)[1]); COL[ths]=2
plot(w1~l1,col=COL,ylim=c(0,100)); abline(h=c((WMD-4*WSD)/1000,(WMD+4*WSD)/1000))
# SOLUTION
a[ths,weights]
a[ths,weights[time]]=NA #  delete that value

# problems with too low weight - lost two seros
ths=which((w1>(WMD-4*WSD)/100)&(w1<((WMD+4*WSD)/100))) # ths=these
COL=rep(1,dim(a)[1]); COL[ths]=2
plot(w1~l1,col=COL,ylim=c(0,5000)); abline(h=c((WMD-4*WSD)/100,(WMD+4*WSD)/100))
# SOLUTION
a[ths,weights]
a[ths,weights[time]]=NA #  delete that value

# problems with too low weight - lost one zero
ths=which((w1>((WMD-4*WSD)/10))&(w1<((WMD+4*WSD)/10))) # ths=these
COL=rep(1,dim(a)[1]); COL[ths]=2
plot(w1~l1,col=COL,ylim=c(0,10000)); abline(h=c((WMD-4*WSD)/10,(WMD+4*WSD)/10))
# SOLUTION
a[ths,weights]
a[ths,weights[time]]=a[ths,weights[time]]*10+5


################################
#  timepoint 6
time=6
w1=a[,weights[time]]
l1=a[,lengths[time]]
WMD=median(w1[(w1<2e4)&(w1>5e3)],na.rm=T); WSD=sd(w1[(w1<2e4)&(w1>5e3)],na.rm=T)
LMD=median(l1[(l1>60)],na.rm=T); LSD=sd(l1[(l1>60)],na.rm=T)
plot(w1~l1)

# problems with very low length (lost one digit)
ths=which((l1>((LMD-4*LSD)/10))&(l1<((LMD+4*LSD)/10))) # ths=these
COL=rep(1,dim(a)[1]); COL[ths]=2
plot(w1~l1,col=COL)
abline(v=c((LMD-4*LSD)/10,(LMD+4*LSD)/10))
abline(v=c((LMD-4*LSD),(LMD+4*LSD)),lty=2)
# SOLUTION
a[ths,lengths]
a[ths,lengths[time]]=a[ths,lengths[time]]*10

# problems with very high (one zero too many) weight
ths=which((w1>((WMD-4*WSD)*10))&(w1<((WMD+4*WSD)*10))) # ths=these
COL=rep(1,dim(a)[1]); COL[ths]=2
plot(w1~l1,col=COL)
abline(h=c((WMD-4*WSD)*10,(WMD+4*WSD)*10))
abline(h=c((WMD-4*WSD),(WMD+4*WSD)),lty=2)
# SOLUTION
a[ths,weights]
a[ths,weights[time]]=floor(a[ths,weights[time]]/10)

# problems with too low weight -  lost three zeros
ths=which((w1>(WMD-4*WSD)/1000)&(w1<((WMD+4*WSD)/1000))) # ths=these
COL=rep(1,dim(a)[1]); COL[ths]=2
plot(w1~l1,col=COL,ylim=c(0,100)); abline(h=c((WMD-4*WSD)/1000,(WMD+4*WSD)/1000))
# SOLUTION
a[ths,weights]
a[ths,weights[time]]=NA #  delete that value

# problems with too low weight - lost two seros
ths=which((w1>(WMD-4*WSD)/100)&(w1<((WMD+4*WSD)/100))) # ths=these
COL=rep(1,dim(a)[1]); COL[ths]=2
plot(w1~l1,col=COL,ylim=c(0,1000)); abline(h=c((WMD-4*WSD)/100,(WMD+4*WSD)/100))
# SOLUTION
a[ths,weights]
a[ths,weights[time]]=a[ths,weights[time]]*100+50

# problems with too low weight - lost one zero
ths=which((w1>((WMD-4*WSD)/10))&(w1<((WMD+4*WSD)/10))) # ths=these
COL=rep(1,dim(a)[1]); COL[ths]=2
plot(w1~l1,col=COL,ylim=c(0,20000))
abline(h=c((WMD-4*WSD)/10,(WMD+4*WSD)/10))
abline(h=c((WMD-4*WSD),(WMD+4*WSD)),lty=2)
# SOLUTION
a[ths,weights]
a[ths,weights[time]]=a[ths,weights[time]]*10


########################################
#  timepoint 7
time=7
w1=a[,weights[time]]
l1=a[,lengths[time]]
WMD=median(w1[(w1>5e3)&(w1<3e4)],na.rm=T); WSD=sd(w1[(w1>5e3)&(w1<3e4)],na.rm=T)
LMD=median(l1[(l1>60)&(l1<200)],na.rm=T); LSD=sd(l1[(l1>60)&(l1<200)],na.rm=T)
plot(w1~l1)

# problems with very high (one zero too many) weight
ths=which((w1>((WMD-4*WSD)*10))&(w1<((WMD+4*WSD)*10))) # ths=these
COL=rep(1,dim(a)[1]); COL[ths]=2
plot(w1~l1,col=COL)
abline(h=c((WMD-4*WSD)*10,(WMD+4*WSD)*10))
abline(h=c((WMD-4*WSD),(WMD+4*WSD)),lty=2)
# SOLUTION
a[ths,weights]
a[ths,weights[time]]=floor(a[ths,weights[time]]/10)

# problems with very HIGH length (lost one digit)
w1=a[,weights[time]]; l1=a[,lengths[time]]  # read re-newed values
ths=which((l1>((LMD-4*LSD)*10))&(l1<((LMD+5*LSD)*10))) # ths=these
COL=rep(1,dim(a)[1]); COL[ths]=2
plot(w1~l1,col=COL)
abline(v=c((LMD-4*LSD)*10,(LMD+4*LSD)*10))
abline(v=c((LMD-4*LSD),(LMD+4*LSD)),lty=2)
# SOLUTION
a[ths,lengths]
a[ths,lengths[time]]=floor(a[ths,lengths[time]]/10)

# problems with very low length (lost one digit)
w1=a[,weights[time]]; l1=a[,lengths[time]]  # read re-newed values
ths=which((l1>((LMD-4*LSD)/10))&(l1<((LMD+5*LSD)/10))) # ths=these
COL=rep(1,dim(a)[1]); COL[ths]=2
plot(w1~l1,col=COL)
abline(v=c((LMD-4*LSD)/10,(LMD+4*LSD)/10))
abline(v=c((LMD-4*LSD),(LMD+4*LSD)),lty=2)
# --->>> none

# problems with too low weight -  lost three zeros
ths=which((w1>(WMD-4*WSD)/1000)&(w1<((WMD+4*WSD)/1000))) # ths=these
COL=rep(1,dim(a)[1]); COL[ths]=2
plot(w1~l1,col=COL,ylim=c(0,100)); abline(h=c((WMD-4*WSD)/1000,(WMD+4*WSD)/1000))
# --->>> none

# problems with too low weight - lost two seros
ths=which((w1>(WMD-4*WSD)/100)&(w1<((WMD+4*WSD)/100))) # ths=these
COL=rep(1,dim(a)[1]); COL[ths]=2
plot(w1~l1,col=COL,ylim=c(0,1000)); abline(h=c((WMD-4*WSD)/100,(WMD+4*WSD)/100))
# --->>> none

# problems with too low weight - lost one zero
ths=which((w1>((WMD-4*WSD)/10))&(w1<((WMD+4*WSD)/10))) # ths=these
COL=rep(1,dim(a)[1]); COL[ths]=2
plot(w1~l1,col=COL,ylim=c(0,10000)); abline(h=c((WMD-4*WSD)/10,(WMD+4*WSD)/10))
# SOLUTION
a[ths,weights]
a[ths,weights[time]]=floor(a[ths,weights[time]]*10)

#########################################
#  timepoint 8
time=8
w1=a[,weights[time]]
l1=a[,lengths[time]]
WMD=median(w1[w1<30000],na.rm=T); WSD=sd(w1[w1<30000],na.rm=T)
LMD=median(l1[l1<200],na.rm=T); LSD=sd(l1[l1<200],na.rm=T)
plot(w1~l1)

# problems with very high length (gained one zero)
ths=which((l1>((LMD-4*LSD)*10))&(l1<((LMD+4*LSD)*10))) # ths=these
COL=rep(1,dim(a)[1]); COL[ths]=2
plot(w1~l1,col=COL)
abline(v=c((LMD-4*LSD)*10,(LMD+4*LSD)*10))
abline(v=c((LMD-4*LSD),(LMD+4*LSD)),lty=2)
# SOLUTION
a[ths,lengths]
a[ths,lengths[time]]=floor(a[ths,lengths[time]]/10)

# problems with very low length (lost one digit)
ths=which((l1>((LMD-4*LSD)/10))&(l1<((LMD+4*LSD)/10))) # ths=these
COL=rep(1,dim(a)[1]); COL[ths]=2
plot(w1~l1,col=COL,xlim=c(0,50)); abline(v=c((LMD-4*LSD)/10,(LMD+4*LSD)/10))
# SOLUTION
a[ths,lengths]
a[ths,lengths[time]]=floor(a[ths,lengths[time]]*10)

# problems with very high (one zero too many) weight
ths=which((w1>((WMD-4*WSD)*10))&(w1<((WMD+4*WSD)*10))) # ths=these
COL=rep(1,dim(a)[1]); COL[ths]=2
plot(w1~l1,col=COL)
abline(h=c((WMD-4*WSD)*10,(WMD+4*WSD)*10))
abline(h=c((WMD-4*WSD),(WMD+4*WSD)),lty=2)
# SOLUTION
a[ths,weights]
a[ths,weights[time]]=floor(a[ths,weights[time]]/10)

# problems with too low weight -  lost three zeros
w1=a[,weights[time]]; l1=a[,lengths[time]]
ths=which((w1>(WMD-4*WSD)/1000)&(w1<((WMD+4*WSD)/1000))) # ths=these
COL=rep(1,dim(a)[1]); COL[ths]=2
plot(w1~l1,col=COL,ylim=c(0,100)); abline(h=c((WMD-4*WSD)/1000,(WMD+4*WSD)/1000))
# --->>> none

# problems with too low weight - lost two zeros
ths=which((w1>(WMD-4*WSD)/100)&(w1<((WMD+4*WSD)/100))) # ths=these
COL=rep(1,dim(a)[1]); COL[ths]=2
plot(w1~l1,col=COL,ylim=c(0,1000)); abline(h=c((WMD-4*WSD)/100,(WMD+4*WSD)/100))
# --->>> none

# problems with too low weight - lost one zero
ths=which((w1>((WMD-4*WSD)/10))&(w1<((WMD+4*WSD)/10))) # ths=these
COL=rep(1,dim(a)[1]); COL[ths]=2
plot(w1~l1,col=COL,ylim=c(0,20000))
abline(h=c((WMD-4*WSD)/10,(WMD+4*WSD)/10))
abline(h=c((WMD-4*WSD),(WMD+4*WSD)),lty=2)
# --->>> none

#######################################
#  timepoint 9
time=9
w1=a[,weights[time]]
l1=a[,lengths[time]]
WMD=median(w1[(w1<30000)&(w1>3000)],na.rm=T); WSD=sd(w1[(w1<30000)&(w1>3000)],na.rm=T)
LMD=median(l1[(l1>60)&(l1<200)],na.rm=T); LSD=sd(l1[(l1>60)&(l1<200)],na.rm=T)
plot(w1~l1)

# problems with very high length (gained one zero)
ths=which((l1>((LMD-4*LSD)*10))&(l1<((LMD+4*LSD)*10))) # ths=these
COL=rep(1,dim(a)[1]); COL[ths]=2
plot(w1~l1,col=COL)
abline(v=c((LMD-4*LSD)*10,(LMD+4*LSD)*10))
abline(v=c((LMD-4*LSD),(LMD+4*LSD)),lty=2)
# SOLUTION
a[ths,lengths]
a[ths,lengths[time]]=floor(a[ths,lengths[time]]/10)

# problems with very low length (lost one digit)
w1=a[,weights[time]]; l1=a[,lengths[time]]
ths=which((l1>((LMD-4*LSD)/10))&(l1<((LMD+4*LSD)/10))) # ths=these
COL=rep(1,dim(a)[1]); COL[ths]=2
plot(w1~l1,col=COL,xlim=c(0,120))
abline(v=c((LMD-4*LSD)/10,(LMD+4*LSD)/10))
abline(v=c((LMD-4*LSD),(LMD+4*LSD)),lty=2)
# SOLUTION
a[ths,lengths]
a[ths,lengths[time]]=floor(a[ths,lengths[time]]*10)

# problems with very high (one zero too many) weight
ths=which((w1>((WMD-4*WSD)*10))&(w1<((WMD+4*WSD)*10))) # ths=these
COL=rep(1,dim(a)[1]); COL[ths]=2
plot(w1~l1,col=COL)
abline(h=c((WMD-4*WSD)*10,(WMD+4*WSD)*10))
abline(h=c((WMD-4*WSD),(WMD+4*WSD)),lty=2)
# SOLUTION
a[ths,weights]
a[ths,weights[time]]=floor(a[ths,weights[time]]/10) # (not one but two values)

# problems with too low weight -  lost three zeros
ths=which((w1>(WMD-4*WSD)/1000)&(w1<((WMD+4*WSD)/1000))) # ths=these
COL=rep(1,dim(a)[1]); COL[ths]=2
plot(w1~l1,col=COL,ylim=c(0,100)); abline(h=c((WMD-4*WSD)/1000,(WMD+4*WSD)/1000))
# --->>> none 

# problems with too low weight - lost two seros
ths=which((w1>(WMD-4*WSD)/100)&(w1<((WMD+4*WSD)/100))) # ths=these
COL=rep(1,dim(a)[1]); COL[ths]=2
plot(w1~l1,col=COL,ylim=c(0,5000)); abline(h=c((WMD-4*WSD)/100,(WMD+4*WSD)/100))
# SOLUTION
a[ths,weights]
a[ths,weights[time]]=floor(a[ths,weights[time]]*100+50) # (not one but two values)

# problems with too low weight - lost one zero
ths=which((w1>((WMD-4*WSD)/10))&(w1<((WMD+4*WSD)/10))) # ths=these
COL=rep(1,dim(a)[1]); COL[ths]=2
plot(w1~l1,col=COL,ylim=c(0,15000)); abline(h=c((WMD-4*WSD)/10,(WMD+4*WSD)/10))
# SOLUTION
a[ths,weights]
a[ths,weights[time]]=floor(a[ths,weights[time]]*10) # (not one but two values)


###########################
#  timepoint 10
time=10
w1=a[,weights[time]]
l1=a[,lengths[time]]
WMD=median(w1[(w1<50e3)&(w1>5e3)],na.rm=T); WSD=sd(w1[(w1<50e3)&(w1>5e3)],na.rm=T)
LMD=median(l1[(l1>60)&(l1<150)],na.rm=T); LSD=sd(l1[(l1>60)&(l1<150)],na.rm=T)
plot(w1~l1)

# problems with very high length (gained one zero)
ths=which((l1>((LMD-4*LSD)*10))&(l1<((LMD+4*LSD)*10))) # ths=these
COL=rep(1,dim(a)[1]); COL[ths]=2
plot(w1~l1,col=COL); abline(v=c((LMD-4*LSD)*10,(LMD+4*LSD)*10))
# --->>> none

# problems with very low length (lost one digit)
ths=which((l1>((LMD-4*LSD)/10))&(l1<((LMD+4*LSD)/10))) # ths=these
COL=rep(1,dim(a)[1]); COL[ths]=2
plot(w1~l1,col=COL,xlim=c(0,120)); abline(v=c((LMD-4*LSD)/10,(LMD+4*LSD)/10))
#   ->>>>>  only couple of values and no clear cluster

# problems with very high (one zero too many) weight
ths=which((w1>((WMD-4*WSD)*10))&(w1<((WMD+4*WSD)*10))) # ths=these
COL=rep(1,dim(a)[1]); COL[ths]=2
plot(w1~l1,col=COL)
abline(h=c((WMD-4*WSD)*10,(WMD+4*WSD)*10))
abline(h=c((WMD-4*WSD),(WMD+4*WSD)),lty=2)
# SOLUTION
a[ths,weights]
a[ths,weights[time]]=floor(a[ths,weights[time]]/10)

# problems with too low weight -  lost three zeros
ths=which((w1>(WMD-4*WSD)/1000)&(w1<((WMD+4*WSD)/1000))) # ths=these
COL=rep(1,dim(a)[1]); COL[ths]=2
plot(w1~l1,col=COL,ylim=c(0,100)); abline(h=c((WMD-4*WSD)/1000,(WMD+4*WSD)/1000))
# --->>>>  none

# problems with too low weight - lost two seros
ths=which((w1>(WMD-4*WSD)/100)&(w1<((WMD+4*WSD)/100))) # ths=these
COL=rep(1,dim(a)[1]); COL[ths]=2
plot(w1~l1,col=COL,ylim=c(0,1000)); abline(h=c((WMD-4*WSD)/100,(WMD+4*WSD)/100))
# --->>>>  none

# problems with too low weight - lost one zero
ths=which((w1>((WMD-4*WSD)/10))&(w1<((WMD+4*WSD)/10))) # ths=these
COL=rep(1,dim(a)[1]); COL[ths]=2
plot(w1~l1,col=COL,ylim=c(0,15000)); abline(h=c((WMD-4*WSD)/10,(WMD+4*WSD)/10))
# SOLUTION
a[ths,weights]
a[ths,weights[time]]=floor(a[ths,weights[time]]*10)

################################
#  timepoint 11
time=11
w1=a[,weights[time]]
l1=a[,lengths[time]]
WMD=median(w1[(w1<50e3)&(w1>5e3)],na.rm=T); WSD=sd(w1[(w1<50e3)&(w1>5e3)],na.rm=T)
LMD=median(l1[(l1>60)&(l1<160)],na.rm=T); LSD=sd(l1[(l1>60)&(l1<150)],na.rm=T)
plot(w1~l1)

# problems with very high length (gained one zero)
ths=which((l1>((LMD-4*LSD)*10))&(l1<((LMD+4*LSD)*10))) # ths=these
COL=rep(1,dim(a)[1]); COL[ths]=2
plot(w1~l1,col=COL); abline(v=c((LMD-4*LSD)*10,(LMD+4*LSD)*10))
# --->>>>  none

# problems with very low length (lost one digit)
ths=which((l1>((LMD-4*LSD)/10))&(l1<((LMD+4*LSD)/10))) # ths=these
COL=rep(1,dim(a)[1]); COL[ths]=2
plot(w1~l1,col=COL,xlim=c(0,100)); abline(v=c((LMD-4*LSD)/10,(LMD+4*LSD)/10))
# --->>>>  none 

# problems with very high (one zero too many) weight
ths=which((w1>((WMD-4*WSD)*10))&(w1<((WMD+4*WSD)*10))) # ths=these
COL=rep(1,dim(a)[1]); COL[ths]=2
plot(w1~l1,col=COL)
abline(h=c((WMD-4*WSD)*10,(WMD+4*WSD)*10))
abline(h=c((WMD-4*WSD),(WMD+4*WSD)),lty=2)
# SOLUTION
a[ths,weights]
a[ths,weights[time]]=floor(a[ths,weights[time]]/10)  # not on the plot

# problems with too low weight -  lost three zeros
ths=which((w1>(WMD-4*WSD)/1000)&(w1<((WMD+4*WSD)/1000))) # ths=these
COL=rep(1,dim(a)[1]); COL[ths]=2
plot(w1~l1,col=COL,ylim=c(0,100)); abline(h=c((WMD-4*WSD)/1000,(WMD+4*WSD)/1000))
# --->>>>  none

# problems with too low weight - lost two seros
ths=which((w1>(WMD-4*WSD)/100)&(w1<((WMD+4*WSD)/100))) # ths=these
COL=rep(1,dim(a)[1]); COL[ths]=2
plot(w1~l1,col=COL,ylim=c(0,1000)); abline(h=c((WMD-4*WSD)/100,(WMD+4*WSD)/100))
# --->>>>  none

# problems with too low weight - lost one zero
ths=which((w1>((WMD-4*WSD)/10))&(w1<((WMD+4*WSD)/10))) # ths=these
COL=rep(1,dim(a)[1]); COL[ths]=2
plot(w1~l1,col=COL,ylim=c(0,30000)); abline(h=c((WMD-4*WSD)/10,(WMD+4*WSD)/10))
# SOLUTION
a[ths,weights]
a[ths,weights[time]]=floor(a[ths,weights[time]]*10)  # most of the points are not on the plot 


##########################################
#  timepoint 12
time=12
w1=a[,weights[time]]
l1=a[,lengths[time]]
WMD=median(w1[(w1<70e3)&(w1>10e3)],na.rm=T); WSD=sd(w1[(w1<70e3)&(w1>10e3)],na.rm=T)
LMD=median(l1[(l1>100)&(l1<170)],na.rm=T); LSD=sd(l1[(l1>100)&(l1<170)],na.rm=T)
plot(w1~l1)

# problems with very high length (gained one zero)
ths=which((l1>((LMD-4*LSD)*10))&(l1<((LMD+4*LSD)*10))) # ths=these
COL=rep(1,dim(a)[1]); COL[ths]=2
plot(w1~l1,col=COL); abline(v=c((LMD-4*LSD)*10,(LMD+4*LSD)*10))
# --->>>>  none

# problems with very low length (lost one digit)
ths=which((l1>((LMD-4*LSD)/10))&(l1<((LMD+4*LSD)/10))) # ths=these
COL=rep(1,dim(a)[1]); COL[ths]=2
plot(w1~l1,col=COL,xlim=c(0,100)); abline(v=c((LMD-4*LSD)/10,(LMD+4*LSD)/10))
# SOLUTION
a[ths,lengths]
a[ths,lengths[time]]=floor(a[ths,lengths[time]]*10)

# problems with very high (one zero too many) weight
ths=which((w1>((WMD-4*WSD)*10))&(w1<((WMD+4*WSD)*10))) # ths=these
COL=rep(1,dim(a)[1]); COL[ths]=2
plot(w1~l1,col=COL)
abline(h=c((WMD-4*WSD)*10,(WMD+4*WSD)*10))
abline(h=c((WMD-4*WSD),(WMD+4*WSD)),lty=2)
# --->>>>  none

# problems with too low weight -  lost three zeros
ths=which((w1>(WMD-4*WSD)/1000)&(w1<((WMD+4*WSD)/1000))) # ths=these
COL=rep(1,dim(a)[1]); COL[ths]=2
plot(w1~l1,col=COL,ylim=c(0,100)); abline(h=c((WMD-4*WSD)/1000,(WMD+4*WSD)/1000))
# --->>>>  none

# problems with too low weight - lost two seros
ths=which((w1>(WMD-4*WSD)/100)&(w1<((WMD+4*WSD)/100))) # ths=these
COL=rep(1,dim(a)[1]); COL[ths]=2
plot(w1~l1,col=COL,ylim=c(0,1000)); abline(h=c((WMD-4*WSD)/100,(WMD+4*WSD)/100))
# --->>>>  none

# problems with too low weight - lost one zero
ths=which((w1>((WMD-4*WSD)/10))&(w1<((WMD+4*WSD)/10))) # ths=these
COL=rep(1,dim(a)[1]); COL[ths]=2
plot(w1~l1,col=COL,ylim=c(0,30000))
abline(h=c((WMD-4*WSD)/10,(WMD+4*WSD)/10))
abline(h=c((WMD-4*WSD),(WMD+4*WSD)),lty=2)
# SOLUTION
a[ths,weights]
a[ths,weights[time]]=floor(a[ths,weights[time]]*10)

###############################################################
#  only for those values that cannot be plotted (no matching weight/height value)

# manually picked lower-upper limits for reasonable values of height and weigth at each timepoint
wup= c( 6e3,10e3,12e3,15e3, 16e3,17e3,20e3,22e3,  30e3,32e3,50e3,60e3 )
wlw= c( 0,0,1e3,1e3,  2e3,3e3,5e3,7e3,  8e3,10e3,12e3,15e3 )
lup= c( 65,80,80,85,  90,100,110,120,  150,150,150,160)
llw= c( 5,20,40,45, 45,50,55,60, 70,80,100,120)

# dealing with nureasonably high values of weigth and height
par(mfrow=c(3,4))
for (time in 1:12){
print(paste(weights[time],sep=""))

w1=a[,weights[time]]; l1=a[,lengths[time]]
#plot(w1~l1,main=paste(weights[time],sep=""))

LMD=median(l1[(l1>llw[time])&(l1<lup[time])],na.rm=T)
LSD=sd(l1[(l1>llw[time])&(l1<lup[time])],na.rm=T)
WMD=median(w1[(w1>wlw[time])&(w1<wup[time])],na.rm=T)
WSD=sd(w1[(w1>wlw[time])&(w1<wup[time])],na.rm=T)

abline(v=c(LMD-4*LSD,LMD+4*LSD),lty=2)
abline(h=c(WMD-4*WSD,WMD+4*WSD),lty=2)

lbad=which((l1>((LMD-4*LSD)*10))&(l1<((LMD+4*LSD)*10)))
wbad=which((w1>((WMD-4*WSD)*10))&(w1<((WMD+4*WSD)*10)))

if (length(lbad)>0) l1[lbad]=l1[lbad]/10
if (length(wbad)>0) w1[wbad]=w1[wbad]/10

a[,weights[time]] <- w1
a[,lengths[time]] <- l1

print(paste(length(lbad)," length corrections",sep=""))
print(paste(length(wbad)," weight corrections",sep=""))


} # end of cycle

.. should it also be implemented for the 10-times lower weight and height?




} # end of function

