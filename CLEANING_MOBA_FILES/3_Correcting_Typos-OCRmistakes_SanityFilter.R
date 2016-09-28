
#  manually correct wrong weight/height entries
#  (lost digits, zeros, decimal points)

#setwd("C:/Users/Jon/Dropbox/GIT/CLEANING_MOBA_FILES")

correctOCRandTypos=function(a,weights,lengths) { 

# manually picked lower/upper limits for reasonable values of height and weigth at each timepoint
wup= c( 7e3,10e3,12e3,15e3, 17e3,20e3,25e3,27e3,  35e3,45e3,55e3,80e3 )
wlw= c( 10,100,1e3,1e3,  2e3,3e3,5e3,6e3,  7e3,10e3,12e3,15e3 )
lup= c( 65,80,80,85,  90,100,110,120,  120,140,160,180)
llw= c( 10,20,30,45, 45,50,55,60, 65,70,90,100)

# collector of flags (insane values)
wmsk1 <<- matrix(0,nr=dim(a)[1] ,nc=length(weights)) # global assignment
lmsk1 <<- matrix(0,nr=dim(a)[1] ,nc=length(lengths)) # global --//--

par(mfrow=c(3,4))
for (time in 1:12){
print(paste(weights[time],sep=""))

w1=a[,weights[time]]
l1=a[,lengths[time]]

#####################################
# FIRST PART:  10xhigher and 10xlower

w1=a[,weights[time]]; l1=a[,lengths[time]]
#plot(w1~l1,main=paste(weights[time],sep=""))

# what is supposed to be a "normal value"?
LMD=median(l1[(l1>llw[time])&(l1<lup[time])],na.rm=T) # lenght (median)
LSD=sd(l1[(l1>llw[time])&(l1<lup[time])],na.rm=T) # lenght (standart deviations)
WMD=median(w1[(w1>wlw[time])&(w1<wup[time])],na.rm=T) # weight (median)
WSD=sd(w1[(w1>wlw[time])&(w1<wup[time])],na.rm=T) # weight (standard deviations)

#abline(v=c(LMD-4*LSD,LMD+4*LSD),lty=2)
#abline(h=c(WMD-4*WSD,WMD+4*WSD),lty=2)

#################################################
# deal with ten-times-higher-than-expected values
# (aka "the lost decimal point")
rm(lbad,wbad)
lbad=which((l1>((LMD-4*LSD)*10))&(l1<((LMD+4*LSD)*10)))
wbad=which((w1>((WMD-4*WSD)*10))&(w1<((WMD+4*WSD)*10)))
if (length(lbad)>0) l1[lbad]=l1[lbad]/10 # correction
if (length(wbad)>0) w1[wbad]=w1[wbad]/10 # correction
# report
print(paste(length(lbad)," 10x higher length corrections",sep=""))
print(paste(length(wbad)," 10x higher weight corrections",sep=""))

################################################
# deal with ten-times-lower-than-expected values
# (aka "the lost end-zero")
rm(lbad,wbad,bad)
terms=which(a$GestationalAge>=257)
# length
lbad=which((l1>((LMD-4*LSD)/10))&(l1<((LMD+4*LSD)/10)))
if (time==1) { bad=lbad[lbad %in% terms] } else { bad=lbad } # to prevent accidental correction of PTDs 
if (length(bad)>0) l1[bad]=l1[bad]*10 # correction
# weight
wbad=which((w1>((WMD-4*WSD)/10))&(w1<((WMD+4*WSD)/10)))
if (time==1) { bad=wbad[wbad %in% terms] } else { bad=wbad } # to prevent accidental correction of PTDs 
if (length(bad)>0) w1[bad]=w1[bad]*10 # correction
# report
print(paste(length(bad)," 10x lower length corrections",sep=""))
print(paste(length(bad)," 10x lower weight corrections",sep=""))

################################################
# deal with the 100-points lower values for length
# (aka "the lost first(digit)-ONE")
rm(lbad,wbad,bad)
if (time %in% c(9:12)) { # not applicable to earlier timepoints
  lbad=which( (l1>=10)&(l1>((LMD-4*LSD)-100))&(l1<100)&(l1<((LMD+4*LSD)-100)))
  if (length(lbad)>0) l1[lbad]=l1[lbad]+100 # correction
  print(paste(length(lbad)," 100cm lower length corrections",sep=""))
  }

if (time %in% c(9:12)) { # not applicable to earlier timepoints
  lbad=which( (l1>=150)&(l1>((LMD-4*LSD)+100))&(l1<200)&(l1<((LMD+4*LSD)+100)))
  if (length(lbad)>0) l1[lbad]=l1[lbad]-100 # correction
  print(paste(length(lbad)," 100cm higher length corrections",sep=""))
}

# no such correction is done for weight (too many digits)

a[,weights[time]] <<- w1 #(global assignment)
a[,lengths[time]] <<- l1 #(global assignment)


####################################################
## SECOND PART: MARK FOR DELETION THE INSANE VALUES
wmsk1[which(w1>wup[time]),time] <<- 1 #(gloabl)
wmsk1[which(w1<wlw[time]),time] <<- 1 #(gloabl)
lmsk1[which(l1>lup[time]),time] <<- 1 #(gloabl)
lmsk1[which(l1<llw[time]),time] <<- 1 #(gloabl)

####################################################
## SECOND PART: SIMPLY DELETE THE INSANE VALUES
a[ which(w1>wup[time]) ,weights[time]] <<- NA # (global)
a[ which(w1<wlw[time]) ,weights[time]] <<- NA # (global)
a[ which(l1>lup[time]) ,lengths[time]] <<- NA # (global)
a[ which(l1<llw[time]) ,lengths[time]] <<- NA # (global)

} # end of cycle


cat("wmsk1 and lmsk1 matrixes were created !")

} # end of function

