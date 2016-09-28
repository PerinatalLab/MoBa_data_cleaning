


library(SDMTools)
library(hexbin)


setwd("C:/Users/Jon/Dropbox/GIT/CLEANING_MOBA_FILES")

# read-in data
name="~/SpiderOak Hive/CLEANING_MOBA_FILES/PDB1555_SV_INFO_Status_MFRsingleton_Q1-8aar_far_koffein_growth_v8_140811.csv"
a=read.csv(name,header=T,sep=",")
colnames(a)


# for weight-length cleaning
weights=c("WeightBirthMFR","Weight6w","Weight3m","Weight5_6m",
          "Weight8m","Weight1y","Weight15_18mQ5","Weight2y",
          "Weight3y","Weight5y","Weight7y","Weight8y")
lengths=c( "LengthBirthMFR","Length6w","Length3m","Length5_6m",
           "Length8m","Length1y","Length15_18mQ5","Length2y",
           "Length3y","Length5y","Length7y","Length8y")

# for head-circumferance cleaning
head=c("HeadBirth","Head6w","Head3m","Head5_6m","Head8m","Head1y")
weig=c("WeightBirthMFR","Weight6w","Weight3m","Weight5_6m","Weight8m","Weight1y")
leng=c("LengthBirthMFR","Length6w","Length3m","Length5_6m","Length8m","Length1y")


# investigate the number of nonmissing values at each timepooint for weight and length
b=a
re=NULL
for (i in 1:12) {
  lpr=round(sum(!is.na(b[,lengths[i]]))/dim(b)[1],3)
  lnr=sum(!is.na(b[,lengths[i]]))
  wpr=round(sum(!is.na(b[,weights[i]]))/dim(b)[1],3)
  wnr=sum(!is.na(b[,weights[i]]))
  re=rbind(re,data.frame(time=i,length=lengths[i],weight=weights[i],lng.nr.prez=lnr,
                         lng.prc.prez=lpr,wgh.nr.prez=wnr,wgh.prc.prez=wpr))
}
re
#write.table(re,"ManyMisisngLengthValuesAtTimepoint11.txt",row.names=F,col.names=T,quote=F,sep="\t")



######################################################################################
###################################################################     CHILDREN FIRST

###############################
## Stage 0) PLOT UNCLEANED DATA
jpeg("./Childs_HghWghCrc/0_Before_cleaning.jpeg",
     width = 30,height = 25, quality=100,units="cm", res=150)
par(mfrow=c(3,4))
for (i in 1:12) plot(a[,weights[i]]~a[,lengths[i]],xlab=lengths[i],ylab=weights[i],main=i)
dev.off()

#################################
## Stage 1) MERGE DUBLICATED DATA

# IN Q5 and Q6
setwd("C:/Users/Jon/Dropbox/GIT/CLEANING_MOBA_FILES")
source("1-1_enrich_15-18mQ5_with_Q6.R") # read-in code and functions
refillQ5withQ6(a) # function imports q6 values when q5 values are wrong or missing

# IN MFR and Q4
setwd("C:/Users/Jon/Dropbox/GIT/CLEANING_MOBA_FILES")
source("1-2_enrich_MFR_with_Q4.R") # read-in code and functions
refillMFRwithQ4(a) # function imports q4 values when MFR values are missing

# PREVIEW CHANGES
jpeg("./Childs_HghWghCrc/1_after_merging_duplicated_date.jpeg",
     width = 30,height = 25, quality=100,units="cm", res=150)
par(mfrow=c(3,4))
for (i in 1:12) plot(a[,weights[i]]~a[,lengths[i]],xlab=lengths[i],ylab=weights[i],main=i)
dev.off()

####################################
## Stage 2) GESTATIONAL AGE PROBLEMS

# cean
bad=which((a$GestationalAge<154)|(a$GestationalAge>307)|(is.na(a$GestationalAge)))
for (i in 1:12) a[bad,lengths[i]]= NA
for (i in 1:12) a[bad,weights[i]]= NA
for (i in 1:6)  a[bad,head[i]]= NA
a[bad,"GestationalAge"]= NA

# preview
jpeg("./Childs_HghWghCrc/2_after_gestational_age_cleaning.jpeg",
     width = 30,height = 25, quality=100,units="cm", res=150)
par(mfrow=c(3,4))
for (i in 1:12) plot(a[,weights[i]]~a[,lengths[i]],xlab=lengths[i],ylab=weights[i],main=i)
dev.off()

#######################################
## Stage 3)    TYPOS, INSANE VALUES

# DEAL WITH TYPOS AND OCR PROBLEMS, INSANE values
setwd("C:/Users/Jon/Dropbox/GIT/CLEANING_MOBA_FILES")
source("3_Correcting_Typos-OCRmistakes_SanityFilter.R") # read-in code and functions
correctOCRandTypos(a,weights,lengths) # function imports q4 values when MFR values are missing

# preview
jpeg("./Childs_HghWghCrc/3_after_insane-deletion_and_typo-correction.jpeg",
     width = 30,height = 25, quality=100,units="cm", res=150)
par(mfrow=c(3,4))
for (i in 1:12) plot(a[,weights[i]]~a[,lengths[i]],xlab=lengths[i],ylab=weights[i],main=i)
dev.off()


################################
## Stage 4) TIME-BASED  cleaning


# FIRST RUN
setwd("C:/Users/Jon/Dropbox/GIT/CLEANING_MOBA_FILES")
source("4-1_time-based_outlier_detection.R") # read-in code and functions
TimeBasedOutlierDetection(a,weights,lengths)

# SECOND RUN (needed)
setwd("C:/Users/Jon/Dropbox/GIT/CLEANING_MOBA_FILES")
source("4-1_time-based_outlier_detection.R") # read-in code and functions
TimeBasedOutlierDetection(a,weights,lengths)

# PREVIEW
jpeg("./Childs_HghWghCrc/4_after_initial_time-besed_cleaning.jpeg",
     width = 30,height = 25, quality=100,units="cm", res=150)
par(mfrow=c(3,4))
for (i in 1:12) plot(a[,weights[i]]~a[,lengths[i]],xlab=lengths[i],ylab=weights[i],main=i)
dev.off()

# see overlaps between flag-files
jpeg("./Childs_HghWghCrc/4_overlaps_between_various_time-based_flags.jpeg",
     width = 30,height = 25, quality=100,units="cm", res=150)
setwd("C:/Users/Jon/Dropbox/GIT/CLEANING_MOBA_FILES")
source("4-2_time-based_LW-overlaps.R") # read-in code and functions
TimeBasedWeightLengthOverlaps(a,weights,lengths)
dev.off()

## Explore problematic remaining values
table(lerrs,lmaxs)
table(werrs,wmaxs)


##########################################
## Stage 5) RESIDUALS AND UNIVAR OUTLIERS

# DEAL WITH outliers USING THE GA-stratified univariate and bivariate outliers (also creates lots of pdfs)
setwd("C:/Users/Jon/Dropbox/GIT/CLEANING_MOBA_FILES")
source("5-1_GA-stratified_residual_and_univar_outliers.R") # read-in code and functions
GAbasedUnivarBivarOutlierDetection(a,weights,lengths)


# PREVIEW ALL CURRENT DELETIONS (how clean is the data)
jpeg("./Childs_HghWghCrc/5_if_all_flags_get_deleted.jpeg",
     width = 30,height = 25, quality=100,units="cm", res=150)
par(mfrow=c(3,4))
  for (time in 1:12) {
  w=a[,weights[time]]
  l=a[,lengths[time]]
wt=c(wmsk1[,time]+wmsk2[,time]+wmsk3[,time]+pmsk3[,time])
lt=c(lmsk1[,time]+lmsk2[,time]+lmsk3[,time]+pmsk3[,time])
t=wt+lt
plot(w[which(t==0)]~l[which(t==0)],main=time,xlab=lengths[time],ylab=weights[time])
  }
dev.off()

#  SEE THE OVERLAPS
jpeg("./Childs_HghWghCrc/5_overlaps_between_flags_length.jpeg",
     width = 30,height = 25, quality=100,units="cm", res=150)
setwd("C:/Users/Jon/Dropbox/GIT/CLEANING_MOBA_FILES")
source("5_OverlappingFlags_length-based.R") # read-in code and functions
PlotOverlapingFlagsLength(a,weights,lengths)
dev.off()

jpeg("./Childs_HghWghCrc/5_overlaps_between_flags_weight.jpeg",
     width = 30,height = 25, quality=100,units="cm", res=150)
setwd("C:/Users/Jon/Dropbox/GIT/CLEANING_MOBA_FILES")
source("5_OverlappingFlags_weight-based.R") # read-in code and functions
PlotOverlapingFlagsWeight(a,weights,lengths)
dev.off()

##  STEP 6A)  CLEANING UNIVAR PROBLEMS
for (time in 1:12) {
  a[,lengths[time]][which(lmsk3[,time]>0)]=NA
  a[,weights[time]][which(wmsk3[,time]>0)]=NA
}


## STEP 6B)  CLEANING ALL TIME-BASED PROBLEMS
for (time in 1:12) {
  a[,lengths[time]][which(lmsk2[,time]>0)]=NA
  a[,weights[time]][which(wmsk2[,time]>0)]=NA
}


## STEP 6C)  CLEANING RESIDUALS PROBL
for (time in 1:12) {
  a[,lengths[time]][which(pmsk3[,time]>0)]=NA
  a[,weights[time]][which(pmsk3[,time]>0)]=NA
}

# save a backup copy: 
write.table(a,"C:/Users/Jon/Desktop/PDB1555/saved_a_file_WghHgh-OK.txt",row.names=F,col.names=T,
            sep="\t",dec=".",quote=F)



## STEP 7)  SAVE UNIVAR HISTOGRAMS FOR LENGTH
jpeg("./Childs_HghWghCrc/7_LENGTH_cleaned_histograms.jpeg",
     width = 30,height = 25, quality=100,units="cm", res=150)
par(mfrow=c(3,4))
for (i in 1:12)  { hist(a[,lengths[i]],breaks=100,col="grey",main=i)
                   abline(v=c(min(a[,lengths[i]],na.rm=T),max(a[,lengths[i]],na.rm=T)),lty=2,col="red")
}
dev.off()

## STEP 8)   SAVE UNIVAR HISTOGRAMS FOR WEIGHT
jpeg("./Childs_HghWghCrc/8_WEIGHT_cleaned_histograms.jpeg",
     width = 30,height = 25, quality=100,units="cm", res=150)
par(mfrow=c(3,4))
for (i in 1:12)  { hist(a[,weights[i]],breaks=100,col="grey",main=i)
                   abline(v=c(min(a[,weights[i]],na.rm=T),max(a[,weights[i]],na.rm=T)),lty=2,col="red")
}
dev.off()




####################################################################
###############    HEAD CIRCUMFERANCE    ###########################

## STEP 9) PREVIEW HEAD CIRCUMFERANCE DAA BEFOR CLEANING
jpeg("./Childs_HghWghCrc/9_head-circumf_before_cleaning.jpeg",
     width = 30,height = 25, quality=100,units="cm", res=150)
par(mfrow=c(2,3))
for (i in 1:6) plot(a[,weig[i]]~a[,head[i]],main=i,xlab=head[i],ylab=weig[i])
dev.off()




## STEP 10) UNREASONABLE VALUES AT HEAD CIRCUMFERANCE
setwd("C:/Users/Jon/Dropbox/GIT/CLEANING_MOBA_FILES")
source("10_imposible_values_Head.R") # read-in code and functions
correctOCRandTyposHEAD(a,weig,leng,head) 
  
## STEP 11) TIEM-BASED
setwd("C:/Users/Jon/Dropbox/GIT/CLEANING_MOBA_FILES")
source("11_time-based_outlier_detection_HEAD.R") # read-in code and functions
TimeBasedOutlierDetectionHEAD(a,head)

## STEP 12) BASED ON UNIVAR AND RESIDUALS
setwd("C:/Users/Jon/Dropbox/GIT/CLEANING_MOBA_FILES")
source("12_GAstratified_residual_and_univar_outliers_HEAD.R") # read-in code and functions
GAbasedUnivarBivarOutlierDetectionHEAD(a,weig,leng,head)




## STEP 13) visualize exclusions
jpeg("./Childs_HghWghCrc/13_HEAD_all_flags_together.jpeg",
     width = 30,height = 25, quality=100,units="cm", res=150)
par(mfrow=c(2,3))
for (i in 1:6) {
 msks=hmsk1[,i]+hmsk2[,i]+hmsk3[,i]+hwmsk3[,i]
 PCH=rep(1,length(msks)); PCH[msks>0]=19
 dada=data.frame(w=a[,weig[i]],h=a[,head[i]],p=PCH,m=msks)
 dada=dada[order(dada$m),]
 plot(dada$w~dada$h,main=i,xlab=head[i],ylab=weig[i],col=dada$m+1,pch=dada$p)
}
dev.off()





##  STEP 14)  CLEANING INCREDIBLE VALUES
for (i in 1:6) print(length(which(hmsk1[,i]>0))) # report the number of problems
for (i in 1:6) a[,head[i]][which(hmsk2[,i]>0)]=NA # make problematic values missing

## STEP 15)  CLEANING UNIVAR PROBLEMS
for (i in 1:6) print(length(which(hmsk3[,i]>0))) # report the number of problems
for (i in 1:6) a[,head[i]][which(hmsk3[,i]>0)]=NA # make problematic values missing

## STEP 16)  CLEANING ALL TIME-BASED PROBLEMS
for (i in 1:6) print(length(which(hmsk2[,i]>0))) # report the number of problems
for (i in 1:6) a[,head[i]][which(hmsk2[,i]>0)]=NA # make problematic values missing

## STEP 17)  CLEANING RESIDUALS PROBL
for (i in 1:6) print(length(which(hwmsk3[,i]>0))) # report the number of problems
for (i in 1:6)   a[,head[i]][which(hwmsk3[,i]>0)]=NA

## STEP 18)  PLOT ALL THE REMAINING VALUES
jpeg("./Childs_HghWghCrc/18_HEAD_cleaned_against_weight.jpeg",
     width = 30,height = 25, quality=100,units="cm", res=150)
par(mfrow=c(2,3))
for (i in 1:6)   plot( a[,weig[i]]~a[,head[i]],main=i )  # with WEIGHT
dev.off()

jpeg("./Childs_HghWghCrc/18_HEAD_cleaned_against_height.jpeg",
     width = 30,height = 25, quality=100,units="cm", res=150)
par(mfrow=c(2,3))
for (i in 1:6)   plot( a[,leng[i]]~a[,head[i]],main=i )  # with LENGTH
dev.off()

jpeg("./Childs_HghWghCrc/18_HEAD_cleaned_histograms.jpeg",
     width = 30,height = 25, quality=100,units="cm", res=150)
par(mfrow=c(2,3))
for (i in 1:6)  { hist(a[,head[i]],breaks=100,col="grey",main=i)
abline(v=c(min(a[,head[i]],na.rm=T),max(a[,head[i]],na.rm=T)),lty=2,col="red")
                }
dev.off()


#################################################################################
###########    DEALING WITH PROBLEMS IN PARENTAL MEASUREMENTS  ##################
#################################################################################




## STEP 20) 
setwd("C:/Users/Jon/Dropbox/GIT/CLEANING_MOBA_FILES")
source("20_maternal-paternal-height-and-weight_problems.R") # read-in code and functions
matpatweigheig(a)


#total misisng values:
sum(is.na(a))  # 1 727 392, while before it was  1 696 274


###   RENAME THE VARIABLES WITH THE "c" prefix

for (clnm in weights) colnames(a)[which(colnames(a) == clnm)] = paste("c_",clnm,sep="")
for (clnm in lengths) colnames(a)[which(colnames(a) == clnm)] = paste("c_",clnm,sep="")
for (clnm in head)    colnames(a)[which(colnames(a) == clnm)] = paste("c_",clnm,sep="")
others=c("PaternalHeight", "MaternalHeight", "MaternalWeight", "GestationalAge")
for (clnm in others)    colnames(a)[which(colnames(a) == clnm)] = paste("c_",clnm,sep="")
colnames(a)[1]="PREG_ID_1555"




# save a backup copy: 
write.table(a,"C:/Users/Jon/Desktop/PDB1555/saved_a_file_WghHghHead_Parents-OK.txt",row.names=F,col.names=T,
            sep="\t",dec=".",quote=F,na="")

dim(a)

colnames(a)







######## BELOW:   SANDBOX  (do not use)



# what is supposed to be a "normal value"
MHMD=median(mh[(mh>mhlw)&(mh<mhup)],na.rm=T) # maternal height median
MHSD=sd(mh[(mh>mhlw)&(mh<mhup)],na.rm=T) # maternal height standard deviations

MWMD=median(mw[(mw>mwlw)&(mw<mwup)],na.rm=T) # maternal weight median
MWSD=sd(mw[(mw>mwlw)&(mw<mwup)],na.rm=T) # maternal weight standard deviations

PHMD=median(ph[(ph>phlw)&(ph<phup)],na.rm=T) # paternal height median
PHSD=sd(ph[(ph>phlw)&(ph<phup)],na.rm=T) # paternal height standard deviations

# visualize
plot(mw~mh)
mh1=MHMD-4*MHSD; mh2=MHMD+4*MHSD
mw1=MWMD-3*MWSD; mw2=MWMD+5*MWSD
ph1=PHMD-4*PHSD; ph2=PHMD+5*PHSD

abline(v=c(mh1,mh2),lty=2)
abline(h=c(mw1,mw2),lty=2)
abline(v=c(mh1*10,mh2*10),lty=2,col="red")
abline(h=c(mw1*10,mw2*10),lty=2,col="red")
abline(v=c(mh1/10,mh2/10),lty=2,col="green")
abline(h=c(mw1/10,mw2/10),lty=2,col="green")
abline(v=c(mh1-100,mh2-100),lty=2,col="blue")
#mh[which((mh<140)&(mh>100))]

# deal with ten-times-higher-than-expected values (aka "the lost decimal point")
mhbad=which((mh>(mh1*10))&(mh<(mh2*10)))
mwbad=which((mw>(mw1*10))&(mw<(mw2*10)))
phbad=which((ph>(ph1*10))&(ph<(ph2*10)))

points(mw[mwbad]~mh[mwbad],col="red",pch=19)
points(mw[mhbad]~mh[mhbad],col="red",pch=19)

if (length(mhbad)>0) mh[mhbad]=mh[mhbad]/10 # correction
if (length(mwbad)>0) mw[mwbad]=mw[mwbad]/10 # correction
if (length(phbad)>0) ph[phbad]=ph[phbad]/10 # correction

print(paste(length(mhbad)," 10x higher maternal height corrections",sep=""))
print(paste(length(mwbad)," 10x higher maternal weight corrections",sep=""))
print(paste(length(phbad)," 10x higher paternal height corrections",sep=""))

rm(mhbad,mwbad,phbad)

# deal with ten-times-lower-than-expected values (aka "the lost end-zero")
mhbad=which((mh>(mh1/10))&(mh<(mh2/10)))
mwbad=which((mw>(mw1/10))&(mw<(mw2/10)))
phbad=which((ph>(ph1/10))&(ph<(ph2/10)))

points(mw[mwbad]~mh[mwbad],col="green",pch=19)
points(mw[mhbad]~mh[mhbad],col="green",pch=19)

if (length(mhbad)>0) mh[mhbad]=mh[mhbad]*10 # correction
if (length(mwbad)>0) mw[mwbad]=mw[mwbad]*10 # correction
if (length(phbad)>0) ph[phbad]=ph[phbad]*10 # correction

print(paste(length(mhbad)," 10x lower maternal height corrections",sep=""))
print(paste(length(mwbad)," 10x lower maternal weight corrections",sep=""))
print(paste(length(phbad)," 10x lower paternal height corrections",sep=""))

rm(mhbad,mwbad,phbad)

# deal with the 100-points lower values for length (aka "the lost first(digit)-ONE")

mhbad=which((mh>(mh1-100))&(mh<(mh2-100)))
phbad=which((ph>(ph1-100))&(ph<(ph2-100)))

points(mw[mhbad]~mh[mhbad],col="blue",pch=19)

if (length(mhbad)>0) mh[mhbad]=mh[mhbad]+100 # correction
if (length(phbad)>0) ph[phbad]=ph[phbad]+100 # correction

print(paste(length(mhbad)," 100cm lower maternal height corrections",sep=""))
print(paste(length(phbad)," 100cm lower paternal height corrections",sep=""))

