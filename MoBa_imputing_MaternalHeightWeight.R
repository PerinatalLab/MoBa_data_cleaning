
############################################################
#######  IMPUTE MATERNAL HEIGHT AND WEIGHT
############################################################

#### This script imputes missing weight (pre-preg and 1st trim)
#### and height values, based on energy intake, age and other things.

# read in the csv data files
## (note the required columns and their order)

## get the date stamp
date_stamp = paste(unlist(strsplit(substr(Sys.time(),1,10),"-")),collapse="")
## get the hash of the current state of the Git folder
setwd("~/Documents/gitrep/MoBa_data_cleaning")
hash = system(paste("git log --pretty=format:'%h' -n 1"),intern=TRUE)

file_dir = "~/Desktop/MoBa_v6/"
file_out = paste(file_dir,"MOBA_PDB1581_IMPUTED_maternalHgh1Wgh1Wgh2_",date_stamp,"_",hash,".txt",sep="")

## this file contains the cleaned and corrected mother height/weight info
q1=read.csv(paste(file_dir,"MOBA_PDB1581_CLEANED_maternalHgh1Wgh1Wgh2_",date_stamp,"_",hash,".txt",sep=""),
            sep="\t", header=T)
names(q1)=c("PREG_ID","WEIGHT","WEIGHT_1stT","HEIGHT","flAA85","flAA86","flAA87")
head(q1); dim(q1)

## this file contains energy info retrieved directly from MoBa without any changes
en=read.csv(paste(file_dir,"output_q2_kilojoules.csv",sep=""), sep=",", header=T)
names(en)=c("PREG_ID","KJ")
head(en); dim(en)

## this file is required only for the mother age info (also taken directly from MoBa)
mfr=read.csv(paste(file_dir,"output_mfr_basicfetalinfo.csv",sep=""),sep=",",header=F)
names(mfr)=c("PREG_ID","CHILDNUM","AGE","GA","SEX","BIRTHWEIGHT","PARITY")
head(mfr); dim(mfr)

## merge the age with mother info
library(sqldf)
M=sqldf("SELECT DISTINCT mfr.AGE, q1.* FROM q1 LEFT JOIN mfr ON mfr.PREG_ID=q1.PREG_ID")
head(M); dim(M)

#####################################################################
# Fistly we checked the weight and height of mothers
######################################################################

# ESTIMATE WHAT IS THE MEDIAN DIFFERENCE BETWEEN PRE and 1st TRIm WEIGHT IN "NORMAL CASES"
# define normal cases
S=M[(!is.na(M$WEIGHT))&(!is.na(M$WEIGHT_1stT))&(M$WEIGHT>=35)&(M$WEIGHT<=200)&(M$WEIGHT_1stT>=35)&(M$WEIGHT_1st<=200),]
S$dif=S$WEIGHT_1stT-S$WEIGHT
hist(S$dif,breaks=100,col="grey")
right=(abs(S$dif))<=15  # arbitrary
S=S[right,]
# what is the median difference between pre and 1st trimester maternal weight
meddif=median(S$dif)  

## attach the weight difference while removing duplicates
en=sqldf("SELECT DISTINCT en.*, WEIGHT_1stT-WEIGHT as DIF FROM en LEFT JOIN S on S.PREG_ID=en.PREG_ID")

## set the upper and lower energy intake caps
en$KJ[en$KJ>25000]=25000
en$KJ[en$KJ<2500]=2500

## create rounded bins for visualizing
en$int=round(en$KJ,-2)

## predict difference for any KJ intake using a linear model
## (note that prediction uses per-individual data, not bin average)
enpredictor=lm(en$DIF~en$KJ)
en$preddif=predict(enpredictor,newdata=data.frame(KJ=en$KJ))
summary(enpredictor)

## plot the mean weight growth for each energy intake bin (circles),
## linear predictor (red) and samples per bin (grey)
enplot=sqldf("SELECT int, COUNT(int) as freq, AVG(DIF) as mean FROM en GROUP BY int")
plot(enplot$int,enplot$mean,xlab="energy intake",ylab="mean weight difference")
abline(enpredictor,col="red")
points(enplot$int,log(enplot$freq,3),type="l",col="grey")
axis(side=4,at=seq(0,10,2),labels=3^seq(0,10,2))


# FIRST ROUND:
# USE THE WEIGHT VALUE FROM 1st TRIMESTER (minus 3kg), WHEN PRE-PREG VALUE IS STRANGE/MISSING
M1=sqldf("SELECT M.*, preddif FROM M LEFT JOIN en ON en.PREG_ID=M.PREG_ID")
M1$preddif[is.na(M1$preddif)]=meddif
head(M1); dim(M1)

## mothers which indicated only 1stTrim weight
group1=which(is.na(M1$WEIGHT) & !is.na(M1$WEIGHT_1stT))
length(group1)
M1[group1,"WEIGHT"]=M1[group1,"WEIGHT_1stT"]-M1[group1,"preddif"]
M1$flAA85[group1]="imputed1"

## mothers which indicated only pre-pregnancy weight
group2=which(is.na(M1$WEIGHT_1stT) & !is.na(M1$WEIGHT))
length(group2)
M1[group2,"WEIGHT_1stT"]=M1[group2,"WEIGHT"]+M1[group2,"preddif"]
M1$flAA86[group2]="imputed2"


# SECOND ROUND:
## M2 will contain BMI calculations, and missing BMI predictions from age
M2=M1

# run the model
## LOESS fitting gives very similar results to linear fitting, except possibly in extreme parts of the range
## residual SE was 4.278 in LOESS and 4.281 in a linear model
M2$BMIobs<-M2$WEIGHT/(M2$HEIGHT/100)^2

BMI.lo<-loess(M2$BMIobs~M2$AGE)
M2$BMIexp<-predict(BMI.lo,M2$AGE)
summary(BMI.lo)

####### M2 now contains BMI predicted and BMI observed columns ##############

# case 1: no height, replace with mean
## (height is a very weak predictor of BMI, hence using the reverse calculation gives clearly wrong results)
noheight<-is.na(M2$HEIGHT)
M2$HEIGHT[noheight]=round(mean(M2$HEIGHT,na.rm=T),0)
M2$flAA87[noheight]="imputed3"
                          
# case 2: no weight, calculate it from mean height and age-based BMI prediction
noweight<-is.na(M2$WEIGHT) & !is.na(M2$HEIGHT)
M2$WEIGHT[noweight]=round(M2$BMIexp[noweight]*(M2$HEIGHT[noweight]/100)^2,1)
M2$flAA85[noweight]="imputed3"

## also, calculate weight at 1st trim. based on newly-predicted pre-preg weight
noweight<-is.na(M2$WEIGHT_1stT) & !is.na(M2$WEIGHT)
M2[noweight,"WEIGHT_1stT"]=M2[noweight,"WEIGHT"]+M2[noweight,"preddif"]
M2$flAA86[noweight]="imputed3"

## since height was predicted just from the mean, there should be no NAs left by now. check:
sum(is.na(M2$WEIGHT)); sum(is.na(M2$WEIGHT_1stT)); sum(is.na(M2$HEIGHT))

## review results
hist(M2$WEIGHT, breaks=50, col="grey")
hist(M2$HEIGHT, breaks=50, col="grey")

## in these tables, numbers refer to changes done by previous script:
## handwriting and OCR corrections, multiparity corrections, multiparity imputations
## 0 - original; 1 - reliable changes; 2 - imputed; 3 - deleted
table(M2$flAA85)
table(M2$flAA86)
table(M2$flAA87)

## write maternal data output
write.table(M2[,c("PREG_ID","WEIGHT","WEIGHT_1stT","HEIGHT","flAA85","flAA86","flAA87")],
            file=file_out,sep="\t",col.names=T,row.names=F,quote=F)