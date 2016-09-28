

################################################
#  this script fills missing WEIGHT and LENGTH
#  values based on same values at different time
#  points. it imputes all present values and 
#  compares them to actual values to predict
#  incorrect data due to, e.g., OCR problems.
################################################


################################################
#  1  ####################   READ AND CLEAN DATA
################################################

DATADIR="C:/Users/lab/Desktop/MoBa/Caffeine_epi/"
DATAFIL="MFR_350_PDB540_v6_singleton_Q1_newFFQb_Caff_EnergyOk_confounder_imputation_n83774_forcorrectingGAandBW.csv"
csv=read.table(paste(DATADIR,"/",DATAFIL,sep=""),sep=";",header=T,dec=",")
  #colnames(csv); dim(csv)

# fix misunderstanding in the erroneously read headline (might be a temporary problem)
colnames(csv)[1]="PREG_ID_540"

# MATERNAL VARIABLES ARE (usually) CAPITALIZED
colnames(csv)[colnames(csv)=="MORS_ALDER"]="AGE"
colnames(csv)[colnames(csv)=="AA85"] = "WEIGHT"
colnames(csv)[colnames(csv)=="AA86"] = "WEIGHT_1stT"
colnames(csv)[colnames(csv)=="AA87"] = "HEIGHT"
colnames(csv)[colnames(csv)=="PARITY_KOMB"]="PARITY"
colnames(csv)[colnames(csv)=="KJONN"]="SEX" # lets pretend this is "maternal variable"

# CHILD VARIABLES ARE IN LOWER CASE
colnames(csv)[colnames(csv)=="DD212"]="weight6w"
colnames(csv)[colnames(csv)=="DD213"]="length6w"
colnames(csv)[colnames(csv)=="DD214"]="headcr6w"
colnames(csv)[colnames(csv)=="DD218"]="weight3m"
colnames(csv)[colnames(csv)=="DD219"]="length3m"
colnames(csv)[colnames(csv)=="DD220"]="headcr3m"
colnames(csv)[colnames(csv)=="DD224"]="weight5_6m"
colnames(csv)[colnames(csv)=="DD225"]="lenght5_6m"
colnames(csv)[colnames(csv)=="DD226"]="headcr5_6m"
colnames(csv)[colnames(csv)=="EE386"]="weight8m"
colnames(csv)[colnames(csv)=="EE387"]="length8m"
colnames(csv)[colnames(csv)=="EE388"]="headcr8m"
colnames(csv)[colnames(csv)=="EE392"]="weight1y"
colnames(csv)[colnames(csv)=="EE393"]="length1y"
colnames(csv)[colnames(csv)=="EE394"]="headcr1y"
colnames(csv)[colnames(csv)=="EE398"]="weight15_18m"
colnames(csv)[colnames(csv)=="EE399"]="length15_18m"
colnames(csv)[colnames(csv)=="GG15"]="height15_18mQ6"
colnames(csv)[colnames(csv)=="GG16"]="weight15_18mQ6"
colnames(csv)[colnames(csv)=="GG20"]="height2y"
colnames(csv)[colnames(csv)=="GG21"]="weight2y"
colnames(csv)[colnames(csv)=="GG25"]="height3y"
colnames(csv)[colnames(csv)=="GG26"]="weight3y"
      #  ... to be continued  ( I have a feeling that some time-points were not exported from SPSS)
      # colnames(csv) # doublecheck

# empty vector to collect gestational ages
GA=rep(NA,dim(csv)[1])

# priority goes to IVF-dating
GA[!is.na(csv$SVLEN_IVF_DG)]=csv$SVLEN_IVF_DG[!is.na(csv$SVLEN_IVF_DG)]
        # sum(!is.na(GA)) # doublecheck

# the rest gets what-is-already-there values
GA[is.na(GA)]=csv$SVLEN_DG[is.na(GA)]
# sum(is.na(GA)) # doublecheck

# attach GA to main file
csv=data.frame(csv,GA)  # do not run it twice!!! (it will attach the GA vector again)
   #head(csv); rm(GA)

# find those with increadible gestational age due to missing dating of pregnancy
excl=which( (csv$GA>=(43*7)&(is.na(csv$SVLEN_IVF_DG))&(is.na(csv$SVLEN_UL_DG)))   )   # change here 43/42***
    # csv[ids,c("SVLEN_DG","SVLEN_UL_DG","SVLEN_IVF_DG")] # preview

# ... and define cleaned dataset
a1=csv[(!seq(dim(csv)[1]) %in% excl),]
    # dim(a1)

# include those with only credible gestational age
incl=which((a1$GA>=154)&(a1$GA<=307))
    # length(incl)
    # sum(is.na(a1[incl,"GA"]))

# ... and create a cleaned new dataset
a2=a1[incl,]
    # head(a2)
    # dim(a1)[1]-dim(a2)[1]  # difference:  


################################################
#  2  ####################   CONVENIENT DATA
################################################

TO=dim(a2)[1]
n=TO # to run smaller test-sets
#n=1000

w1=a2$VEKT[1:TO]
l1=a2$LENGDE[1:TO]
w2=a2$weight6w[1:TO]
l2=a2$length6w[1:TO]
w3=a2$weight3m[1:TO]
l3=a2$length3m[1:TO]
w4=a2$weight5_6m[1:TO]
l4=a2$lenght5_6m[1:TO]
w5=a2$weight8m[1:TO]
l5=a2$length8m[1:TO]
w6=a2$weight1y[1:TO]
l6=a2$length1y[1:TO]
w7=a2$weight15_18m[1:TO]
l7=a2$length15_18m[1:TO]
#a2$height15_18mQ6
#a2$weight15_18mQ6
   # there are more...




################################################
#  3  ########  GENERATE SIMPLE TIME-BASED FLAGS
################################################
   #(does not take into accopunt that diff 
   #time points have diff scale)

# empty data frame to collect data
  REZ=NULL

# for each inter-time step do:
  for( step in 1:6) { #      adjust based on total time-points (minus one) ***

# collector of flags
  rez=matrix(NA,n,ncol=4)
  colnames(rez)=c(paste("fld",step,sep=""),paste("flb",step,sep=""),
                paste("flw",step,sep=""),paste("fll",step,sep=""))

# input data
  wbfr=get(paste("w",step,sep=""))  # weight before (=bfr)
  wftr=get(paste("w",step+1,sep="")) # weight after (=ftr)
  lbfr=get(paste("l",step,sep="")) # length before
  lftr=get(paste("l",step+1,sep="")) # length after

# distance of two points on w/l coordinates (not adjusted to strech)
  dtemp= sqrt((wftr-wbfr)^2 + (lftr-lbfr)^2 ) # distance
  mndtemp=mean(dtemp,na.rm=T)
  thrdtemp=5*sd(dtemp,na.rm=T)   # arbitrary threshold in SDs  ***

# direction (beta) of two points on w/l coordinates (not adjusted to strech)
  btemp=(wftr-wbfr)/(lftr-lbfr) # (non-robust to quadrants, i.e. 3 quadrant = 1 quadrant)
  btemp[btemp==Inf]=NA # get rid of division by zero problem
  btemp[btemp==-Inf]=NA #  ---- ||--------
  mnbtemp=mean(btemp,na.rm=T)
  thrbtemp=5*sd(btemp,na.rm=T)  # arbitrary threshold in SDs  ***

# for each individual
  rez[,1] = abs(dtemp-mndtemp)>thrdtemp # flag for deviant distance change
  rez[,2] = abs(btemp-mnbtemp)>thrbtemp # flag for deviant direction change
  rez[,3] = wftr<wbfr # flag for weight loss
  rez[,4] = lftr<lbfr # flag for length loss

# rename reasonably (f=flag, d=distance, b=beta, w=weightloss, l=lengthloss,nmbrs=inter-time-jump)
  colnames(rez)=c(paste("fd",step,step+1,sep=""),paste("fb",step,step+1,sep=""),
                paste("fw",step,step+1,sep=""),paste("fl",step,step+1,sep=""))
  #head(rez)
  REZ=cbind(REZ,rez)
  rm(dtemp,btemp,thrdtemp,thrbtemp,mndtemp,mnbtemp,rez)
  
  }  # end of inter-time-steps

  head(REZ)
  dim(REZ)

################################################
#  3  ########  EXPLORE SIMPLE TIME-BASED FALGS
################################################


# option a) based on the total number of flags
  SUM=function(x){sum(x,na.rm=T)}
  col=apply(REZ,1,SUM)
  table(col)
  ids=which(col==4)

# option b) based on clusters in the plot
  ids=which((l7<60)&(l7>40))

# some settings for plotting
  library("RColorBrewer")
  palette <- colorRampPalette(c("white",brewer.pal(5,"Greens")), space = "Lab")
  ymax=15000 # for visual purposes (outliers will be max'ed)
  xmax=100 # ---||---

# draw all suspicious data points at diff times based (on questionnaire #)
  for (id in ids){
    JPGOUTDIR="C:/Users/lab/Desktop/GrowthQC"
    #jpeg(paste(JPGOUTDIR,"/Anomaly_id_",id,".jpeg",sep=""),width = 1420, height = 795) # mute if necessaary
    par(mfrow=c(2,4)) # each time refreshes
  
  # for each time point:
  for( plot in 1:7){
    weight=get(paste("w",plot,sep=""))
    length=get(paste("l",plot,sep=""))

    smoothScatter(length,weight,colramp=palette,pch=1,cex=0.7,nbin=100,nrpoints=150,
              main=paste("timepoint ",plot,sep=""),xlim=c(0,xmax),ylim=c(0,ymax))

    # deal with points outside the plots
    xx=ifelse(length[id]<xmax,length[id],xmax) # coordinate X
    yy=ifelse(weight[id]<ymax,weight[id],ymax) # coordinate Y
    zz=ifelse( (weight[id]>=ymax)|(length[id]>xmax),4,2) # size parameter
    points(xx,yy,col="red",pch=20,cex=zz)
    rm(xx,yy)
  }
#readkey<-function()
cat ("Press [enter] to continue")
line <- readline()
  
    #dev.off() # mute if necessary
  }

###################################################################################


#####################################################
#  5  ############## MANUALY CLEANUP OF NONSENSE DATA
#####################################################
        # it is clear from above that W/L clusters are OCR artifacts. filter them out:

  plot=1
  par(mfrow=c(1,1))
  weight=get(paste("w",plot,sep=""))
  length=get(paste("l",plot,sep=""))
  smoothScatter(length,weight,colramp=palette,pch=1,cex=0.7,nbin=100,nrpoints=150,
                  main=paste("timepoint ",plot,sep=""))

# LEFT SIDE OF POPULATION
      #experiment with exclusion lines:
      #abline(-1000,65)
      #abline(-14300,400)
      #abline(-5800,200)

    idsl=which((weight>(length*65-1000))&(weight>(length*400-14300))&(weight>(length*200-5800))  )
    points(length[idsl],weight[idsl],col="red",pch=20,cex=1.5)

# RIGHT SIDE OF POPULATION
      #experiment with exclusion lines:

    #abline(-19200,400)
    #abline(-3100,100)
    #abline(-7800,200)

    idsr=which((weight<(length*400-19200))|(weight<(length*100-3100))|(weight<(length*200-7800)) )
    points(length[idsr],weight[idsr],col="red",pch=20,cex=1.5)

#  ADD ALL POINTS (no smoothscatter)
points(weight~length) # black
points(length[idsl],weight[idsl],col="red",pch=20,cex=1.5) # left
points(length[idsr],weight[idsr],col="red",pch=20,cex=1.5) # right



#####################################################
#  INVESTIGATE, WHERE CLUSTER MEMBERS END UP IN...
#####################################################

# find the ids of the cluster members:
plot=1
weight=get(paste("w",plot,sep=""))
length=get(paste("l",plot,sep=""))
ids=which((weight>2000)&(weight>(length*400-12300))&(length>=30))

par(mfrow=c(2,3))
for (plot in 1:6){
weight=get(paste("w",plot,sep=""))
length=get(paste("l",plot,sep=""))
#smoothScatter(length,weight,colramp=palette,pch=1,cex=0.7,nbin=100,nrpoints=150,
#              main=paste("timepoint ",plot,sep=""),xlim=c(0,100),ylim=c(0,15000))
plot(length,weight,main=paste("timepoint ",plot,sep=""),xlim=c(0,100),ylim=c(0,15000))
points(length[ids],weight[ids],col="red",pch=20,cex=1.5)
}

... continue with other time points...


######################################################
#  6  ############## STATISTICAL PROOF OF OCR FAILURES
######################################################

weight=get(paste("w",1,sep=""))
length=get(paste("l",1,sep=""))
# inspect the 30-39 interval
par(mfrow=c(1,1))
 plot(weight~length) # black
idsu=which((weight>(length*65-1000)&(length<=39)&(length>=30)) ) #upper
points(length[idsu],weight[idsu],col="green",pch=20,cex=1.5)
idsl=which((weight<=(length*65-1000)&(length<=39)&(length>=30)) ) #lower
points(length[idsl],weight[idsl],col="red",pch=20,cex=1.5)
   #abline(v=50)


lup=length[which(weight>(length*65-1000))]
llw=length[which(weight<=(length*65-1000))]
lup=as.character(lup[which((lup>=30)&(lup<=39)&(!is.na(lup)))])
llw=as.character(llw[which((llw>=30)&(llw<=39)&(!is.na(llw)))])
par(mfrow=c(1,2))
hist(as.numeric(substr(lup,2,2)),breaks=10,col="grey",main="outliers (green)",xlab="length: 30+x")
hist(as.numeric(substr(llw,2,2)),breaks=10,col="grey",main="normal (red)",xlab="length: 30+x")
dat=data.frame(up=as.numeric(table(lup)),lw=as.numeric(table(llw)))
chisq.test(dat)

# compare with other intervals:
reg40=as.character(length[which((length>=40)&(length<=49))])
reg50=as.character(length[which((length>=50)&(length<=59))])

par(mfrow=c(2,2))
hist(as.numeric(substr(lup,2,2)),breaks=10,col="grey",main="outliers 30-39",xlab="length: 30+x")
hist(as.numeric(substr(llw,2,2)),breaks=10,col="grey",main="normal 30-39",xlab="length: 30+x")
hist(as.numeric(substr(reg40,2,2)),breaks=10,col="grey",main="normal 40-49",xlab="length: 40+x")
hist(as.numeric(substr(reg50,2,2)),breaks=10,col="grey",main="normal 50-59",xlab="length: 50+x")

# proove that 40-49 and 50-59 intervals had separate effects (more detailed proof)
lupup=length[which( (weight> (length*65-1000))&(weight> 3100) )] # upper-upper
luplw=length[which( (weight<=(length*65-1000))&(weight<=3100) )] # upper-lower
lupup=as.character(lupup[which((lupup>=30)&(lupup<=39)&(!is.na(lupup)))])
luplw=as.character(luplw[which((luplw>=30)&(luplw<=39)&(!is.na(luplw)))])

# collored plot
library("RColorBrewer")
orng=brewer.pal(5,"Oranges")[c(2,5)]
gree=brewer.pal(5,"Greens")[c(2,5)]
par(mfrow=c(1,1))
plot(weight~length) # black
idsu=which((weight>(length*65-1000)&(length<=39)&(length>=30)&(weight>3100)) ) #upper
points(length[idsu],weight[idsu],col=orng[1],pch=20,cex=1.5)
idsl=which((weight>(length*65-1000)&(length<=39)&(length>=30)&(weight<=3100)) ) #lower
points(length[idsl],weight[idsl],col=gree[1],pch=20,cex=1.5)
# color the dots that are >39 in length
frts=which( (weight<=3100)&(weight>1000)&(length>=40)&(length<=49)) #upper
points(length[frts],weight[frts],col=gree[2],pch=20,cex=1.5)
ffts=which( (weight>3100)&(weight<5200)&(length>=50)&(length<=59)) #lower
points(length[ffts],weight[ffts],col=orng[2],pch=20,cex=1.5)

# second digit distribution histograms
par(mfrow=c(2,2))
hist(as.numeric(substr(lupup,2,2)),breaks=10,col="grey",main="outliers (upper) 30-39",xlab="length: 30+x")
hist(as.numeric(substr(luplw,2,2)),breaks=10,col="grey",main="outliers (lower) 30-39",xlab="length: 30+x")
hist(as.numeric(substr(reg50,2,2)),breaks=10,col="grey",main="normal 50-59",xlab="length: 50+x")
hist(as.numeric(substr(reg40,2,2)),breaks=10,col="grey",main="normal 40-49",xlab="length: 40+x")





#####################################################
#  4  ########  NONPARAMETRIC ROBUST TIME-BASED FALGS
#####################################################
    # (a tool to impute missing value based on other time points)


## CREATE REFERENCE TABLES BASED ON NONMISSING, CLEANED VALUES

lines=13  # MUST BE ODD
timepoints=c(1,2,3)
references=matrix(NA,nrow=dim(a2)[1],ncol=length(timepoints)*lines*2)
dim(references)

for( time in timepoints) {
  
  # read-in
  weight=get(paste("w",time,sep=""))
  length=get(paste("l",time,sep=""))
  
  # clean-up
  validids=which((!is.na(weight))&(!is.na(length)))
  weight=weight[validids] # BOTH VALUES MUST BE PRESENT
  length=length[validids] # BOTH VALUES MUST BE PRESENT
  
  # outlier clean-up step is stil missing
  #...                                 
  #...

  # for timepoint 1
  if (time==1) {
#abline(-8500,300)  # no above
#abline(-1300,70)   # no above
#abline(-14000,300) # no below
#abline(-1800,60)   # no below
bad=which( ((weight>length*300-8500)&(weight>length*70-1300)) | ((weight<length*300-14000)|(weight<length*60-1800)) )
good=which(! ((weight>length*300-8500)&(weight>length*60-1000)) | ((weight<length*300-14000)|(weight<length*60-1800)) )
plot(weight~length)
points(weight[bad]~length[bad],col="red",pch=19)
  }

if (time==2) {
#  abline(-18000,500)  # no above
#  abline(-3300,150)   # no above
#  abline(-30000,500) # no below
#  abline(-6000,150)   # no below
  bad=which( ((weight>length*500-18000)&(weight>length*150-3300)) | ((weight<length*500-30000)|(weight<length*150-6000)) )
  good=which(! ((weight>length*500-18000)&(weight>length*150-3300)) | ((weight<length*500-30000)|(weight<length*150-6000)) )
  plot(weight~length)
  points(weight[bad]~length[bad],col="red",pch=19)
}

if (time==3) {
#  abline(-19000,500)  # no above
#  abline(-3300,150)   # no above
#  abline(-30000,500) # no below
#  abline(-6000,150)   # no below
  bad=which( ((weight>length*500-19000)&(weight>length*150-3300)) | ((weight<length*500-30000)|(weight<length*150-6000)) )
  good=which(! ((weight>length*500-19000)&(weight>length*150-3300)) | ((weight<length*500-30000)|(weight<length*150-6000)) )
  plot(weight~length)
  points(weight[bad]~length[bad],col="red",pch=19)
}

if (time>3) {
  print("THE OUTLIERS WERE NOT MANUALLY DEFINED FOR THIS TIMEPOINT YET")
}
  
    # midpoint
  wmd=median(weight[good],na.rm=T) 
  lmd=median(length[good],na.rm=T)
    
  # CENTER TO THE ZERO
  weight=weight[good]-wmd
  length=length[good]-lmd
  smoothScatter(length,weight,colramp=palette,pch=1,cex=0.7,nbin=100,nrpoints=150,
                main=paste("timepoint ",time,sep=""))
  
  # scale
  weight=weight/sd(weight,na.rm=T)
  length=length/sd(length,na.rm=T)
  smoothScatter(length,weight,colramp=palette,pch=1,cex=0.7,nbin=100,nrpoints=150,
                main=paste("timepoint ",time,sep=""))
  
  range(weight)
  range(length)
  
  # JITTER ALL VALUES SO THAT THERE WOULD NOT BE ZERO in x and y
  weight=weight+rnorm(length(weight),0,(sd(weight,na.rm=T)*0.0001)) # noise with sd = 50ug
  length=length+rnorm(length(length),0,(sd(length,na.rm=T)*0.0001)) # noise with sd = 2um
    smoothScatter(length,weight,colramp=palette,pch=1,cex=0.7,nbin=100,nrpoints=150,
                main=paste("timepoint ",time,sep=""))
  
  
  # make sure that there are no values that are ==0
  sum(length==0,na.rm=T)
  sum(weight==0,na.rm=T)
  
    # define a quadrants
  quad=rep(NA,length(weight))
  quad[ which((weight>0)&(length>0))] = 1
  quad[ which((weight>0)&(length<0))] = 2
  quad[ which((weight<0)&(length<0))] = 3
  quad[ which((weight<0)&(length>0))] = 4
  table(quad,useNA="a") # THINK OF SOLUTION TO GET RID OF NAs
  
  abline(v=0)
  abline(h=0)
  
  # give each individual a value in degrees:
  ideg=atan(weight/length)/(2*pi)*360 # not processed yet (belongs to only one quadrant)
  
  # correct degrees based on a quadrant
  ideg[which(quad==1)]=ideg[which(quad==1)]
  ideg[which(quad==2)]=ideg[which(quad==2)]+180
  ideg[which(quad==3)]=ideg[which(quad==3)]+180
  ideg[which(quad==4)]=90+ideg[which(quad==4)]+270

  # check ideg content
  sum(is.na(ideg))
  hist(ideg,breaks=100,col="grey")
  range(ideg)  # strange shape...***
  
  # create angular segments
  qntls=quantile(ideg,probs=seq(0,1,0.02)) # assuming 50 output segments
  lbs=seq(50)
  segcol=cut(ideg,qntls,labels=lbs)
  table(cut(ideg,qntls))
  
  
  # visually check whether all segments ar ein order
  for(seg in 1:50) {
  smoothScatter(length,weight,colramp=palette,pch=1,cex=0.7,nbin=100,nrpoints=10,
                main=paste("timepoint ",time,sep=""))
  points(length[segcol==seg],weight[segcol==seg],pch=20)
  }
  
  ### to universalize plots at different time points:
  # 1) find the middle degree for each segment
  middeg=NULL
  for(seg in 1:50) middeg=c(middeg,median(ideg[segcol==seg],na.rm=T))
  
  # 2) find population growth direction (beta -> degrees)
  popdir=atan(coef(summary(lm(weight~length)))[2,1])/(2*pi)*360  # assuming beta is positive (biology)
  abline(0,coef(summary(lm(weight~length)))[2,1]) # make sure outliers do not affect beta too much
  
  # 3) which segment faces same direction as growth curve
  segid=which(abs(middeg-popdir)==min(abs(middeg-popdir)))
  
  # 4) transpoze all degrees so that population faces the zero degree
    ideg=ideg-middeg[segid]
    range(ideg)  
    ideg[ideg<0]= 360+ideg[ideg<0] # quite tricky logic here
    hist(ideg,breaks=100,col="grey")
    range(ideg)    
  
  
  # RENAME angular segments
  qntls=quantile(ideg,probs=seq(0,1,0.02)) # assuming 50 output segments
  lbs=seq(50)
  #table(cut(ideg,qntls,include.lowest = TRUE))
  segcol=cut(ideg,qntls,labels=lbs,include.lowest = TRUE) # as in segment column

  
  # DISTANCES FROM TEH CENTER
  DISTZ=sqrt(weight^2 + length^2) # Z=zero
  hist(DISTZ,breaks=100,col="grey")
  
  ## quite poor predictor of outliers (does not take into account shape)
  smoothScatter(length,weight,colramp=palette,pch=1,cex=0.7,nbin=100,nrpoints=150,
                main=paste("timepoint ",time,sep=""))
  points(weight[DISTZ>5]~length[DISTZ>5],pch=20) # pur evil
  
  
    # ESTIMATE PERCENTILES OF DISTANCE IN EACH SEGMENT
    table(segcol,useNA="a")  # shoudl be no NA's
  
  ### VISUALIZE THE CONCEPT
  qnt90=NULL
  for(seg in 1:50){
    extr=DISTZ[which(segcol==seg)]
    extr=extr[order(extr)]
    qnt90=c(qnt90,extr[floor(length(extr)*0.90)])
  }
  
  # visually check validity
  smoothScatter(length,weight,colramp=palette,pch=1,cex=0.7,nbin=100,nrpoints=10,
                main=paste("timepoint ",time,sep=""))
  abline(v=0); abline(h=0)
  #plot(length,weight)
      for(seg in 1:50) {
      points(length[which((DISTZ>qnt90[seg])&(segcol==seg))],
      weight[which((DISTZ>qnt90[seg])&(segcol==seg))],pch=20,col="red")
    }
  
  
  
  #sandbox
  # ARTEFACT AT VERTICAL LINE::::::: ??
  sum(segcol==seg,na.rm=T)  # 3 why??
  table(segcol)
  
  
  # go back and exclude outliers
  # go back and make sure that direction of growth is correct visually
  # create REF MATRIX (dimensions(rows)= ceiling(dim(a2)[1]/50) =1668
  # with ordered distances for each segment, for each timepoint
  # CREATE TWO MATRIXES: (dimensions= dim(a2)[1] = 83355)
  # one with DISTZ for each individual at each timepoint + column with segment id
  # one with DEGREES for each individual at each timepoint  + column with segment id
  # find out how to evaluate the change in degrees at junction 0-360
  
  
  
  
  
  
  
  
  # areas of plot
  lines=lines  # MUST BE ODD
  segments=lines*2    
  strech=30000/60  # strech factor (visual)
  tns=tan(((2*pi)/(lines*2))*(seq(lines)))
  for(d in 1:lines)  abline(0,1*tns[d]*strech)
    
  for (d in 1:segments) {
    if (d==1){  # segment No 1 is always on bottom (6 oclock)
    marked = which( (weight<(length*tns[floor(lines/2)]*strech)) & (weight<(length*tns[ceiling(lines/2)]*strech)) )
      } else {
      if (d==(segments/2+1)) {
      marked = which( (weight>=(length*tns[floor(lines/2)]*strech)) & (weight>=(length*tns[ceiling(lines/2)]*strech))  )
    } else {
      if (d<=segments/2) {
        lineids=c(NA,ceiling(lines/2):lines,1:(floor(lines/2)))
        #length(lineids)
        marked = which( (weight>=(length*tns[lineids[d]]*strech)) & (weight<(length*tns[lineids[d+1]]*strech)) )
      #length(marked)
        } else {
        if (d>(segments/2+1)) {
          lineids=c(rep(NA,segments/2+1),ceiling(lines/2):lines,1:(floor(lines/2)))
          marked = which( (weight<(length*tns[lineids[d]]*strech)) & (weight>=(length*tns[lineids[d+1]]*strech)) )
        } else {
          print("problem")
        } 
      }
     }
    }    
      # save the reference of distances from center for each segment
  distz=DISTZ[marked]
  distz=distz[order(distz)]
  addnas=dim(a2)[1]-length(distz)
  references[,((time-1)*segments)+d]=c(distz,rep(NA,addnas))
  rm(distz)
  
    smoothScatter(length,weight,colramp=palette,pch=1,cex=0.7,nbin=100,nrpoints=150,
                main=paste("timepoint ",time,sep=""),xlim=c(-30,30),ylim=c(-15000,15000))
  points(weight[marked]~length[marked])
    rm(marked)
  }
  
  # next time estimate distance for each direction vector
    
}

head(references)
dim(references)
references[1:200,66]
#












# second pair of timepoints
for ( bd in 1:nbads ){
  idss=ids[bd]
  par(mfrow=c(1,2))
  ymax=max(c(a2$weight6w,a2$weight3m),na.rm=T)
  xmax=max(c(a2$length6w,a2$length3m),na.rm=T)
  plot(a2$weight6w~a2$length6w,xlim=c(0,xmax),ylim=c(0,ymax))
  points(a2$length6w[idss],a2$weight6w[idss],col="red",pch=20)
  plot(a2$weight3m~a2$length3m,xlim=c(0,xmax),ylim=c(0,ymax))
  points(a2$length3m[idss],a2$weight3m[idss],col="red",pch=20)
  readkey <- function()
    cat ("Press [enter] to continue")
  line <- readline()
}





display.brewer.all()
par(mfrow=c(1,1))

Lab.palette <- colorRampPalette(c("white",brewer.pal(5,"Greens")), space = "Lab")
smoothScatter(a2$LENGDE,a2$VEKT, colramp = Lab.palette,pch=1,cex=0.7)

## WHY DID NOT WE EXPORT VALUES FOR HEIGHT AT 2+ years?













###########################################################
############   CHECK WHETHER CHILD VARIABLES ARE REASONABLE
###########################################################

colnames(a2)
hist(a2$LENGDE,breaks=100,col="grey")
hist(a2$HODE,breaks=100,col="grey")


plot(a2$VEKT~a2$weight6w); abline(0,1)
plot(a2$VEKT~a2$weight3m); abline(0,1)
plot(a2$VEKT~a2$weight5_6m); abline(0,1)

plot(a2$VEKT~a2$LENGDE)

library("RColorBrewer")
display.brewer.all()




### COLOR BASED ON GESTATIONAL AGE
groups=6
step=round((range(a2$GA)[2]-range(a2$GA)[1])/groups)
MIN=min(a2$GA)
MAX=max(a2$GA)
THRS=NULL
for (i in 1:(groups)) THRS=c(THRS,(MIN+step*i))
CAT=rep(NA,dim(a2)[1])
CAT[ (a2$GA>=MIN)&(a2$GA<THRS[1])]=1
CAT[ (a2$GA>=THRS[1])&(a2$GA<THRS[2])]=2
CAT[ (a2$GA>=THRS[2])&(a2$GA<THRS[3])]=3
CAT[ (a2$GA>=THRS[3])&(a2$GA<THRS[4])]=4
CAT[ (a2$GA>=THRS[4])&(a2$GA<THRS[5])]=5
CAT[ (a2$GA>=THRS[5])&(a2$GA<=MAX)]=6
table(CAT)
display.brewer.pal(8,"RdYlGn")
cols=brewer.pal(8,"RdYlGn")[c(1,2,3,6,7,8)]
sub=a2[rev(order(a2$SVLEN_DG)),]  
par(mfrow=c(1,3))
plot(a2$VEKT~a2$LENGDE,col=cols[CAT],pch=20,xlab="Birthlength",ylab="Birthweight"); abline(v=40)
plot(a2$weight6w~a2$length6w,col=cols[CAT],pch=20,xlab="length6w",ylab="weight6w"); abline(v=40)
plot(a2$weight3m~a2$length3m,col=cols[CAT],pch=20,xlab="length3m",ylab="weight3m"); abline(v=40)



abline(v=20)
abline(v=65)

# some strange stufff...
sum(((is.na(a2$VEKT))|(is.na(a2$LENGDE)))==1)
sub=a2[ ((a2$LENGDE<10)|(a2$LENGDE>64))&(!is.na(a2$VEKT))&(!is.na(a2$LENGDE)) ,]
dim(sub)
sub[,c("SVLEN_DG","SVLEN_UL_DG","SVLEN_SM_DG","SVLEN_IVF_DG","VEKT","LENGDE","weight6w","length6w")]


lengthdif6w=(a2$length6w-a2$LENGDE)
hist(lengthdif6w,breaks=100,col="grey",xlim=c(-5,20))
sum((lengthdif6w<0)|(lengthdif6w>15),na.rm=T)



#########   COLOR BY MATERNAL HEIGHT  #########
library("RColorBrewer")
display.brewer.all()
display.brewer.pal(8,"RdYlGn")
cols=brewer.pal(8,"RdYlGn")[c(1,2,3,6,7,8)]
groups=6
step=round((range(a2$HEIGHT,na.rm=T)[2]-100) /groups)
MIN=min(a2$HEIGHT,na.rm=T)
MAX=max(a2$HEIGHT,na.rm=T)
THRS=NULL
for (i in 1:(groups)) THRS=c(THRS,(MIN+step*i))
CAT=rep(NA,dim(a2)[1])
CAT[ (a2$HEIGHT>=MIN)&(a2$HEIGHT<THRS[1])]=1
CAT[ (a2$HEIGHT>=THRS[1])&(a2$HEIGHT<THRS[2])]=2
CAT[ (a2$HEIGHT>=THRS[2])&(a2$HEIGHT<THRS[3])]=3
CAT[ (a2$HEIGHT>=THRS[3])&(a2$HEIGHT<THRS[4])]=4
CAT[ (a2$HEIGHT>=THRS[4])&(a2$HEIGHT<THRS[5])]=5
CAT[ (a2$HEIGHT>=THRS[5])&(a2$HEIGHT<=MAX)]=6
sum(is.na(CAT))
sum(is.na(a2$HEIGHT))

sub=data.frame(GA=a2$GA,CW=a2$VEKT,CH=a2$LENGDE,MH=a2$HEIGHT,CAT=CAT)
sub=sub[rev(order(sub$MH)),]
table(sub$CAT)
plot(sub$CW~sub$CH,col=cols[sub$CAT],pch=20)


sub1=sub[sub$GA>260,]
plot(sub1$CW~sub1$CH,col=cols[sub1$CAT],pch=20)


##  COLOR BY COMPLICATIONS
display.brewer.all()
display.brewer.pal(11,"Paired")
cols=brewer.pal(11,"Paired")
cols=c(cols,"black")
#hist(a2$GA[a2$DODKAT==10])
sub=data.frame(GA=a2$GA,CW=a2$VEKT,CH=a2$LENGDE,MH=a2$HEIGHT,CAT=(a2$DODKAT+1))
sub=sub[order(sub$CAT),]
head(sub)
table(sub$CAT)
plot(sub$CW~sub$CH,col=cols[sub$CAT],pch=20)


##  COLOR BY MALFORMATIONS:
table(a2$C00_MALF_ALL)
cols=c("blue","red")
sub=data.frame(GA=a2$GA,CW=a2$VEKT,CH=a2$LENGDE,MH=a2$HEIGHT,CAT=(a2$C00_MALF_ALL+1))
sub=sub[order(sub$CAT),]
head(sub)
table(sub$CAT)
plot(sub$CW~sub$CH,col=cols[sub$CAT],pch=20,xlab="Birthlength",ylab="Birthweight")
abline(v=40)


# color by CLUSTERING ON THE PLOT  (based on regression)
model=loess(a2$VEKT~a2$LENGDE)
CW.exp=predict(model,a2$LENGDE)
dif=a2$VEKT-CW.exp
hist(dif,breaks=100,col="grey")
SD=sd(dif,na.rm=T)
CAT=abs(dif)>4*SD # ***
table(CAT)
CAT=CAT+1
cols=c("green","red")
sub=data.frame(GA=a2$GA,CW=a2$VEKT,CH=a2$LENGDE,MH=a2$HEIGHT,CAT=CAT)
sub=sub[order(sub$CAT),]
table(sub$CAT)
plot(sub$CW~sub$CH,col=cols[sub$CAT],pch=20,xlab="Birthlength",ylab="Birthweight")

# based on clustering (K-means etc)

library(cluster)
dat=a2[ (!is.na(a2$LENGDE))&(!is.na(a2$VEKT)),c("LENGDE","VEKT")]
head(dat)
plot(dat); dim(dat)
badass=which(dat$LENGDE<40)
godass=which(dat$LENGDE>=40)
dat1=dat[c(badass,sample(godass,5000)),]
dat2=dat[badass,]

m1=as.numeric(dat$LENGDE); m1=(m1-mean(m1))/sd(m1)
m2=as.numeric(dat$VEKT); m2=(m2-mean(m2))/sd(m2)
m=matrix(c(m1,m2),ncol=2)
head(m);dim(m)
#CAT=pam(m,3)$clustering

CAT=kmeans(m, 2,iter.max = 10000)$cluster
plot(m,col=CAT,pch=20)



1-pnorm(5)

cage=rep(0.1,dim(a2)[1])
a3=data.frame(a2,cage=cage)
colnames(a3)

##############################
#  based on WHO

plot(a3$VEKT~a3$LENGDE)
a4=a3[sample(seq(dim(a3)[1]),1000,replace=F),]
sum(is.na(a4$LENGDE))
a4=a4[  (!is.na(a4$VEKT))&(!is.na(a4$LENGDE)),]
head(a4)
colnames(a4)

a4$cage=4
  
igrowup.restricted(FilePath="C:/Users/lab/Desktop/WHO", 
                   FileLab="MySurvey",mydf=a4,sex=SEX,age=cage,age.month=F, 
                   weight=weight3y, lenhei=height3y)




rez=read.csv("MySurvey_z_rc.csv",sep=",")
head(rez)
hist(rez$zlen,breaks=20,col="grey")
hist(rez$zbmi,breaks=20,col="grey")
hist(rez$zwei,breaks=20,col="grey")





# example of clustering
head(dat)
x <- rbind(cbind(rnorm(1000,0,1), rnorm(1000,0,1)),
           cbind(rnorm(500,3,0.5), rnorm(500,3,1)),
           cbind(rnorm(200,0,0.1), rnorm(200,3,0.5)),
           cbind(rnorm(100,4,0.5), rnorm(100,0,2)))
class(x)
head(x)
plot(x)
CAT=pam(x, 3)$clustering
plot(x,col=CAT,pch=20)

CAT=kmeans(x, 4)$cluster
plot(x,col=CAT,pch=20)








################################################
#  x  ####################   CELLULAR AGE FACTOR
################################################

### a vector to correct child's age  (childs age in days minus correction = childs corrected age)
agecor= rep(NA,dim(a2)[1])
agecor [a2$GA>=259]= 0
agecor[a2$GA<259]=280 - a2$GA[a2$GA<259]
hist(agecor,breaks=100,col="grey")
#sum(agecor>0)
