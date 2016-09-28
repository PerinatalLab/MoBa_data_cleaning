library("RColorBrewer")
palette <- colorRampPalette(c("white",brewer.pal(5,"Greens")), space = "Lab")

######################################################################################
# defines distribution of weight-gain and length-gain in different between-time-points 
######################################################################################



################################################
#  1  ####################   READ AND CLEAN DATA
################################################

DATADIR="C:/Users/lab/Desktop/MoBa/Caffeine_epi/"
DATAFIL="MFR_350_PDB540_v6_singleton_Q1_newFFQb_Caff_EnergyOk_confounder_imputation_n83774_forcorrectingGAandBW.csv"
csv=read.table(paste(DATADIR,"/",DATAFIL,sep=""),sep=";",header=T,dec=",")
#colnames(csv); dim(csv)

# RENAME COLUMNS TO CONVENIENT CODING
colnames(csv)[1]="PREG_ID_540"
colnames(csv)[colnames(csv)=="VEKT"]="w1"
colnames(csv)[colnames(csv)=="LENGDE"]="l1"
colnames(csv)[colnames(csv)=="DD212"]="w2" #"weight6w"
colnames(csv)[colnames(csv)=="DD213"]="l2" #"length6w"
colnames(csv)[colnames(csv)=="DD218"]="w3" #"weight3m"
colnames(csv)[colnames(csv)=="DD219"]="l3" #"length3m"

# empty vector to collect gestational ages
GA=rep(NA,dim(csv)[1])
# priority goes to IVF-dating
GA[!is.na(csv$SVLEN_IVF_DG)]=csv$SVLEN_IVF_DG[!is.na(csv$SVLEN_IVF_DG)]
# the rest gets what-is-already-there values
GA[is.na(GA)]=csv$SVLEN_DG[is.na(GA)]
# attach GA to main file
csv=data.frame(csv,GA)  # do not run it twice!!! (it will attach the GA vector again)
# find those with increadible gestational age due to missing dating of pregnancy
excl=which( (csv$GA>=(43*7)&(is.na(csv$SVLEN_IVF_DG))&(is.na(csv$SVLEN_UL_DG)))   )   # change here 43/42***
# ... and define cleaned dataset
a1=csv[(! seq(dim(csv)[1]) %in% excl),]
# include those with only credible gestational age
incl=which((a1$GA>=154)&(a1$GA<=307))
# ... and create a cleaned new dataset
## NOTE: it is very important to retain this structure: ID, GA, w1, l1, w2, l2.... due to automatic identification of columns further down
a2=a1[incl,c("PREG_ID_540","GA","w1","l1","w2","l2","w3","l3")] ##  next time expand to more timepoints
#head(a2); colnames(a2)



###############################################
#########    DELET OUTLIERS
###############################################

clean=a2 # a2 will remain WITH outliers and "clean" will be cleaned
timepoints=c(1,2,3) # manually define how many timepoints are ready to be analysed 

POINTS=matrix(NA,nrow=length(timepoints),ncol=4) # COLLECT THE VALUES OF OUTLIER-FREE DADATSET
colnames(POINTS)=c("wght_med","wght_sd","lngth_med","lngth_sd")
rownames(POINTS)=paste("time_",timepoints,sep="")

par(mfrow=c(1,3))
for( time in timepoints) {

  # read-in
  weight=clean[, (2+1+(time-1)*2) ] ## AUTOMATIC ASSUMPTION THAT w and l are ordered in a particular way! ***
  length=clean[, (2+2+(time-1)*2) ]
  
  ###  identify outliers using a manually predifined mask
  if (time==1) {
    #abline(-8500,300)  # no above
    #abline(-1300,70)   # no above
    #abline(-14000,300) # no below
    #abline(-1800,60)   # no below
    bad=which( ((weight>length*300-8500)&(weight>length*70-1300)) | ((weight<length*300-14000)|(weight<length*60-1800)) )
    good=which(! ((weight>length*300-8500)&(weight>length*60-1000)) | ((weight<length*300-14000)|(weight<length*60-1800)) )
    plot(weight~length,main=paste("timepoint ",time,sep=""))
    points(weight[bad]~length[bad],col="red",pch=19)
  }
  
  if (time==2) {
    #  abline(-18000,500)  # no above
    #  abline(-3300,150)   # no above
    #  abline(-30000,500) # no below
    #  abline(-6000,150)   # no below
    bad=which( ((weight>length*500-18000)&(weight>length*150-3300)) | ((weight<length*500-30000)|(weight<length*150-6000)) )
    good=which(! ((weight>length*500-18000)&(weight>length*150-3300)) | ((weight<length*500-30000)|(weight<length*150-6000)) )
    plot(weight~length,main=paste("timepoint ",time,sep=""))
    points(weight[bad]~length[bad],col="red",pch=19)
  }
  
  if (time==3) {
    #  abline(-19000,500)  # no above
    #  abline(-3300,150)   # no above
    #  abline(-30000,500) # no below
    #  abline(-6000,150)   # no below
    bad=which( ((weight>length*500-19000)&(weight>length*150-3300)) | ((weight<length*500-30000)|(weight<length*150-6000)) )
    good=which(! ((weight>length*500-19000)&(weight>length*150-3300)) | ((weight<length*500-30000)|(weight<length*150-6000)) )
    plot(weight~length,main=paste("timepoint ",time,sep=""))
    points(weight[bad]~length[bad],col="red",pch=19)
  }
  
  if (time>3) {
    print("THE OUTLIERS WERE NOT MANUALLY DEFINED FOR THIS TIMEPOINT YET")
  }

  clean[bad,(2+1+(time-1)*2)]= NA ## AUTOMATIC ASSUMPTION THAT w and l are ordered in a particular way! ***
  clean[bad,(2+2+(time-1)*2)]= NA ## AUTOMATIC ASSUMPTION THAT w and l are ordered in a particular way! ***
  
  # collect the MEDIAN and SD  from the cleaned (outlier-free) dataset
  POINTS[time,1]= median( clean[,(2+1+(time-1)*2)], na.rm=T)
  POINTS[time,2]= sd( clean[,(2+1+(time-1)*2)], na.rm=T)
  POINTS[time,3]= median( clean[,(2+2+(time-1)*2)], na.rm=T)
  POINTS[time,4]= sd( clean[,(2+2+(time-1)*2)], na.rm=T)
  
  # save these outliers for the future use (when defining equally sized angular segments of the plot)
  assign(paste("a2bad_",time,sep=""),bad)
} # end of cycling through timepoints


###################################################################################
#  3  #####   DEFINE DISTRIBUTION OF delta-w/l for different in-between-timepoints
#             using a cleaned (outlier-free) dataset
###################################################################################

head(clean)
nt=(dim(clean)[2]-2)/2  # how many timepoints are there in a dataset
ng=nt-1  # how many gaps (inter-timepoints)  are there in a dataset

GAPS=matrix(NA,nrow=ng,ncol=4)
colnames(GAPS)=c("wght_dif_med","wght_dif_sd","lngth_dif_med","lngth_dif_sd")
rownames(GAPS)=paste("time_",1:ng,"-",2:(ng+1),sep="")
par(mfrow=c(2,2))
for (o in 1:ng) {
  wb= clean[,(2+1+(o-1)*2)]   # Weight Before
  wf= clean[,(2+1+(o)*2)]   # Weight aFter
  lb= clean[,(2+2+(o-1)*2)]   # Length Before
  lf= clean[,(2+2+(o)*2)]   # Length aFter
wd=wf-wb # weight difefrence
ld=lf-lb # length difference
 
### the sanity filter
#sum( wd<median(wb,na.rm=T)*0.1*-1, na.rm=T) # the 10% decrease in weight is increadible
#sum( ld<0, na.rm=T) # the 10% decrease in weight is increadible
wdb=which(wd<median(wb,na.rm=T)*0.1*-1) # Weight Diff Bad (10% decrease ini weight is increadible)
ldb=which(ld<0) # Length Dif Bad ( decrease in length is increadible)

hist(wd[-wdb],breaks=100,main=paste("weight dif ",o,"-",o+1," md=",median(wd,na.rm=T), sep=""))
hist(ld[-ldb],breaks=50,main=paste("length dif ",o,"-",o+1," md=",median(ld,na.rm=T),sep=""),col="grey")
GAPS[o,1]=median(wd[-wdb],na.rm=T)
GAPS[o,2]=sd(wd[-wdb],na.rm=T)
GAPS[o,3]=median(ld[-ldb],na.rm=T)
GAPS[o,4]=sd(ld[-ldb],na.rm=T)
}
GAPS

###################################################################################
#  4  #####   STANDARTIZE w/l-values at each timepoint
###################################################################################

#POINTS
#GAPS

# create a backbone for a STANDARTIZED dataset (of w/l values)
st=matrix(nrow=dim(a2)[1],ncol=dim(a2)[2]-2) # st= standartized
colnames(st)=colnames(a2)[-c(1,2)]
dim(st); head(st)

# fill the backbone with standartized values (while trying to make 0-tick not to coincide with a cluster of values)
par(mfrow=c(1,3))
for (u in 1:3) {  #  THIS VALUE IS SET MANUALLY, should be equal do timepoints in "POINTS", "GAPS", "clean" dataframes..)
  if(u==3) {dumindx=2} else {dumindx=u} # since there are less gaps than timepoints
  # standartized weight = ((weight) - (median in weight)) / (median in weight gain) - non-zero-match-coef
  nzmc=(1/GAPS[dumindx,1])/2 # non-zero-match-coef (zero coordnate should not coincide with cluster of data points)
  st[,1+(u-1)*2]= (a2[,2+1+(u-1)*2]-POINTS[u,1]) / GAPS[dumindx,1] - nzmc #/ GAPS[dumindx,2]
  # standartized length = ((length) - (median in length)) / (median in length gain)
  nzmc=(1/GAPS[dumindx,3])/2 # non-zero-match-coef (zero coordnate should not coincide with cluster of data points)
  st[,2+(u-1)*2]= (a2[,2+2+(u-1)*2]-POINTS[u,3]) / GAPS[dumindx,3] - nzmc #/ GAPS[dumindx,4]
#plot(st[,1+(u-1)*2]~st[,2+(u-1)*2])
#smoothScatter(st[,2+(u-1)*2],st[,1+(u-1)*2],colramp=palette,pch=1,cex=0.7,nbin=100,nrpoints=150,
#              main=paste("timepoint ",time,sep=""))
# make sure there are no value sthat are ==0  (due to tangent function downstream)
st[which(st[,1+(u-1)*2]==0),1+(u-1)*2] = rnorm( sum(st[,1+(u-1)*2]==0,na.rm=T),0,0.001)
st[which(st[,2+(u-1)*2]==0),2+(u-1)*2] = rnorm( sum(st[,2+(u-1)*2]==0,na.rm=T),0,0.001)
}

# check whether there still are any values equal zero
for(g in 1:6)   print(sum(st[,g]==0,na.rm=T))


for (h in 1:3) {....... #for each timepoint NOT FINISHED
  #h=1
  weight=st[,1+(h-1)*2]
  length=st[,2+(h-1)*2]



par(mfrow=c(1,1))
wgh=weight[-a2bad_1]  # a2bad_1 was defined earlier automatically when creating a "cleaned" file
lng=length[-a2bad_1]
plot(wgh~lng)


# create equaly spaced very small angular mini-beams
dgs=seq(0,360,360/39999) # degree-steps: the goal is to form "beams" with as equal number of individuals as possible
dgs[(dgs==90)|(dgs==180)|(dgs==270)]=NA # these might cause problems
dgs=dgs[!is.na(dgs)]; length(dgs); head(dgs)

# collector of counts of individuals in each mini-beam
qnt=rep(NA,length(dgs)) # the number (quantity) of individuals in each "beam"

# count individuals in each beam
for(v in 1:length(dgs)){
  if(v!=length(dgs)){
if ((dgs[v]<90)&(dgs[v+1]<90)) ids=which( (wgh>=lng*tan((dgs[v]/360)*2*pi)) & (wgh<lng*tan((dgs[v+1]/360)*2*pi)) )
if ( (dgs[v]<90)&(dgs[v+1]>90) ) ids=which( (wgh>=lng*tan((dgs[v]/360)*2*pi)) & (wgh>lng*tan((dgs[v+1]/360)*2*pi)) )
if ( (dgs[v]>90)&(dgs[v+1]<270)) ids=which( (wgh<=lng*tan((dgs[v]/360)*2*pi)) & (wgh>lng*tan((dgs[v+1]/360)*2*pi)) )
if ( (dgs[v]>90)&(dgs[v+1]>270)) ids=which( (wgh<=lng*tan((dgs[v]/360)*2*pi)) & (wgh<lng*tan((dgs[v+1]/360)*2*pi)) )
if ( (dgs[v]>270)&(dgs[v+1]>270)) ids=which( (wgh>=lng*tan((dgs[v]/360)*2*pi)) & (wgh<lng*tan((dgs[v+1]/360)*2*pi)) )
} else {
  ids=which( (wgh>=lng*tan((dgs[v]/360)*2*pi)) & (wgh<0) )
}
#points(wgh[ids]~lng[ids],pch=19,col=v)
qnt[v]=length(ids)
rm(ids)
}

# check whether all individuals were classified
#sum(qnt)
#sum((!is.na(wgh))&(!is.na(lng)))
#length(qnt)

qntid=NULL  # how many equal-sized degree-steps are there in one "beam"
ninds=NULL  # number of individuals in each "beam"  (the last "beam" might contain low number of individuals)
marks=NULL  # the index of each degree which ends the "beam" (counter-clocwise most distant edge)

mark=0
repeat {
if ((mark+5000)>length(qnt)) {mx=length(qnt)} else {mx=mark+5000} # you might need to increase this one in number of individuals in beams are too low
print(mark)
#mx
cll=NULL
for(x in (mark+1):mx) { cll=c(cll,sum(qnt[(mark+1):x])) }
if (sum(cll>1000,na.rm=T) >0) { 
  mark.new=min(which(cll>1000))+mark 
  qntid=c(qntid,min(which(cll>1000)))
  ninds=c(ninds,sum(qnt[(mark+1):mark.new]))
  marks=c(marks,mark.new)
}
if (sum(cll>1000,na.rm=T)==0) { 
  mark.new=length(cll)+mark
  marks=c(marks,mark.new)
  qntid=c(qntid,length(cll))
  ninds=c(ninds,sum(qnt[(mark+1):mark.new]))
}
if (sum(is.na(cll))==length(cll)) {break}
#print(mark.new)

if(mark.new==mark) {
  break
} else {
  mark=mark.new
}
}

length(marks)
length(qntid)
length(ninds)

# manually trim adjust the values in "marks"
marks=c(0,marks) # the first value shoud be zero
marks=marks[-((length(marks)-2):(length(marks)))] # the last value should be tripple-trimmed  # INSPECT MANUALLY !!!
qntid=qntid[-((length(ninds)-1):(length(ninds)))] # should be double-trimmed     INSPECT MANUALLY!!!
ninds=ninds[-((length(ninds)-1):(length(ninds)))] # should be double-trimmed
length(marks)
length(qntid)
length(ninds)

######  explore the generated "beams" and individuals in them
qntid  # how many equal-sized degree-steps are there in one "beam" (avoid the value 2000 - the cap defined previously)
ninds  # number of individuals in each "beam"  (the last "beam" might contain low number of individuals)
marks  # the index of each degree which ends the "beam" (counter-clocwise most distant edge)
dgs[marks]

length(qntid)
length(ninds)
length(marks)
hist(ninds,breaks=20,col="grey") #  all values should be very equal! if not - split into more mini-chunks
dgs[marks]


# draw
plot(wgh~lng)
abline(h=0)
abline(v=0)
# identify individuals who belong to each of the beams
NPDI=matrix(NA,nrow=2000,ncol=length(marks)) # NPDI = NonParametric DIstance table
for(k in 1:(length(marks))) {
deg1=dgs[marks[k+0]+1]
if(k==length(marks)) {deg2=360} else {deg2=dgs[marks[k+1]+0]}
#deg2
if ( (deg1<90)&(deg2<90) ) ids=which( (wgh>=lng*tan((deg1/360)*2*pi)) & (wgh<lng*tan((deg2/360)*2*pi)) )
if ( (deg1<90)&(deg2>90) ) ids=which( (wgh>=lng*tan((deg1/360)*2*pi)) & (wgh>lng*tan((deg2/360)*2*pi)) )
if ( (deg1>90)&(deg2<270) ) ids=which( (wgh<=lng*tan((deg1/360)*2*pi)) & (wgh>lng*tan((deg2/360)*2*pi)) )
if ( (deg1<270)&(deg2>270) ) ids=which( (wgh<=lng*tan((deg1/360)*2*pi)) & (wgh<lng*tan((deg2/360)*2*pi)) )
if ( (deg1>270)&(deg2>270) ) {
  if (k==length(marks)) { #assumption - marks contain one entry too many
    ids=which( (wgh>=lng*tan((deg1/360)*2*pi)) & (wgh<0) ) # in case of the last beam start - continue to end
  } else {
  ids=which( (wgh>=lng*tan((deg1/360)*2*pi)) & (wgh<lng*tan((deg2/360)*2*pi)) )
  }
}

NPDI[,k]=c(ids,rep(NA,2000-length(ids)))
points(wgh[ids]~lng[ids],pch=19,col=k)  # it is visible that some small amount of individuals are not classified (around lng=0)
print(paste(round(deg1,1),"-",round(deg2,1),sep=""))
}







... problem.. not all individuals get classified....




# preview
dim(NPDI)
NPDI[1:10,75:77]
NPDI[,77]
test=NULL; for(x in 1:77) {test=c(test,sum(!is.na(NPDI[,x])))}
sum(test)

valid=which(((!is.na(wgh))&(!is.na(lng))))
classified=NULL; for(x in 1:77) {classified=c(classified,NPDI[,x])}
classified=classified[!is.na(classified)]
length(classified)
unclasified=which(!valid %in% classified)
plot(wgh~lng)
abline(h=0)
abline(v=0)
points(wgh[unclasified]~lng[unclasified],pch=19,col="red")
head(unclasified)
use=unclasified[1:100]


... here above are unclassified indivduals ...





# find the direction toward which the cluster is pointing..
beta=coef(summary(lm(wgh~lng)))[2,1] # beta coefficient (aka tangent)
angl=atan(beta)/(2*pi)*360  # the angle at which the dataset is pointing
abline(0,beta)

# which beam should be the referrence beam
hbix=max(which(dgs[marks]<angl)) # the index of the head beam (assumption tha tthere are 78 beams)
length(marks)
dgs[marks]

# for each beam save the distance from the center
dists=(wgh^2+lng^2)^(1/2)
hist(dists,breaks=100,col="grey")

# see how the parametric distance looks like on plot 
plot(wgh~lng)
points(wgh[dists>1]~lng[dists>1],col="red",pch=20)


# see how the nonparametric distance looks like on plot
thr=0.999
plot(wgh~lng)
for(dd in 1:dim(NPDI)[2]){
  #print(dd)
  bid=NPDI[,dd] # ids from one beam
  bid=bid[!is.na(bid)]
  bdi=dists[bid] # distances from beam
  bwg=wgh[bid] # weights in beam
  bln=lng[bid] # lengths in beam
  chu=data.frame(d=bdi,w=bwg,l=bln) # a chunk of data that belongsto a beam
  #head(chu)
  chu=chu[order(chu$d),]
  maxd= chu[floor(dim(chu)[1]*thr),"d"] # maximum distance from the center
  y=chu[chu$d>maxd,"w"]
  x=chu[chu$d>maxd,"l"]
  points(y~x,col="red",pch=19)
}


# 







dim(NPDI)
NPDI[1:10,1:10]











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







