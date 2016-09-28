

matpatweigheig=function(a) {

  cat(" these libraries are needed:
library(SDMTools)
library(hexbin)
      ")

mh= a$MaternalHeight
mw= a$MaternalWeight
ph= a$PaternalHeight

# descriptives
#sum(is.na(mh))
#sum(is.na(mh))/length(mh)
#sum(is.na(mw))
#sum(is.na(mw))/length(mw)
#sum(is.na(ph))
#sum(is.na(ph))/length(ph)

# preview the problems
jpeg("./Childs_HghWghCrc/20_PARENTAL_before_cleaning_univar.jpeg",
     width = 30,height = 25, quality=100,units="cm", res=150)
par(mfrow=c(1,3))
plot(ph,main="paternal height",ylab="cm")
plot(mh,main="maternal height",ylab="cm")
plot(mw,main="maternal weight",ylab="kg")
dev.off()

# preview the problems
jpeg("./Childs_HghWghCrc/20_PARENTAL_before_cleaning_bivar.jpeg",
     width = 30,height = 10, quality=100,units="cm", res=150)
par(mfrow=c(1,3))
plot(ph,mh,xlab="paternal height",ylab="maternal height")
plot(mh,mw,xlab="maternal height",ylab="maternal weight")
plot(ph,mw,xlab="paternal height",ylab="maternal weight")
dev.off()


# preview the problems
jpeg("./Childs_HghWghCrc/20_PARENTAL_before_cleaning_univar_histogram.jpeg",
     width = 30,height = 10, quality=100,units="cm", res=150)
par(mfrow=c(1,3))
hist(ph,main="paternal height",xlab="cm",breaks=100,col="grey")
hist(mh,main="maternal height",xlab="cm",breaks=100,col="grey")
hist(mw,main="maternal weight",xlab="kg",breaks=100,col="grey")
dev.off()



# the thresholds for increadible values 
mhup=210 # maternal height upper threshold
mhlw=130 # maternal height lower threshold
phup=220 # paternal height upper threshold
phlw=130 # paternal height lower threshold
mwup=200 # maternal weight upper threshold
mwlw=30 # maternal weight lower threshold


# preview the problems
jpeg("./Childs_HghWghCrc/20_PARENTAL_before_cleaning_increadible-values_bivar.jpeg",
     width = 30,height = 10, quality=100,units="cm", res=150)
par(mfrow=c(1,3))
plot(ph,mh,xlab="paternal height",ylab="maternal height",xlim=c(0,phup),ylim=c(0,mhup))
abline(v=c(phlw,phup),col="red"); abline(h=c(mhlw,mhup),col="red") 

plot(mh,mw,xlab="maternal height",ylab="maternal weight",xlim=c(0,240))
abline(v=c(mhlw,mhup),col="red"); abline(h=c(mwlw,mwup),col="red") 

plot(ph,mw,xlab="paternal height",ylab="maternal weight",xlim=c(0,phup))
abline(v=c(phlw,phup),col="red"); abline(h=c(mwlw,mwup),col="red") 
dev.off()


################################################
##########   BIVARIATE METHOD (only for mothers)
################################################

# only used for plotting, does not perform corrections
# is redundand with the univariate method below

bin <- hexbin(mh,mw,xbins=100)
x=attributes(bin)$xcm # coordinates of hexagons
y=attributes(bin)$ycm # coordinates of hexagons
cc=attributes(bin)$count # count of individuals in each hexagon
par(mfrow=c(1,1))
plot(bin)

THR=.99

# find out what minimum number of individuals in one bin would include 99% of the population
ccc=cc[rev(order(cc))] ; head(ccc) # numbers of individuals in hexbins, ordered from the highest
te=rep(NA,length(ccc)); for (i in 1:length(ccc)) te[i]=sum(ccc[1:i]) # cumulative amount
MI=min(which(te> sum(cc)*THR)) # when is the earliest when 99% is reached (MI=MIN ID of Hexbin)
NI=ccc[MI] # the number of individuals in that threshold-hexbin (NI= number of individuals)
PCH=rep(1,length(x)); PCH[cc>NI]=19
plot(x,y,pch=PCH)

xx=x[cc>NI] # only the interesting hexagons
yy=y[cc>NI] # only the interesting hexagons

fun=function(xx,yy,mh,mw) {
  
  # define convex hull
  xtr=chull(xx, yy); xtr=c(xtr,xtr[1])
  #define the points and polygon
  pnts=data.frame(heigth=mh,weight=mw)
  poly=data.frame(heigth=xx[xtr],weight=yy[xtr])
  polygon(poly,col="grey")
  #create check which points fall within the polygon
  insiders = pnt.in.poly(pnts,poly)
  #identify points not in the polygon with an X
  points(insiders[which(insiders$pip==1),1:2])
}


jpeg("./Childs_HghWghCrc/20_PARENTAL_maternal_typos-OCRs_bivar.jpeg",
     width = 30,height = 20, quality=100,units="cm", res=150)
par(mfrow=c(1,1))
plot(x,y,pch=PCH,xlim=c(0,max(x)),ylim=c(0,max(y)),cex=0.1)
fun((xx/10),(yy/10),mh,mw)
fun((xx/10),yy,mh,mw)
fun(xx,(yy/10),mh,mw)
fun(xx,(yy*10),mh,mw)
fun((xx*10),yy,mh,mw)
fun(xx-100,yy,mh,mw)
fun(xx-100,yy/10,mh,mw)
fun(xx-100,yy*10,mh,mw)
fun(xx/10,yy*10,mh,mw)
dev.off()



###########################################
##########   UNIVARIATE METHOD
## only to be applied after bivariate method
## (after checking weight-height plots)

vals=c("mh","ph","mw")
vups=c("mhup","phup","mwup")
vlws=c("mhlw","phlw","mwlw")
nms=c("maternal height","paternal height","maternal weight")



jpeg("./Childs_HghWghCrc/20_PARENTAL_typos-OCRs_univar.jpeg",
     width = 30,height = 20, quality=100,units="cm", res=150)
par(mfrow=c(3,1))
for (r in 1:3) {
# GET THE VALUES
val=get(vals[r]) 
vup=get(vups[r]) 
vlw=get(vlws[r]) 

valp=val[which((val>vlw)&(val<vup))] # pure values (without insane outliers)
valpmd=median(valp)
valp=valp-valpmd # centered to zero
valp1=valp[valp<=0]*(-1) # left side
valp2=valp[valp>=0]      # rihgt side
valp1=valp1[order(valp1)] 
valp2=valp2[order(valp2)]

# THRESHOLDS FOR "MOST NORMAL" VALUES
thr1 =valpmd-valp1[floor(length(valp1)*0.99)] # nonparametric center      # ARBITRARY HERE !!
thr2 =valpmd+valp2[floor(length(valp2)*0.99)] # nonparametric center      # ARBITRARY HERE !!

# PLOTTING PART
hist(val,breaks=100,col="grey",main=nms[r])
abline(v=c(thr1,thr2),col="red")
abline(v=c(thr1*10,thr2*10),col="red")
abline(v=c(thr1/10,thr2/10),col="red")
abline(v=c(thr1-100,thr2-100),col="red")

# SELECTING PART
use1=which((val>thr1*10)&(val<thr2*10))
use2=which((val>thr1/10)&(val<thr2/10))
use3=which((val>thr1-100)&(val<thr2-100)&(val>=10)&(val<100))

# REPORTING PART
print(paste(length(use1)," 10xhigher values",sep=""))
print(paste(length(use2)," 10xlower values",sep=""))
print(paste(length(use3)," 100 units values",sep=""))

# FIXING PART
val[use1]=val[use1]/10
val[use2]=val[use2]*10
val[use3]=val[use3]+100

# DELETE NONSENSES
val[which((val<vlw)|(val>vup))]=NA

# SAVE
assign(paste(vals[r],sep=""),val)
}
dev.off()




# FINAL CHECKUP: VALIDITY OF CORRECTIONS


jpeg("./Childs_HghWghCrc/20_PARENTAL_after_cleaning_univar.jpeg",
     width = 30,height = 20, quality=100,units="cm", res=150)
par(mfrow=c(3,1))
nms=c("maternal height","paternal height","maternal weight")
for (r in 1:3) {
  val=get(vals[r]) 
  hist(val,breaks=100,col="grey",main=nms[r])
  }
dev.off()

jpeg("./Childs_HghWghCrc/20_PARENTAL_after_cleaning_bivar.jpeg",
     width = 30,height = 20, quality=100,units="cm", res=150)
par(mfrow=c(1,2))
plot(mh,mw,xlab="maternal height",ylab="maternal weight")
plot(mh,ph,xlab="maternal height",ylab="paternal height")
dev.off()

a$MaternalHeight <<- mh # global level
a$MaternalWeight <<- mw
a$PaternalHeight <<- ph

#sum(is.na(mh))
#sum(is.na(mh))/length(mh)
#sum(is.na(mw))
#sum(is.na(mw))/length(mw)
#sum(is.na(ph))
#sum(is.na(ph))/length(ph)




} # end of function


###########
### SANDBOX





