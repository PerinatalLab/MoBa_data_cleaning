

GAbasedUnivarBivarOutlierDetection=function(a,weights,lengths)  {

  
  wmsk3 <<- matrix(0,nr=dim(a)[1],nc=length(weights)) # mask for weights
  lmsk3 <<- matrix(0,nr=dim(a)[1],nc=length(lengths)) # mask for lengths
  pmsk3 <<- matrix(0,nr=dim(a)[1],nc=length(lengths))  # paired-value mask (residuals)

for (time in 1:length(weights)) { # cycle through various time points
  print(time)
  
  w1=a[,weights[time]]
  l1=a[,lengths[time]]
  
  GAW=floor(a$GestationalAge/7)
  GAW[GAW<25]=NA  # 25 is where the survival is possible
  GAW[GAW>43]=NA  # upper limit to GA values
  
    
#######################################
###  estimate residuals for each individual 
###  while stratifying on gestational age week
###  using a clean dataset (no outliers or high-impact values)
  
  # (in the future: include the week 20-24 in the plots ?)
  
  #plot(a$WeightBirthQ4~a$LengthBirthQ4)
  
  # the matrix of thresholds (as standard deviations) 
  # in case we want to apply individual thresholds for each week or each timepoint of age (better)
  mthrs=matrix(NA,nr=(43-25+1),nc=3)
  colnames(mthrs)=c("resid","weigh","lengt")
  #mthrs[,1]=seq(from=2,to=6,by=1/(dim(mthrs)[1]-1)*(6-2))
  #mthrs[,2]=seq(from=3,to=6,by=1/(dim(mthrs)[1]-1)*(6-3))
  #mthrs[,3]=seq(from=3,to=6,by=1/(dim(mthrs)[1]-1)*(6-3))
  mthrs[,1]=5  # threshold for "residuals" method (bivariate)
  mthrs[,2]=5  # threshold for weight as univariate
  mthrs[,3]=5  # threshold for length as univariate
  
  LOCATION="C:/Users/Jon/Dropbox/GIT/CLEANING_MOBA_FILES/Childs_HghWghCrc/5_PDFs/"
  pdf(file=paste(LOCATION,"intermediate_time",time,".pdf",sep=""))
  par(mfrow=c(3,3))
  
  DAT=NULL
  betas=NULL

  for (w in 25:43) {  # for each gestational week
    sel=which(GAW==w) # sel=selection
    
    # gest.age weeks with too many individuals (not time efficient)
    n.valid.data.points=sum((!is.na(w1[sel]))&(!is.na(l1[sel])))
    if (n.valid.data.points>100) {
      
      bad1=which(abs(w1[sel]-mean(w1[sel],na.rm=T))>2*sd(w1[sel],na.rm=T)) # exclude weight outliers
      bad2=which(abs(l1[sel]-mean(l1[sel],na.rm=T))>2*sd(l1[sel],na.rm=T)) # exclude length outliers
      bad=unique(c(bad1,bad2))
      
      # prevent zero-length exclusion list
      if (length(bad)==0) {y=w1[sel];x=l1[sel]} else {y=w1[sel][-bad];x=l1[sel][-bad]}
      
      # estimate the regression model and its coeficients
      mod=lm(y~x) # model (for that particular gestational week)
      beta=round(coef(summary(mod))[2,1],0)
      
      # visualize the regression line estimated without outliers
      COL=rep(1,length(sel)); COL[bad]=2 # colors
      plot(w1[sel]~l1[sel],main=paste(w," (beta=",beta,")",sep=""),col=COL)
      abline(mod)
      
      # weeks with low number of individuals (time efficient to estimate impact of each datapoint)
    } else {
      
      if (n.valid.data.points>=4) { # the minimum count of datapoints needed for estimations
                    
      # estimate individual impact factor (the effect on the beta, when value is excluded)
      impct=rep(NA,length(sel)) # impact
      for (z in 1:length(sel))  {
        rm(mod)
        mod=lm(w1[sel][-z]~l1[sel][-z])
        impct[z]=coef(summary(mod))[2,1]
      }
      # when paired values are missing - impact shoud be non-existent (redundant rule?)
      impct[is.na(w1[sel])|is.na(l1[sel])]=NA
      
      bad0=which(abs(impct-mean(impct,na.rm=T))>(2*sd(impct,na.rm=T))) # exclude residual outliers
      bad1=which(abs(w1[sel]-mean(w1[sel]))>2*sd(w1[sel])) # exclude weight outliers
      bad2=which(abs(l1[sel]-mean(l1[sel]))>2*sd(l1[sel])) # exclude length outliers
      bad=unique(c(bad0,bad1,bad2))
      
      # visualize the high-impact values (based on residuals)
      COL=rep(1,length(sel)); COL[bad0]=2
      plot(impct~seq(length(sel)),main=paste(w," (impact)",sep=""),col=COL,
           ylab="Beta (weight~length)",xlab="Individual") #impact factors
     abline(h=mean(impct,na.rm=T))
      
      # visualize the regression line estimated without outliers
      if (length(bad)==0) {y=w1[sel];x=l1[sel]} else {y=w1[sel][-bad];x=l1[sel][-bad]}
      mod=lm(y~x)
      beta=round(coef(summary(mod))[2,1],0)
      COL=rep(1,length(sel)); COL[bad]=2
      plot(w1[sel]~l1[sel],main=paste(w," (beta=",beta,")",sep=""),col=COL,
           xlab="Length",ylab="Weight")
      abline(mod)
    }
    }
    
    
    if (n.valid.data.points>=4) {
    par(mfrow=c(1,1))
    plot(w1[sel]~l1[sel])
    
    # save beta values for plotting
    betas=c(betas,beta)
    
    rm(bad,bad0,bad1,bad2)
    
    ## RESIDUAL FILTER
    resid=w1[sel]-as.numeric(predict(mod,data.frame(x=l1[sel])))
    bad0=which(abs(resid-mean(resid,na.rm=T))>(mthrs[time,"resid"]*sd(resid,na.rm=T))) 
    
    ## WEIGHT FILTER
    bad1=which(abs(w1[sel]-mean(w1[sel],na.rm=T))>mthrs[time,"weigh"]*sd(w1[sel],na.rm=T))
    
    ##  LENGTH FILTER
    bad2=which(abs(l1[sel]-mean(l1[sel],na.rm=T))>mthrs[time,"lengt"]*sd(l1[sel],na.rm=T))
    
    ## only unique bad datapoints 
    bad=unique(c(bad0,bad1,bad2))
    
    # visualize exlusions
    COL=rep(1,length(sel)); COL[bad]=2
    plot(w1[sel]~l1[sel],main=paste(w," (final exclusions)",sep=""),col=COL,
         xlab="Length",ylab="Weight")
    
    # for saving the flags
    if (length(bad0)>0)  pmsk3[sel[bad0],time] <<- 1 #(global)
    if (length(bad1)>0)  wmsk3[sel[bad1],time] <<- 1 #(global)
    if (length(bad2)>0)  lmsk3[sel[bad2],time] <<- 1 #(global)
    
    # for plotting the flags
    dat=data.frame(w=w1[sel],l=l1[sel],t1=1,t2=1,t3=1)  
    dat$t1[bad0]=2
    dat$t2[bad1]=2
    dat$t3[bad2]=2
    DAT=rbind(DAT,dat)
  
    } else {
    dat=data.frame(w=w1[sel],l=l1[sel],t1=1,t2=1,t3=1)  
    DAT=rbind(DAT,dat)
    betas=c(betas,NA)
    }
  
  }
  
  dev.off()
 


  #head(DAT)
  #table(DAT$t1)
  
  LOCATION="C:/Users/Jon/Dropbox/GIT/CLEANING_MOBA_FILES/Childs_HghWghCrc/5_PDFs/"
  jpeg(file=paste(LOCATION,"exclusions_time",time,".jpeg",sep=""),
       width = 30,height = 25, quality=100,units="cm", res=150)

  par(mfrow=c(2,3))
  plot(DAT$w~DAT$l,col=DAT$t1,main="residuals",xlab="Length",ylab="Weight")
  plot(DAT$w~DAT$l,col=DAT$t2,main="weight SD",xlab="Length",ylab="Weight")
  plot(DAT$w~DAT$l,col=DAT$t3,main="length SD",xlab="Length",ylab="Weight")
  plot(DAT$w~DAT$l,col=(DAT$t1+DAT$t2+DAT$t3)-2,main="combined filters",
       xlab="Length",ylab="Weight")
  bbb=which((DAT$t1+DAT$t2+DAT$t3)>3)
  plot(DAT$w[-bbb]~DAT$l[-bbb],main="cleaned",xlab="Length",ylab="Weight")
  plot(betas~seq(length(betas)),type="l",xlab="Gestational week",ylab="Beta (Weight~Length)",
       main="quality check for thresholds")
  dev.off()
  
} # end of cycling through time points

} # end of function


