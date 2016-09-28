


cat("
    creates two flag-matrixes:

hmsk3 = a mask based on GAweek-stratified distribution of 
        head circumferance and its outliers, based on variable 
        number of standard deviations

hwmsk3 = a mask based on residuals in the 
          weight ~  head circumferance plot 
          (also stratified by th egestational week)
          (also based on variable number of standard deviations)

also saves intermediate plots that visualizes flags and their overlaps

in the future it might be more reasonable to use the 
length variable for residuals instead of weight, since
the length has less missing values after the previous cleaning steps
    ")

GAbasedUnivarBivarOutlierDetectionHEAD=function(a,weig,leng,head)  {
  
  hmsk3 <<- matrix(0,nr=dim(a)[1],nc=length(weig)) # mask for head circumferance (univar, standard deviations)
  hwmsk3 <<- matrix(0,nr=dim(a)[1],nc=length(leng))  # paired-value mask (residuals)  head-weight
  #hlmsk3 <<- matrix(0,nr=dim(a)[1],nc=length(leng))  # paired-value mask (residuals)   head-length
  
  
  
  for (time in 1:length(weig)) { # cycle through various time points
  
    print(time)
    h1=a[,head[time]]
    w1=a[,weig[time]]
    l1=a[,leng[time]]
    
    GAW=floor(a$GestationalAge/7)
    range(GAW,na.rm=T)
    GAW[GAW<25]=NA  # 25 is where the survival is possible
    GAW[GAW>43]=NA  # upper limit to GA values
    
    
    #######################################
    ###  estimate residuals for each individual 
    ###  while stratifying on gestational age week
    ###  using a clean dataset (no outliers or high-impact values)
    
    
    # the matrix of thresholds (as standard deviations) 
    # in case we want to apply individual thresholds for each week or each timepoint of age (better)
    
    mthrs=array(data=NA,dim=c(length(head),(43-25+1),2))
    
    # threshold for "residuals" method (bivariate)
    mthrs[1,,1] = seq(from=3,to=5,by=(5-3)/(19-1))
    mthrs[2,,1] = seq(from=3,to=5,by=(5-3)/(19-1))
    mthrs[3,,1] = seq(from=3,to=5,by=(5-3)/(19-1))
    mthrs[4,,1] = seq(from=3,to=5,by=(5-3)/(19-1))
    mthrs[5,,1] = seq(from=3,to=5,by=(5-3)/(19-1))
    mthrs[6,,1] = seq(from=3,to=5,by=(5-3)/(19-1))
    
    # threshold for weight as univariate
    mthrs[1,,2] = seq(from=3,to=5,by=(5-3)/(19-1))
    mthrs[2,,2] = seq(from=3,to=5,by=(5-3)/(19-1))
    mthrs[3,,2] = seq(from=3,to=5,by=(5-3)/(19-1))
    mthrs[4,,2] = seq(from=3,to=5,by=(5-3)/(19-1))
    mthrs[5,,2] = seq(from=3,to=5,by=(5-3)/(19-1))
    mthrs[6,,2] = seq(from=3,to=5,by=(5-3)/(19-1))
    
    
    LOCATION="C:/Users/Jon/Dropbox/GIT/CLEANING_MOBA_FILES/Childs_HghWghCrc/5_PDFs/"
    pdf(file=paste(LOCATION,"HEAD_intermediate_time",time,".pdf",sep=""))
    par(mfrow=c(3,3))
    
    DAT=NULL
    betas=NULL
    
    wID=0 # WEEK id
    for (w in 25:43) {  # for each gestational week
    wID=wID+1
    
      sel=which(GAW==w) # sel=selection
      n.valid.data.points=sum((!is.na(w1[sel]))&(!is.na(h1[sel])))
      
      
      
      if (n.valid.data.points<4) {
      dat=data.frame(w=w1[sel],l=l1[sel],t1=1,t2=1)  
      DAT=rbind(DAT,dat)
      betas=c(betas,NA)
      plot(1, type="n", axes=F, xlab="", ylab="")
      plot(1, type="n", axes=F, xlab="", ylab="")
      plot(1, type="n", axes=F, xlab="", ylab="")
        
      } else {
      
    
        # gest.age weeks with too many individuals (not time efficient)
      if (n.valid.data.points>100) {
       bad1=unique(which(abs(h1[sel]-mean(h1[sel],na.rm=T))>3*sd(h1[sel],na.rm=T))) # exclude head-circum outliers
       bad2=unique(which(abs(w1[sel]-mean(w1[sel],na.rm=T))>3*sd(w1[sel],na.rm=T))) # exclude weight outliers
       bad=unique(c(bad1,bad2)) 
        # prevent zero-length exclusion list
      if (length(bad)==0) {y=w1[sel];x=h1[sel]} else {y=w1[sel][-bad];x=h1[sel][-bad]}
        # estimate the regression model and its coeficients
      mod=lm(y~x) # model (for that particular gestational week)
      beta=round(coef(summary(mod))[2,1],0)
        # visualize the regression line estimated without outliers
      COL=rep(1,length(sel)); COL[bad]=3 # colors
      PCH=rep(1,length(sel)); PCH[bad]=19 # colors
      plot(1, type="n", axes=F, xlab="", ylab="") # since there is no residual-based plot (too time consuming to generate)
      plot(w1[sel]~h1[sel],main=paste(w," (beta=",beta,")",sep=""),col=COL,pch=PCH)
      abline(mod) # NOTE!  model (mod) will be used later on 
                                    }
     
     
     
      # gest.age weeks with low number of individuals (time efficient to estimate impact of each datapoint)
        if (n.valid.data.points<=100) { # the minimum count of datapoints needed for estimations
      # estimate individual impact factor (the effect on the beta, when value is excluded)
        impct=rep(NA,length(sel)) # impact
              rm(mod)
          for (z in 1:length(sel))  {
            mod=lm(w1[sel][-z]~h1[sel][-z])
            impct[z]=coef(summary(mod))[2,1]
                                    }
          # when paired values are missing - impact shoud be non-existent (redundant rule?)
          impct[is.na(w1[sel])|is.na(h1[sel])]=NA
          
          bad0=which(abs(impct-mean(impct,na.rm=T))>(3*sd(impct,na.rm=T))) # exclude residual-based outliers
          bad1=which(abs(h1[sel]-mean(h1[sel],na.rm=T))>3*sd(h1[sel],na.rm=T)) # exclude univar (head circumf) outliers
          bad2=which(abs(w1[sel]-mean(w1[sel],na.rm=T))>3*sd(w1[sel],na.rm=T)) # exclude univar (weight) outliers
          bad=unique(c(bad0,bad1,bad2))
          
          # visualize the high-impact values (based on residuals)
          COL=rep(1,length(sel));  COL[bad0]=2  # color residual-based outliers
          COL[unique(c(bad1,bad2))]=3  # color univar-based outliers
          COL[bad0[ which(bad0 %in% unique(c(bad1,bad2)))]]=4  # when both residuals and univar outliers give flags
          PCH=rep(1,length(sel)); PCH[bad]=19
          plot(impct~seq(length(sel)),main=paste(w," (impact)",sep=""),col=COL, pch=PCH,
               ylab="Beta (head circumf ~ weight)",xlab="Individual") #impact factors
          abline(h=mean(impct,na.rm=T))
          
          # visualize the regression line estimated without outliers
          if (length(bad)==0) {y=w1[sel];x=h1[sel]} else {y=w1[sel][-bad];x=h1[sel][-bad]}
          mod=lm(y~x); beta=round(coef(summary(mod))[2,1],0)
          plot(w1[sel]~h1[sel],main=paste(w," (beta=",beta,")",sep=""),col=COL,pch=PCH,xlab="Head circumferance",ylab="Weight")
          abline(mod)  # NOTE!  mod (model) will be used later on
                                          }
      
      
        
        # save beta values for plotting
        betas=c(betas,beta)
        
        rm(bad,bad0,bad1,bad2)
        
        ## RESIDUAL FILTER
        resid=w1[sel]-as.numeric(predict(mod,data.frame(x=h1[sel])))
        bad0=which(abs(resid-mean(resid,na.rm=T))>(mthrs[time,wID,1]*sd(resid,na.rm=T)))
                
        ## HEAD-UNIVAR FILTER
        bad1=which(abs(h1[sel]-mean(h1[sel],na.rm=T))>mthrs[time,wID,2]*sd(h1[sel],na.rm=T))
        
        ## only unique bad datapoints 
        bad=unique(c(bad0,bad1))
        
        # visualisation parameters
        COL=rep(1,length(sel)); COL[bad]=2
        PCH=rep(1,length(sel)); PCH[bad]=19
        
        # visualize exclusions
        plot(w1[sel]~h1[sel],main=paste(w," (final exclusions)",sep=""),col=COL,pch=PCH,
             xlab="Head circumferance",ylab="Weight")
        
        # visualize thresholds based on residuals
        rs=as.numeric(mthrs[time,wID,1]*sd(resid,na.rm=T))
        mi=min(h1[sel],na.rm=T); ma=max(h1[sel],na.rm=T)
        xv=seq(mi,ma,by=1) # x-axes values
        up=as.numeric(predict(mod,data.frame(x=xv)))+rs
        lw=as.numeric(predict(mod,data.frame(x=xv)))-rs
        lines(up~xv,lty=2)
        lines(lw~xv,lty=2)
        
        # visualize thresholds based on head circumferance univar
        up=mean(h1[sel],na.rm=T) + as.numeric(mthrs[time,wID,2]*sd(h1[sel],na.rm=T))
        lw=mean(h1[sel],na.rm=T) - as.numeric(mthrs[time,wID,2]*sd(h1[sel],na.rm=T))
        abline(v=c(up,lw),lty=2)
        
        # for saving the flags
        if (length(bad0)>0)  hwmsk3[sel[bad0],time] <<- 1 #(global)
        if (length(bad1)>0)  hmsk3[sel[bad1],time] <<- 1 #(global)
        
        # for plotting the flags
        dat=data.frame(w=w1[sel],h=h1[sel],t1=1,t2=1)  
        dat$t1[bad0]=2
        dat$t2[bad1]=2
        #dat$t3[bad2]=2
        DAT=rbind(DAT,dat)
      
          
      } # end of else
    } # end of cycling through gestational weeks
    dev.off()
    

    
    LOCATION=("C:/Users/Jon/Dropbox/GIT/CLEANING_MOBA_FILES/Childs_HghWghCrc/5_PDFs/")
    jpeg(file=paste(LOCATION,"HEAD_exclusions_time",time,".jpeg",sep=""),
         width = 30,height = 25, quality=100,units="cm", res=150)
    
    
    par(mfrow=c(2,3))
    plot(DAT$w~DAT$h,col=DAT$t1,main="residuals",xlab="Head",ylab="Weight")
    plot(DAT$w~DAT$h,col=DAT$t2,main="weight SD",xlab="Head",ylab="Weight")
    plot(DAT$w~DAT$h,col=(DAT$t1+DAT$t2)-1,main="combined filters",
         xlab="Head",ylab="Weight")
    bbb=which((DAT$t1+DAT$t2)>2)
    plot(DAT$w[-bbb]~DAT$h[-bbb],main="cleaned",xlab="Head",ylab="Weight")
    plot(betas~seq(length(betas)),type="l",xlab="ID of Gestational week",ylab="Beta (Head~Weight)",
         main="quality check for thresholds")
    dev.off()
    
  } # end of cycling through time points
  
} # end of function

