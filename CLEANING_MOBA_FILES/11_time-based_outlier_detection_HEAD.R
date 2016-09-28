

# time-based outlier detection (flagging)


TimeBasedOutlierDetectionHEAD=function(a,head) {
  

  print("finding errors for head")
  h=a[,head]
  nr=dim(h)[1]
  nc=dim(h)[2]
  m1=matrix(NA,nr=nr,nc=nc)
  
  # fill the matrix with error flags
  for (cl in 1:nc) {
    print(paste("running timepoint ",cl,sep=""))  
    print(Sys.time())
    if (cl==1)  { for (rw in 1:nr) { m1[rw,1]=  sum( (h[rw,1] >= as.numeric(h[rw,(2:nc)])),na.rm=T) } }
    if ((cl>1)&(cl<nc)) {  
      for (rw in 1:nr) {
        m1[rw,cl]= sum( (h[rw,cl] >= as.numeric(h[rw,(cl+1):nc])),na.rm=T) + 
          sum( (h[rw,cl] <= as.numeric(h[rw,1:(cl-1)])),na.rm=T)
      }
    }
    if (cl==nc) { for (rw in 1:nr) { m1[rw,nc]= sum( (h[rw,nc] <= as.numeric(h[rw,1:(nc-1)])),na.rm=T) } }
  }
  
  
  hmaxs <<- apply(m1,1,max)
  herrs <<- apply(m1,1,sum)
  # before the deletion 
  print(table(herrs,hmaxs))
  
  par(mfrow=c(4,5));  for (i in 1:20) {   y=as.numeric(a[which(herrs==8),head][i,]); plot(y~seq(6)) }
  
  # which ones to delete
  for (i in which(hmaxs>1) ) h[i,which(m1[i,]==hmaxs[i])]=NA
  
  # which ones to delete
  hmsk2=matrix(0,nr=nr,nc=nc)  # first version
  for (i in which(hmaxs==1) ) hmsk2[i,]=as.numeric(m1[i,]==hmaxs[i])
  rm(m1)
  
  hmsk2  <<- hmsk2
  a[,head] <<- h   # not sure whether this line works
  rm(h)
    
  cat("hmsk2 was created")
}

