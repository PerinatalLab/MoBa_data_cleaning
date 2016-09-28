

cat("
    creating the matrix with flags for increadible values of head-circumferance
    ")
# for HEAD circumferance

correctOCRandTyposHEAD=function(a,weig,leng,head) { 
  
  # manually picked lower/upper limits for reasonable values of height and weigth at each timepoint
  hlw=c(10,20,20,30,30,30) # upper for head circumferance
  hup=c(45,50,55,60,60,60) # lower for head circumferance
  
  # collector of flags (insane values)
  hmsk1 <<- matrix(0,nr=dim(a)[1] ,nc=length(head)) # global assignment
  
  par(mfrow=c(2,3))
  for (time in 1:6){
    print(paste(head[time],sep=""))
    h1=a[,head[time]]
    w1=a[,weig[time]]
    l1=a[,leng[time]]
    
#    plot(w1~h1,main=time,xlab=head[time],ylab=weig[time])
     plot(l1~h1,main=time,xlab=head[time],ylab=leng[time])
      abline(v=c(hup[time],hlw[time]),lty=2)
   
####################################################
## SECOND PART: MARK FOR DELETION THE INSANE VALUES
hmsk1[which(h1>hup[time]),time] <<- 1 #(gloabl)
hmsk1[which(h1<hlw[time]),time] <<- 1 #(gloabl)

####################################################
## SECOND PART: SIMPLY DELETE THE INSANE VALUES
#  a[ which(h1>hup[time]) ,head[time]] <<- NA # (global)
#  a[ which(h1<hlw[time]) ,head[time]] <<- NA # (global)

  } # end of cycling through time-points
  
cat("hmsk1  matrix was created")

} # end of function
