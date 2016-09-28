

#time-based outlier detection (flagging)



TimeBasedOutlierDetection=function(a,weights,lengths) {

########################################################
#################   LENGTH   ###########################
print("finding errors for length")
  lng=a[,lengths]
  nr=dim(lng)[1]
  nc=dim(lng)[2]
  m1=matrix(NA,nr=nr,nc=nc)

# fill the matrix with error flags
for (cl in 1:nc) {
  print(paste("running timepoint ",cl,sep=""))  
  print(Sys.time())
    if (cl==1)  { for (rw in 1:nr) { m1[rw,1]=  sum( (lng[rw,1] >= as.numeric(lng[rw,(2:nc)])),na.rm=T) } }
    if ((cl>1)&(cl<nc)) {  
       for (rw in 1:nr) {
        m1[rw,cl]= sum( (lng[rw,cl] >= as.numeric(lng[rw,(cl+1):nc])),na.rm=T) + 
        sum( (lng[rw,cl] <= as.numeric(lng[rw,1:(cl-1)])),na.rm=T)
                        }
                        }
    if (cl==nc) { for (rw in 1:nr) { m1[rw,nc]= sum( (lng[rw,nc] <= as.numeric(lng[rw,1:(nc-1)])),na.rm=T) } }
                        }


lmaxs <<- apply(m1,1,max)
lerrs <<- apply(m1,1,sum)
# before the deletion 
print(table(lerrs,lmaxs))

# which ones to delete
for (i in which(lmaxs>1) ) lng[i,which(m1[i,]==lmaxs[i])]=NA

# which ones to delete
lmsk2=matrix(0,nr=nr,nc=nc)  # first version
for (i in which(lmaxs==1) ) lmsk2[i,]=as.numeric(m1[i,]==lmaxs[i])
rm(m1)

lmsk2  <<- lmsk2
a[,lengths] <<- lng   # not sure whether this line works
rm(lng)


############################################################################
#################   WEIGHT   #######################


print("finding errors for weight")
  wgh=a[,weights]
  nr=dim(wgh)[1]
  nc=dim(wgh)[2]
  m2=matrix(NA,nr=nr,nc=nc)



for (cl in 1:nc) {
  print(paste("running timepoint ",cl,sep=""))  
  print(Sys.time())
  if (cl==1)  { for (rw in 1:nr) { m2[rw,1]=  sum( (wgh[rw,1] >= as.numeric(wgh[rw,(2:nc)])),na.rm=T) } }
  if ((cl>1)&(cl<nc)) {  
    for (rw in 1:nr) {
      m2[rw,cl]= sum( (wgh[rw,cl] >= as.numeric(wgh[rw,(cl+1):nc])),na.rm=T) + 
        sum( (wgh[rw,cl] <= as.numeric(wgh[rw,1:(cl-1)])),na.rm=T)
    }
  }
  if (cl==nc) { for (rw in 1:nr) { m2[rw,nc]= sum( (wgh[rw,nc] <= as.numeric(wgh[rw,1:(nc-1)])),na.rm=T) } }
}

wmaxs <<- apply(m2,1,max)
werrs <<- apply(m2,1,sum)
print(table(werrs,wmaxs))

# which ones to delete
for (i in which(wmaxs>1) ) wgh[i,which(m2[i,]==wmaxs[i])]=NA

# which ones to flag as candidates
wmsk2=matrix(0,nr=nr,nc=nc)  # first version
for (i in which(wmaxs==1) ) wmsk2[i,]=as.numeric(m2[i,]==wmaxs[i])
rm(m2)

wmsk2  <<- wmsk2
a[,weights] <<- wgh   # not sure whether this line works
rm(wgh)

cat("wmsk2 lmsk2 matrixes were created")
cat("wmaxs werrs matrixes were created")
cat("lmaxs lerrs matrixes were created")



} # END OF FUNCTION


