

TimeBasedWeightLengthOverlaps=function(a,weights,lengths) {
# see overlaps between weight and length time-based flags
par(mfrow=c(3,4))

cat("
      color
RED:   length-based
GREEN: weight-based
BLUE:  both length and weight -based
    ")

  for (i in 1:12) {
  l=a[,lengths[i]]
  w=a[,weights[i]]
  colL=rep(0,dim(a)[1])
  colL[which(lmsk2[,i]==1)]=1
  colW=rep(0,dim(a)[1])
  colW[which(wmsk2[,i]==1)]=2
  COL=colL+colW+1
  print(table(COL))
  daa=data.frame(l=l,w=w,COL=COL)
  daa=daa[order(daa$COL),]
  plot(daa$w~daa$l,col=daa$COL,xlab=lengths[i],ylab=weights[i],main=i)
}

}