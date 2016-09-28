

# use this code as a function

# priority is given to he Q5 since this one was filled closer-in-time to the
# age in question


#weights=c("WeightBirthQ4","Weight6w","Weight3m","Weight5_6m","Weight8m",
#          "Weight1y","Weight15_18mQ5","Weight15_18mQ6","Weight2y",
#          "Weight3y","Weight5y","Weight7y","Weight8y","WeightBirthMFR")


#lengths=c( "LengthBirthQ4","Length6w","Length3m","Length5_6m","Length8m",
#           "Length1y","Length15_18mQ5","Length15_18mQ6","Length2y",
#           "Length3y","Length5y","Length7y","Length8y","LengthBirthMFR")


refillQ5withQ6=function(a) {


# WEIGHT data
wq5=a[,"Weight15_18mQ5"]
wq6=a[,"Weight15_18mQ6"]
#hist(wq5,breaks=100,col="grey")
#hist(wq6,breaks=100,col="grey")

# fix the weight
#hist(wq5,breaks=100,col="grey",ylim=c(0,20))
wbad=which((wq5<4000)|(wq5>=20000)|(is.na(wq5)))
#hist(wq5[wbad],breaks=100,col="grey")
wq5[wbad]=wq6[wbad]
#hist(wq5[wbad],breaks=100,col="grey")
#hist(wq5,breaks=100,col="grey")


# LENGTH data
lq5=a[,"Length15_18mQ5"]
lq6=a[,"Length15_18mQ6"]
#hist(lq5,breaks=100,col="grey")
#hist(lq6,breaks=100,col="grey")

# fix the length
#hist(lq5,breaks=100,col="grey",ylim=c(0,20))
lbad=which((lq5<40)|(is.na(lq5)))
#length(lbad)
#data.frame(lq5=lq5[lbad],lq6=lq6[lbad])
lq5[lbad]=lq6[lbad]

# SAVE ALL CHANGES THAT WERE MADE INTO A MAIN FILE
a[,"Weight15_18mQ5"] <<- wq5  #  note the global assignment symbol (<<-)
a[,"Length15_18mQ5"] <<- lq5  #  note the global assignment symbol (<<-)

}

# SANDBOX
#wthrlw=mean(wq5[which((wq5>5e3)&(wq5<17e3))],na.rm=T)-6*sd(wq5[which((wq5>5e3)&(wq5<17e3))],na.rm=T)
#wthrup=mean(wq5[which((wq5>5e3)&(wq5<17e3))],na.rm=T)+6*sd(wq5[which((wq5>5e3)&(wq5<17e3))],na.rm=T)
#hist(wq5,breaks=100,col="grey"); abline(v=c(wthrlw,wthrup))
#wbad=which((wq5<wthrlw)|(wq5>wthrup))

#lthrlw=mean(lq5[which((lq5>60)&(lq5<90))],na.rm=T)-8*sd(lq5[which((lq5>60)&(lq5<90))],na.rm=T)
#lthrup=mean(lq5[which((lq5>60)&(lq5<90))],na.rm=T)+8*sd(lq5[which((lq5>60)&(lq5<90))],na.rm=T)
#hist(lq5,breaks=100,col="grey"); abline(v=c(lthrlw,lthrup))
