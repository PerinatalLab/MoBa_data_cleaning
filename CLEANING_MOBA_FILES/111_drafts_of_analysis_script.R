

use=which((!is.na(a[,lengths[1]]))&(!is.na(a[,lengths[2]]))&(!is.na(a[,lengths[3]]))&
            (!is.na(a[,lengths[4]]))&(!is.na(a[,lengths[5]]))&(!is.na(a[,lengths[6]]))&
            (!is.na(a[,lengths[7]]))&(!is.na(a[,lengths[8]]))&(!is.na(a[,lengths[9]])))
      #&(!is.na(a[,lengths[10]]))&(!is.na(a[,lengths[11]]))&(!is.na(a[,lengths[12]])))

use=which((!is.na(a[,lengths[7]]))&(!is.na(a[,lengths[8]]))&(!is.na(a[,lengths[9]]))&
          (!is.na(a[,lengths[10]]))&(!is.na(a[,lengths[12]])))

l=a[use,lengths[c(7:10,12)]]
l=a[use,lengths[1:9]]

# preview all timepoints
par(mfrow=c(1,1))
plot(as.numeric(l[1,])~seq(dim(l)[2]),type="l",ylim=c(40,150))
for (i in 2:100) {
lines(as.numeric(l[i,])~seq(dim(l)[2]),lwd=0.1)
}


dis=matrix(NA,nr=dim(l)[1],nc=dim(l)[2]-1)
for (i in 1:dim(l)[1]) {
  for (j in 1:(dim(l)[2]-1)) {
  dis[i,j]=l[i,j+1]-l[i,j]
}
}

head(dis)
fu=function(x) which( x==max(x) )
wer=as.numeric(unlist(apply(dis,1,fu)))
table(wer)
plot(table(wer))
hist(di,breaks=100,col="grey")


# select non-deviants
ll=l[ which(l[,1] %in% 50:55) ,]
lll=ll[which(lll[,2] >= lll[,1]),]
dim(lll)
par(mfrow=c(1,1))
plot(as.numeric(lll[1,1:2])~seq(2),type="l",ylim=c(40,60))
for (i in 2:6624) {
  lines(as.numeric(lll[i,1:2])~seq(2))
}


##  split into two groups based on caffeine consumption:

head(lll)
j1=which(lll$kof_snitt<60); length(j1)
j2=which(lll$kof_snitt>60); length(j2)

par(mfrow=c(1,1))
plot(as.numeric(lll[1,1:10])~seq(10),type="l",ylim=c(50,120))
for (i in j1) lines(as.numeric(lll[i,1:10])~seq(10),col="green")
for (i in j2) lines(as.numeric(lll[i,1:10])~seq(10),col="red")
s1=as.numeric(apply(lll[j1,1:10],2,mean))
s2=as.numeric(apply(lll[j2,1:10],2,mean))
plot(s1~s2,type="l")
points(s1~s2)
abline(0,1)

lll=a[,vars]
j1=which(lll$kof_snitt<60); length(j1)
j2=which(lll$kof_snitt>60); length(j2)
par(mfrow=c(3,3))
for (dd in 1:9) {
di=lll[,dd+1]-lll[,dd]
#hist(di,breaks=100,col="grey")
n=min(length(j1),length(j2))
no=di[j1[1:n]]; no=no[order(no)]
ye=di[j2[1:n]]; ye=ye[order(ye)]
plot(no,ye,xlab="nondrinkers",ylab="drinkers")
abline(0,1)
}

# explor eeffect of the first increase in length
di=lll[,2]-lll[,1]
hist(di[lll==52])
hist(di)
j1=which(di>median(di)+1.5*sd(di))
j2=which(di<median(di)-1.5*sd(di))
par(mfrow=c(1,1))
plot(as.numeric(lll[1,1:10])~seq(10),type="l",ylim=c(40,120))
for (i in j1) lines(as.numeric(lll[i,1:10])~seq(10),col="green")
for (i in j2) lines(as.numeric(lll[i,1:10])~seq(10),col="red")
s1=as.numeric(apply(lll[j1,1:10],2,median))
s2=as.numeric(apply(lll[j2,1:10],2,median))
plot(s1~s2,type="l")
points(s1~s2)
abline(0,1)

colnames(lll)


# explore effect of gestatinal age
s1=as.numeric(apply(ll[which(ll$GestationalAge<260),1:10],2,mean))
s2=as.numeric(apply(ll[which(ll$GestationalAge>=260),1:10],2,mean))
plot(s1~s2,type="l")
points(s1~s2)
abline(0,1)
