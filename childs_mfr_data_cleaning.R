#!/usr/bin/Rscript

# WARNING! this is a very quick-and-dirty way of cleaning mfr data

#########  CLEAN DATA
# clean weight values
hist(dat$AA85,breaks=100,col="grey")
dat$AA85[which(dat$AA85<37)] = NA  # too small
dat$AA85[which(dat$AA85>370)] = round(dat$AA85[which(dat$AA85>370)]/10,0) # decimal mistake

# clean height values
hist(dat$AA87,breaks=100,col="grey")
dat$AA87[which(dat$AA87<50)] = NA  # too small
dat$AA87[which(dat$AA87 %in% 50:90)] = 100 + dat$AA87[which(dat$AA87 %in% 50:90)] # one ("1") meter was lost
dat$AA87[which(dat$AA87<140)] = NA  # too small

# clean HeadCircumf and BirthLength values
hist(dat$HODE,breaks=100,col="grey")
hist(dat$LENGDE,breaks=100,col="grey")
plot(dat$LENGDE ~ dat$HODE); abline(6,1.6,col="red"); abline(-11,1.4,col="red") # outliers
ix = which((dat$LENGDE> (dat$HODE*1.6 + 6)) |(dat$LENGDE< (dat$HODE*1.4 - 11)))
dat$LENGDE[ix] = NA
dat$HODE[ix] = NA
dat$HODE[which(dat$HODE<5)] = NA
dat$LENGDE[which(dat$LENGDE<10)] = NA

# clean HeadCircumf and BirthWeight values
plot(dat$VEKT ~ dat$HODE)
# all OK

# clean BirthWeight and BirthLength values
hist(dat$VEKT,breaks=100,col="grey")
plot(dat$VEKT ~ dat$LENGDE); abline(-9000,300,col="red"); abline(-2000,70,col="blue"); abline(h=1200)
ix = which(((dat$VEKT>1200) & (dat$VEKT > (dat$LENGDE*300 - 9000))) | (dat$VEKT < (dat$LENGDE*70 - 2000)))
dat$LENGDE[ix] = NA
dat$VEKT[ix] = NA

# Placent weight based on BirthWeight
hist(dat$PLACENTAVEKT,breaks=100,col="grey") # ,xlim=c(0,2000)
ix = sample(nrow(dat),1e4,replace=F); plot(dat$PLACENTAVEKT[ix],dat$VEKT[ix]); abline(0,3.2)  # two clusters: twins!
bad_ix = which((dat$VEKT<dat$PLACENTAVEKT*3.2)&(dat$PLACENTAVEKT>499))
table(dat$FLERFODSEL[bad_ix])
dat$PLACENTAVEKT[bad_ix] = NA

# clean HeadCircumf based on GestAge
plot(dat$HODE ~ dat$SVLEN_DG); abline(-5,0.2,col="red"); abline(0,0.1,col="red")
ix = which((dat$HODE> (dat$SVLEN_DG*0.2 - 5)) |(dat$HODE< (dat$SVLEN_DG*0.1 - 1)))
dat$HODE[ix] = NA

# clean BirthLength based on GestAge
plot(dat$LENGDE ~ dat$SVLEN_DG); abline(8,0.2,col="red"); abline(-5,0.15,col="red")
ix = which((dat$LENGDE>(dat$SVLEN_DG*.2+8))|((dat$LENGDE<(dat$SVLEN_DG*.15-5))&(dat$SVLEN_DG>150)))
points(dat$SVLEN_DG[ix],dat$LENGDE[ix],pch=19,col="red")
dat$LENGDE[ix] = NA

# clean Birthweight based on GestAge
plot(dat$VEKT ~ dat$SVLEN_DG); abline(h=1500); abline(-12000,70,col="red"); abline(-4000,20,col="red")
ix = which(((dat$VEKT>1500)&(dat$VEKT>(dat$SVLEN_DG*70-12e3)))|(dat$VEKT<(dat$SVLEN_DG*20-4e3)))
points(dat$SVLEN_DG[ix],dat$VEKT[ix],pch=19,col="red")
dat$VEKT[ix] = NA
