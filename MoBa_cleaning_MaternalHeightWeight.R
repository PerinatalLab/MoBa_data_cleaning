# a script for cleaning height and weight variables in MoBa projects
# by Jonas Bacelis
# 2015 09 16

options(stringsAsFactors=F)

# load the phenotypes (height and weight) taken directly from MoBa
#### NOTE THE COLUMNS REQUIRED FOR THE INPUT CSV FILE
q1=read.csv("~/Desktop/MoBa_v6/output_q1a_basicsizeinfo.csv",h=F,sep=",")
names(q1)=c("PREG_ID","AA85","AA86","AA87")
head(q1); dim(q1)

# save a primal state for describing all changes in the end
q1_backup = q1

# load maternal ID and PregID data
#### NOTE THE COLUMNS REQUIRED FOR THE INPUT CSV FILE
svi = read.csv("~/Desktop/MoBa_v5/SV_INFO_converted.csv",h=F,",",stringsAsFactors=F)
colnames(svi) = c("PREG_ID","M_ID","YEAR")
head(svi); dim(svi)

# set output file and folder

## get the date stamp
date_stamp = paste(unlist(strsplit(substr(Sys.time(),1,10),"-")),collapse="")
## get the hash of the current state of the Git folder
hash = system(paste("git log --pretty=format:'%h' -n 1"),intern=TRUE)

file_dir = "~/Desktop/MoBa_v6/"
file_name = paste(file_dir,"MOBA_PDB1581_CLEANED_maternalHgh1Wgh1Wgh2_",date_stamp,"_",hash,".txt",sep="")


# check! there should be no twin pregnancies in Q1 and SVI file
sum(duplicated(q1$PREG_ID))
sum(duplicated(svi$PREG_ID))

# define the function of plotting large quantities of values values
library(hexbin)
fun_plot = function(X,Y,LOG,XLAB,YLAB) {
        # LOG = 100 works best
        h=hexbin(Y~X); x=h@xcm; y=h@ycm; s=h@count
        plot(x,y,pch=1,cex=log(s,LOG),xlab=XLAB,ylab=YLAB,
             xlim=range(X,na.rm=T),ylim=range(Y,na.rm=T))
}


# define a function for fixing unreasonable HEIGHT values (or deleting otherwise)
fun_heightCuration = function(height,low,upp,thr) {
        # low = lowest allowed resonable value
        # upp = highest allowed reasonable value
        # thr  = threshold for normal values (in standard deviations)
        flags = rep("original",length(height)) # vector of flags with a code of curation
        h = hist(height,breaks=50,col="grey")
        tmp = height[ which( (height>low)&(height<upp) ) ]
        # estimate the median
        mdn=median(tmp,na.rm=T)
        # estimate two standard deviations
        tmp1 = tmp-mdn
        lft = c( tmp1[which(tmp1<=0)],abs(tmp1[which(tmp1<0)]) )
        rgh = c( tmp1[which(tmp1>=0)],-(tmp1[which(tmp1>0)]) )
        sdl = sd(lft,na.rm=T)  # sd of the left side
        sdr = sd(rgh,na.rm=T)  # sd of the right side
        
        abline(v=c(mdn-thr*sdl,mdn+thr*sdr))
        # what fraction of individuals falls outside the thresholds?
        print("fraction of samples inside thresholds:")
        print(sum((tmp>(mdn-thr*sdl))&(tmp<(mdn+thr*sdr)),na.rm=T) / sum(!is.na(tmp)))
        
        # saving private Ryan ("1" was lost at the front)
        left_lim = floor(mdn - thr*sdl - 100)
        rght_lim = floor(mdn + thr*sdr - 100)+1
        if(rght_lim >= 100) rght_lim = 99
        segments(y0=max(h$counts)/2, y1=max(h$counts)/2, 
                 x0=left_lim,x1=rght_lim,col="turquoise3",lwd=10)
        bad1=which( (height>=left_lim)&(height<=rght_lim) )
        
        # saving lost units ("end" was lost)
        left_lim = floor((mdn-thr*sdl)/10)
        rght_lim = floor((mdn+thr*sdr)/10 )
        segments(y0=max(h$counts)/2, y1=max(h$counts)/2, 
                 x0= left_lim,x1=rght_lim,col="magenta",lwd=10)
        bad2=which((height>=left_lim)&(height<=rght_lim)) # note >= and <=
        
        
        points(height[bad1],rep(max(rep(h$counts))*.4,length(bad1)))
        points(height[bad2],rep(max(rep(h$counts))*.4,length(bad2)))
        # fix what is fixable
        height[bad1]=height[bad1]+100
        height[bad2]=height[bad2]*10+5 ## no reason to add exactly 5
        flags[bad1]="changed" # best-guess curation
        flags[bad2]="changed" # best-guess curation
        
        # eliminate outliers
        bad3 = which( (height<low)|(height>upp) )
        points(height[bad3],rep(max(rep(h$counts))*.3,length(bad3)),col="red")
        height[which(height<low)]=NA  # ***
        flags[bad3]="deleted" # best-guess curation
        out = data.frame(height,flag = flags,stringsAsFactors = F)
        out
}  # end of function




# fix unreasonable WEIGHT values or clean otherwise
fun_weightCuration = function(weight,low,upp,thr,xlim) {
        # low = lowest allowed resonable value
        # upp = highest allowed reasonable value
        # thr  = threshold in standard deviations
        flags = rep("original",length(weight)) # vector of flags with a code of curation
        h = hist(weight,breaks=200,col="grey")
        if (sum(is.na(xlim))==0) h = hist(weight,breaks=200,col="grey",xlim=xlim) # better resolution
        tmp = weight[ which( (weight>low)&(weight<upp) ) ]
        # estimate the median
        mdn=median(tmp,na.rm=T)
        # estimate two standard deviations
        tmp1 = tmp-mdn
        lft = c( tmp1[which(tmp1<=0)],abs(tmp1[which(tmp1<0)]) )
        rgh = c( tmp1[which(tmp1>=0)],-(tmp1[which(tmp1>0)]) )
        sdl = sd(lft,na.rm=T)  # sd of the left side
        sdr = sd(rgh,na.rm=T)  # sd of the right side
        
        abline(v=c(mdn-thr*sdl,mdn+thr*sdr))
        # what fraction of individuals falls outside the thresholds?
        print("fraction of samples inside thresholds:")
        print(sum((tmp>(mdn-thr*sdl))&(tmp<(mdn+thr*sdr)),na.rm=T) / sum(!is.na(tmp)))
        
        # lost last digit
        left_lim = floor( (mdn - thr*sdl)/10)
        rght_lim = floor( (mdn + thr*sdr)/10) +1
        #sort(weight)[1:100]
        segments(y0=max(h$counts)/2, y1=max(h$counts)/2, 
                 x0=left_lim,x1=rght_lim,col="turquoise3",lwd=10)
        bad1=which( (weight>=left_lim)&(weight<=rght_lim) )
        
        # lost decimal point
        left_lim = floor((mdn-thr*sdl)*10)
        rght_lim = floor((mdn+thr*sdr)*10 )
        segments(y0=max(h$counts)/2, y1=max(h$counts)/2, 
                 x0= left_lim,x1=rght_lim,col="magenta",lwd=10)
        bad2=which((weight>=left_lim)&(weight<=rght_lim)) # note >= and <=
        
        points(weight[bad1],rep(max(rep(h$counts))*.4,length(bad1)))
        points(weight[bad2],rep(max(rep(h$counts))*.4,length(bad2)))
        # fix what is fixable
        weight[bad1]=weight[bad1]*10+5 ### no reason to add +5, but whatever
        weight[bad2]=weight[bad2]/10
        flags[bad1]="changed" # best-guess curation
        flags[bad2]="changed" # best-guess curation
        
        # eliminate outliers
        bad3 = which( (weight<low)|(weight>upp) )
        points(weight[bad3],rep(max(rep(h$counts))*.3,length(bad3)),col="red")
        weight[which(weight<low)]=NA  # ***
        flags[bad3]="deleted" # best-guess curation
        
        #table(flags,useNA="a")
        out = data.frame(weight,flag = flags,stringsAsFactors = F)
        out
} # end of function



###  comparing variable vs similar variable
fun_weightVSweight = function(weight1,weight2,upp,low) {
        # upp = 30        # low = -30
        if(length(weight1)!=length(weight2)) warning("lengths don't match")
        flags = rep("original",length(weight1)) # vector of flags with a code of curation
        plot(weight1~weight2)
        abline(upp,1,col="red"); abline(low,1,col="red"); abline(0,1,col="blue")
        flags[which(weight1 > (upp+weight2))] = "deleted"
        flags[which(weight1 < (low+weight2))] = "deleted"
        
        weight1[which(flags=="deleted")] = NA
        weight2[which(flags=="deleted")] = NA
        lst = data.frame(weight1,weight2,flag = flags,stringsAsFactors = F)
        lst
}



###  comparing variable vs diferent variable
fun_weightVSheight = function(weight,height,WthrL,WthrR) {
        if(length(weight)!=length(height)) warning("lengths do not match")
        fun_plot(height,weight,100,"height","weight") 
        col=NULL
        for (hgh in sort(unique(height))) {
                wgh = weight[which(height %in% (hgh-2):(hgh+2))]
                hhh = height[which(height %in% (hgh-2):(hgh+2))]
                if(sum(!is.na(wgh))>20) {
                        d = density(wgh,na.rm=T)
                        mod = d$x[which(d$y==max(d$y))]
                        wgh0 = wgh-mod
                        wghL = c(wgh0[which(wgh0<0)], abs(wgh0[which(wgh0<=0)]))
                        wghR = c(wgh0[which(wgh0>0)], -wgh0[which(wgh0>=0)])
                        sdL = sd(wghL,na.rm=T)
                        sdR = sd(wghR,na.rm=T)
                        
                        ixL = which(wgh < (mod-sdL*WthrL))
                        ixR = which(wgh > (mod+sdR*WthrR))
                        ixs = unique(c(ixL,ixR))
                        
                        xs  = hhh[ixs]
                        ys  = wgh[ixs]
                        points(xs,ys,pch=19,col="red",cex=0.8)
                        
                        col = c(col, paste(xs,ys,sep="_")) # height_weight
                        points(xs,ys,pch=19,col="red",cex=0.8)
                        rm(xs,ys,ixL,ixR,ixs,sdL,sdR,wgh0,wghL,wghR,mod,d)
                } else {
                        if(sum(!is.na(wgh))>2) { 
                                sdA = sd(wgh,na.rm=T)
                                mnA = mean(wgh,na.rm=T)
                                bad = unique(c(which(wgh < (mnA-sdA*2) ), which(wgh > (mnA+sdA*2)) ))
                                xs = hhh[bad]
                                ys = wgh[bad]
                                points(xs,ys,pch=19,col="blue",cex=1)
                                col = c(col, paste(xs,ys,sep="_")) # height_weight
                                rm(xs,ys,sdA,mnA,bad)
                                
                        }
                }
        }
        
        col=unique(col); length(col)
        ind = paste(height,weight,sep="_")
        flags = rep("original",length(weight))
        flags[which(ind %in% col)] = "deleted"
        out = data.frame(weight=weight,height=height,flag=flags,stringsAsFactors = F)
        out$weight[which(out$flag=="deleted")]=NA
        out$height[which(out$flag=="deleted")]=NA
        out
} # end of function



# define the function that checks whether bmi values are reasonable
fun_bmi = function(weight,height,thrL,thrR) {
        # thrL = threshold in SDs for low values, thrR = threshold in SDs for high values
        if(length(height)!=length(weight)) warning("Height and Weight lengths do not match!")
        bmi = weight / ((height/100)^2)
        flags = rep("original",length(weight))
        flags[which((bmi>100)|(bmi<10))] = "deleted" # sanity cleaning
        
        hist(bmi,breaks=100,col="grey",xlim=c(0,100))
        d = density(bmi,na.rm=T)
        mod = d$x[which(d$y==max(d$y))]
        bmi0 = bmi-mod
        bmiL = c(bmi0[which(bmi0<0)],  abs(bmi0[which(bmi0<=0)]))
        bmiR = c(bmi0[which(bmi0>0)],-(bmi0[which(bmi0>=0)]))
        sdL = sd(bmiL,na.rm=T)
        sdR = sd(bmiR,na.rm=T)
        abline(v=mod,col="red",lwd=4)
        abline(v=c(mod-thrL*sdL,mod+thrR*sdR),col="blue",lwd=2)
        abline(v=c(mod-sdL,mod+sdR),col="blue",lwd=1)
        
        flags[which(bmi<(mod-thrL*sdL))]="deleted"
        flags[which(bmi>(mod+thrR*sdR))]="deleted"
        fun_plot(height,weight,100,"height","weight") 
        points(weight[which(flags=="deleted")]~height[which(flags=="deleted")],
               pch=19,cex=0.7,col="red")
        
        out = data.frame(weight=weight,height=height,flag=flags,stringsAsFactors = F)
        out$weight[which(out$flags=="deleted")]=NA
        out$height[which(out$flags=="deleted")]=NA
        out
} # end of function



########################
########################

### first round of cleaning: CURATION of hand-writing problems and OCR problems
# play with constants, and when you are satisfied - proceed further
w1 = fun_weightCuration(weight=q1$AA85,low=30,upp=200,thr=3,xlim=NA)
w2 = fun_weightCuration(weight=q1$AA86,low=30,upp=200,thr=3,xlim=NA)
h1 = fun_heightCuration(q1$AA87,130,210,3)

# summary of changes
table(w1$flag)
table(w2$flag)
table(h1$flag)

# assign new (corrected) values
q1$AA85 = w1$weight
q1$AA86 = w2$weight
q1$AA87 = h1$height

########################
########################
########################

# second round of cleaning:  various deletions
# first, play with constants, view summaries, and when happy - proceed
w1w2_1 = fun_weightVSweight(q1$AA85,q1$AA86,30,-30); table(w1w2_1$flag)
w1h1_1 = fun_weightVSheight(q1$AA85,q1$AA87,4,5); table(w1h1_1$flag)
w2h1_1 = fun_weightVSheight(q1$AA86,q1$AA87,4,5); table(w2h1_1$flag)
bm1_1 = fun_bmi(q1$AA85,q1$AA87,thrL=4,thrR=5); table(bm1_1$flag)
bm2_1 = fun_bmi(q1$AA86,q1$AA87,thrL=4,thrR=5); table(bm2_1$flag)        

# preview selected deletions in one plot:
fun_plot(q1$AA87,q1$AA85,100,"height","weight")
points(AA85~AA87,data = q1[which(bm1_1$flag=="deleted"),],pch=19,col="green",cex=1.2)
points(AA85~AA87,data = q1[which(bm2_1$flag=="deleted"),],pch=3,col="black")
points(AA85~AA87,data = q1[which(w1w2_1$flag=="deleted"),],pch=19,col="blue",cex=0.8)
points(AA85~AA87,data = q1[which(w1h1_1$flag=="deleted"),],pch=19,col="red",cex=0.6)
points(AA85~AA87,data = q1[which(w2h1_1$flag=="deleted"),],pch=1,col="black",cex=2)

fun_plot(q1$AA87,q1$AA86,100,"height","weight (preg)")
points(AA86~AA87,data = q1[which(bm1_1$flag=="deleted"),],pch=19,col="green",cex=1.2)
points(AA86~AA87,data = q1[which(bm2_1$flag=="deleted"),],pch=3,col="black")
points(AA86~AA87,data = q1[which(w1w2_1$flag=="deleted"),],pch=19,col="blue",cex=0.8)
points(AA86~AA87,data = q1[which(w2h1_1$flag=="deleted"),],pch=19,col="red",cex = 0.6)
points(AA86~AA87,data = q1[which(w1h1_1$flag=="deleted"),],pch=1,col="black",cex=2)


# delete selected ("untrustworthy") values
q1$AA85[which(w1w2_1$flag=="deleted")]=NA
q1$AA85[which(w1h1_1$flag=="deleted")]=NA
q1$AA85[which(bm1_1$flag=="deleted")]=NA
q1$AA86[which(w1w2_1$flag=="deleted")]=NA
q1$AA86[which(w2h1_1$flag=="deleted")]=NA
q1$AA86[which(bm2_1$flag=="deleted")]=NA
q1$AA87[which(w1h1_1$flag=="deleted")]=NA
q1$AA87[which(w2h1_1$flag=="deleted")]=NA
q1$AA87[which(bm1_1$flag=="deleted")]=NA
q1$AA87[which(bm2_1$flag=="deleted")]=NA


########################
########################
########################
library(sqldf)
q1m=sqldf("SELECT q1.PREG_ID as PREG_ID, M_ID, AA85, AA86, AA87 FROM q1 INNER JOIN svi ON q1.PREG_ID=svi.PREG_ID")

## this stage will remove the smallest (largest) value of measurements per each mother,
## if it differs from the second smallest (largest) value for the same mother
## by an amount greater than the thresholds set here
e85=20
e86=20
e87=5

## select second largest value for every mother
secondMax85=sqldf("SELECT M_ID, max(AA85) AS m FROM q1m
        WHERE PREG_ID NOT IN ( SELECT PREG_ID FROM q1m
                INNER JOIN (
                        SELECT M_ID, max(AA85) as m FROM q1m GROUP BY M_ID
                ) t2 
                ON q1m.AA85=t2.m AND q1m.M_ID=t2.M_ID
          ) GROUP BY M_ID")
secondMax86=sqldf("SELECT M_ID, max(AA86) AS m FROM q1m
        WHERE PREG_ID NOT IN ( SELECT PREG_ID FROM q1m
                  INNER JOIN (
                  SELECT M_ID, max(AA86) as m FROM q1m GROUP BY M_ID
                  ) t2 
                  ON q1m.AA86=t2.m AND q1m.M_ID=t2.M_ID
        ) GROUP BY M_ID")
secondMax87=sqldf("SELECT M_ID, max(AA87) AS m FROM q1m
                  WHERE PREG_ID NOT IN ( SELECT PREG_ID FROM q1m
                  INNER JOIN (
                  SELECT M_ID, max(AA87) as m FROM q1m GROUP BY M_ID
                  ) t2 
                  ON q1m.AA87=t2.m AND q1m.M_ID=t2.M_ID
          ) GROUP BY M_ID")

## select second smallest value for every mother
secondMin85=sqldf("SELECT M_ID, min(AA85) AS m FROM q1m
        WHERE PREG_ID NOT IN ( SELECT PREG_ID FROM q1m
                INNER JOIN (
                        SELECT M_ID, min(AA85) as m FROM q1m GROUP BY M_ID
                ) t2 
                ON q1m.AA85=t2.m AND q1m.M_ID=t2.M_ID
          ) GROUP BY M_ID")
secondMin86=sqldf("SELECT M_ID, min(AA86) AS m FROM q1m
        WHERE PREG_ID NOT IN ( SELECT PREG_ID FROM q1m
                  INNER JOIN (
                  SELECT M_ID, min(AA86) as m FROM q1m GROUP BY M_ID
                  ) t2 
                  ON q1m.AA86=t2.m AND q1m.M_ID=t2.M_ID
        ) GROUP BY M_ID")
secondMin87=sqldf("SELECT M_ID, min(AA87) AS m FROM q1m
                  WHERE PREG_ID NOT IN ( SELECT PREG_ID FROM q1m
                  INNER JOIN (
                  SELECT M_ID, min(AA87) as m FROM q1m GROUP BY M_ID
                  ) t2 
                  ON q1m.AA87=t2.m AND q1m.M_ID=t2.M_ID
          ) GROUP BY M_ID")

## retrieves pregIDs where value is above second largest + a threshold
## or where value is below second smallest - a threshold
remove85=sqldf(paste("SELECT q1m.PREG_ID FROM q1m
        INNER JOIN secondMin85 min ON q1m.M_ID=min.M_ID
        INNER JOIN secondMax85 max ON q1m.M_ID=max.M_ID
        WHERE AA85 > max.m +",e85,"OR AA85 < min.m -",e85))
remove86=sqldf(paste("SELECT q1m.PREG_ID FROM q1m
        INNER JOIN secondMin86 min ON q1m.M_ID=min.M_ID
        INNER JOIN secondMax86 max ON q1m.M_ID=max.M_ID
        WHERE AA86 > max.m +",e86,"OR AA86 < min.m -",e86))
remove87=sqldf(paste("SELECT q1m.PREG_ID FROM q1m
        INNER JOIN secondMin87 min ON q1m.M_ID=min.M_ID
        INNER JOIN secondMax87 max ON q1m.M_ID=max.M_ID
        WHERE AA87 > max.m +",e87,"OR AA87 < min.m -",e87))

### this query can generate a table useful for reports
#sqldf(paste("SELECT q1m.*, min.m as SecondMin, max.m as SecondMax FROM q1m
#        INNER JOIN secondMin87 min ON q1m.M_ID=min.M_ID
#        INNER JOIN secondMax87 max ON q1m.M_ID=max.M_ID
#        WHERE AA87 > max.m +",e87,"OR AA87 < min.m - ",e87))

## set the detected outliers as missing
q1m[which(q1m$PREG_ID %in% remove85$PREG_ID),"AA85"]=NA
q1m[which(q1m$PREG_ID %in% remove86$PREG_ID),"AA86"]=NA
q1m[which(q1m$PREG_ID %in% remove87$PREG_ID),"AA87"]=NA


## save the flags in a data frame
mpflags=data.frame(PREG_ID=q1m$PREG_ID, AA85=c(q1m$PREG_ID %in% remove85$PREG_ID),
                   AA86=c(q1m$PREG_ID %in% remove86$PREG_ID), AA87=c(q1m$PREG_ID %in% remove87$PREG_ID))
mpflags[mpflags==FALSE]="original"
mpflags[mpflags==TRUE]="deleted"

## attach per-mother averages
imp=sqldf("SELECT * FROM q1m
        INNER JOIN (SELECT M_ID, AVG(AA85) as a85, AVG(AA86) as a86, AVG(AA87) as a87 FROM q1m GROUP BY M_ID) as avg
        ON q1m.M_ID=avg.M_ID")
imp$mp=duplicated(imp$M_ID)

## select mothers with parity=2
mps=sqldf("SELECT M_ID FROM ( SELECT M_ID,count(M_ID) as c FROM q1m GROUP BY M_ID ) as counts WHERE counts.c=2")
mps=imp[which(imp$M_ID %in% mps$M_ID),]
mps=mps[order(mps$M_ID),]

## calculate average difference between parities
mpst=cbind(mps[seq(1,nrow(mps),2),"M_ID"], mps[seq(1,nrow(mps),2),"AA85"], mps[seq(2,nrow(mps),2),"AA85"])
mpst=mpst[which(!is.na(mpst[,2]) & !is.na(mpst[,3])),]
d85=mean(as.numeric(mpst[,2])-as.numeric(mpst[,3]))
mpst=cbind(mps[seq(1,nrow(mps),2),"M_ID"], mps[seq(1,nrow(mps),2),"AA86"], mps[seq(2,nrow(mps),2),"AA86"])
mpst=mpst[which(!is.na(mpst[,2]) & !is.na(mpst[,3])),]
d86=mean(as.numeric(mpst[,2])-as.numeric(mpst[,3]))
mpst=cbind(mps[seq(1,nrow(mps),2),"M_ID"], mps[seq(1,nrow(mps),2),"AA87"], mps[seq(2,nrow(mps),2),"AA87"])
mpst=mpst[which(!is.na(mpst[,2]) & !is.na(mpst[,3])),]
d87=mean(as.numeric(mpst[,2])-as.numeric(mpst[,3]))

## weight increases by d85 (d86) between pregnancies, which should be taken into account while imputing
## height increase is negligible
## weight increase between 2nd and 3rd pregnancies is also observed, but there's few of them

imp$AA85[is.na(imp$AA85)]=round(imp$a85[is.na(imp$AA85)]+d85*(-1)^imp$mp[is.na(imp$AA85)],1)
imp$AA86[is.na(imp$AA86)]=round(imp$a86[is.na(imp$AA86)]+d86*(-1)^imp$mp[is.na(imp$AA86)],1)
imp$AA87[is.na(imp$AA87)]=round(imp$a87[is.na(imp$AA87)],0)

## successful imputes:
i85=which(imp$AA85 != q1$AA85)
i86=which(imp$AA86 != q1$AA86)
i87=which(imp$AA87 != q1$AA87)

mpflags[i85,"AA85"]="imputed"
mpflags[i86,"AA86"]="imputed"
mpflags[i87,"AA87"]="imputed"

# inspect the flags
table(mpflags$AA85)
table(mpflags$AA86)
table(mpflags$AA87)

q1=imp[,c("PREG_ID","AA85","AA86","AA87")]

############################################ 

###  create a table of all flags 
out = data.frame(bm1_1=(bm1_1$flag=="deleted"),bm2_1=(bm2_1$flag=="deleted"),
                 w1w2_1=(w1w2_1$flag=="deleted"),w1h1_1=(w1h1_1$flag=="deleted"),
                 w2h1_1=(w2h1_1$flag=="deleted"), stringsAsFactors=F)
head(out); dim(out)

# vectors that decide what should be deleted
fun_sum = function(x) sum(x,na.rm=T)
h1t = apply(out[,c(1,2,4,5)],1,fun_sum) # for height
w1t = apply(out[,c(1,3,4)],1,fun_sum) # for prepreg weight
w2t = apply(out[,c(2,3,5)],1,fun_sum) # for pregnancy weight
table(h1t); table(w1t); table(w2t)

###  delete what is selected for deletion
q1$AA85[which(w1t>0)]=NA
q1$AA86[which(w2t>0)]=NA
q1$AA87[which(h1t>0)]=NA

### (note: multiparity-based outliers were already deleted earlier)

# check how many are missing now
sum(is.na(q1$AA85))
sum(is.na(q1$AA86))
sum(is.na(q1$AA87))

# how many pregnancies have (and how many do not) a full set of data
print("missing data:")
sum( (is.na(q1$AA85))|(is.na(q1$AA87)))
print("full data:")
sum( (!is.na(q1$AA85))&(!is.na(q1$AA87)))

# create a dataframe of flags describing what happened with each Pregnancy
## 0 - original
## 1 - corrected
## 2 - imputed from other pregnancies
## 3 - deleted
flags = data.frame(PREG_ID=q1$PREG_ID,AA85=NA,AA86=NA,AA87=NA)
for (col_name in c("AA85","AA86","AA87")) {
        changes = as.numeric(q1[,col_name] != q1_backup[,col_name])
        changes[is.na(changes)] = 0
        deletions = (is.na(q1[,col_name]))&(!is.na(q1_backup[,col_name]))
        # table(changes); table(deletions)
        #since changes were done first and only then -  deletions :
        changes[which(deletions)]=3
        imputeds=(mpflags[,col_name]=="imputed")
        
        # produce an error if a row was set as imputed, but not as a change
        changes[imputeds & (changes!=1)]="ERROR"
        changes[imputeds & (changes==1)]=2
        flags[,col_name] = changes
        rm(changes,deletions)
}

# preview
head(flags)
table(flags$AA85)
table(flags$AA86)
table(flags$AA87)


# rename columns so that they would not be confused with original data columns
ex_names = colnames(flags)[grep("^AA",colnames(flags))]
new_names = paste("fl",ex_names,sep="")
colnames(flags)[grep("^AA",colnames(flags))] = new_names
head(flags)

# create the object that will be saved
output_obj = data.frame(q1,flags[,-grep("PREG_ID",colnames(flags))])
head(output_obj)

# write output
write.table(output_obj,file_name,row.names=F,col.names=T,sep="\t",quote=F)
