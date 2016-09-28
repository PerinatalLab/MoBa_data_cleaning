

# use this code as a function

# priority is given to the MFR


refillMFRwithQ4=function(a) {

  par(mfrow=c(1,2))
  # WEIGHT data
  wqM=a[,"WeightBirthMFR"]
  wq4=a[,"WeightBirthQ4"]
  
  #hist(wqM,breaks=100,col="grey")
  #hist(wq4,breaks=100,col="grey")
  
  # fix the weight
  wbad=which( (is.na(wqM) | (wqM<(mean(wqM,na.rm=T)-4*sd(wqM,na.rm=T))) |
              (wqM>(mean(wqM,na.rm=T)+4*sd(wqM,na.rm=T)))))
  
  #length(wbad)
  plot(wqM[wbad]~wq4[wbad])
  wqM[wbad]=wq4[wbad]
    
  # LENGTH data
  lqM=a[,"LengthBirthMFR"]
  lq4=a[,"LengthBirthQ4"]
  
  #hist(lqM,breaks=100,col="grey")
  #hist(lq4,breaks=100,col="grey")
  
  # fix the length
  lbad=which( (is.na(lqM) | (lqM<(mean(lqM,na.rm=T)-4*sd(lqM,na.rm=T))) |
                 (lqM>(mean(lqM,na.rm=T)+4*sd(lqM,na.rm=T)))))
  #length(lbad)
  plot(lqM[lbad]~lq4[lbad])
    lqM[lbad]=lq4[lbad]
  
  # SAVE ALL CHANGES THAT WERE MADE INTO A MAIN FILE
  a[,"WeightBirthMFR"] <<- wqM  #  note the global assignment symbol (<<-)
  a[,"LengthBirthMFR"] <<- lqM  #  note the global assignment symbol (<<-)
  
}

# script should also include Head Circumferance values