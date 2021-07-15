
#TO USE
#update the armmat lines to give the months each arm is in the trial
#you can change the number of months and arms
#Then use
#TrialInfo(armmat,numpermonth,ARMTOESTIMATE)
#    where ARMTOESTIMATE is the arm of interest
#    e.g. TrialInfo(armmat,numpermonth,5) to estimate arm 5

#INPUTS HERE
#arm matrix...rows are arms, columns are interims
#    entry 1=arm in trial at that interim, 0=arm not in trial at that interim
#code currently only support 10 arm, unlimited number of months
#    might be very slow or crash with too many months
armmat=matrix(0,nrow=5,ncol=12)
armmat[1,]=c(1,1,1,0,0,0,0,0,0,0,0,0)
armmat[2,]=c(1,1,1,1,1,1,1,0,0,0,0,0)
armmat[3,]=c(0,0,0,1,1,1,1,1,1,0,0,0)
armmat[4,]=c(0,0,0,0,0,0,0,1,1,1,1,1)
armmat[5,]=c(0,0,0,0,0,0,0,0,0,1,1,1)
#INPUTS FINISHED

#DERIVED QUANTITIES
#    this is a bit of a hack...assume a number
#    of subjects per month that lets every possible
#    number of active arms divide evenly
#    (least common multiple)
#    would be better to derive matrices more cleverly
numpermonth=c(1,2,6,12,60,60,420,840,2520,2520)
numpermonth=numpermonth[nrow(armmat)+1]
#END DERIVED QUANTITIES

ExpandArmMat=function(armmat,numpermonth) {
  #There are numpermonth*(number of months) patients in trial
  #this function creates indicator variables for all of the
  #arms and months for future use
  nummonths=ncol(armmat)
  numarms=nrow(armmat)   #number of active arms
  numpatients=nummonths*numpermonth
  ctrlvec=rep(0,numpatients)
  xarm=matrix(0,nrow=numpatients,ncol=numarms)
  xmonth=matrix(0,nrow=numpatients,ncol=nummonths)
  #tally control and arm indicators
  cur=0
  for (j in 1:nummonths) {
    #calculate number of patients per arm per month
    curarmmonth=numpermonth/(sum(armmat[,j])+1)
    #mark control patients for month j
    ctrlvec[cur+(1:curarmmonth)]=1
    cur=cur+curarmmonth
    #go through active arms
    for (i in 1:numarms) {
      #mark active arm patients in appropriate column
      #only includes patients for arms actively enrolling
      #in that month
      if (armmat[i,j]==1) {
        xarm[cur+(1:curarmmonth),i]=1
        cur=cur+curarmmonth
      }
    }
  }
  #tally month indicators
  cur=0
  for (j in 1:nummonths) {
    xmonth[cur+(1:numpermonth),j]=1
    cur=cur+numpermonth
  }
  #return values
  return(list(ctrlvec=ctrlvec,xarm=xarm,xmonth=xmonth))
}

GetDesignMat=function(armmat,numpermonth) {
  #get indicator variables for control, arms, months
  indlist=ExpandArmMat(armmat,numpermonth)
  #create design matrix
  xmat=matrix(0,nrow=length(indlist$ctrlvec),ncol=nrow(armmat)+ncol(armmat))
  #intercept
  xmat[,1]=1
  #arm indicators
  xmat[,1+(1:nrow(armmat))]=indlist$xarm
  #month indicators (no month 1 for identifiability)
  xmat[,1+nrow(armmat)+(1:(ncol(armmat)-1))]=indlist$xmonth[,-1]
  return(xmat)
}

#Model = Intercept + Arm Effects + Month Effects + error
#Can get design matrix X
#the arm 5 effect is (X^t X)^{-1} X^t y
#can derive out this effect as a weighted average of the means in each arm/month

TrialInfo=function(armmat,numpermonth,armtoestimate=5) {
  #armtoestimate is the arm effect of interest
  #get indicator variables and design matrix
  indlist=ExpandArmMat(armmat,numpermonth)
  xmat=GetDesignMat(armmat,numpermonth)
  numarms=nrow(armmat)
  nummonths=ncol(armmat)
  #We know betahat = (X^t X)^{-1} X^t y
  #    coefmat is the first part of that matrix
  #    if we take a row of that matrix, we get
  #    the coefficients for a weighted combination
  #    of the data that is used to estimate the
  #    treatment effect
  coefmat=solve(t(xmat)%*%xmat)%*%t(xmat)
  coefvec=coefmat[1+armtoestimate,]
  
  #The estimate for armtoestimate is
  #     sum(coefvec*data)
  #     we aren't computing the estimate here (no data!)
  #     but you can see how the estimate is constructed
  #The coefficients are constant for all data from
  #     a particular arm and interim (e.g. all control data
  #     in interim 6). Thus, the treatment effect estimate is
  #sum_{arm a, month m} c_{a,m} ybar_{a,m}
  #     we can find the c_{a,m} coefficients by adding
  #     coefvec over all the observations in arm a, month m
  
  combomat=matrix(0,nrow=1+numarms,ncol=nummonths)
  for (j in 1:nummonths) {
    #ctrl first
    use=which((indlist$ctrl==1)&(indlist$xmonth[,j]==1))
    combomat[1,j]=sum(coefvec[use])
    #all active arms
    for (i in 1:numarms) {
      use=which((indlist$xarm[,i]==1)&(indlist$xmonth[,j]==1))
      combomat[i+1,j]=sum(coefvec[use])
    }
  }
  rownames(combomat)=c("Control","Arm A","Arm B","Arm C","Arm D","Arm E")
  
  print("The following matrix gives the weight given to data in each arm/month")
  print(round(combomat,4))
  
  print("To verify the the time effects drop, note the columns sum to 0")
  print("The desired treatment effect is unbiased because the row sum for")
  print("Control is (-1) and the row sum for the arm of interest is 1")
  print("All other row sums are 0")
  
  print("The variance reduction using all data is")
  varmat=solve(t(xmat)%*%xmat)
  varall=varmat[armtoestimate+1,armtoestimate+1]
  nconcurrentperarm=sum(indlist$xarm[,armtoestimate])
  varconcurrent=2/nconcurrentperarm
  varratio=varall/varconcurrent
  #print(c(varall,varconcurrent,varratio))
  print(round(varratio,3))
  print("of the variance from using only concurrent controls")
}

