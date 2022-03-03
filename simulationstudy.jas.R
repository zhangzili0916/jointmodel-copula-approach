library("copula")
library("evd")
library(CVTuningCov)
library(mvtnorm)
library("abind")
library("magic")
library("matrixcalc")
library(tikzDevice)
library(joineR)

#Simulation study
#For up to m longitudinal measurements generate data 
Y.lmmsurGaucopulaNOPLOT=function(beta1,beta2,s,n,r,sigma,rhoty,rhoy,rate,trprob,genprob,ageprob,lmmcor,surlmmcor)
{
  m=length(s)                                                                   #s is the prespecified measurement time
  Data=NULL
  Dataobs=NULL
  R11=1
  if(lmmcor=="ex")  R21=rep(rhoty,m)                                            #lmmcor is correlation (Ry) within the longitudinal process 
  if(lmmcor=="ar")  R21=rhoty^(m:1)
  if(surlmmcor=="ex") R22=matrix(c(rep(c(1,rep(rhoy,m)),(m-1)),1),ncol=m)       #surlmmcor is correlation (Ry) between the two processes 
  if(surlmmcor=="ar")  R22=AR1(m,rhoy)
  R12=t(R21)
  R=rbind(cbind(R11,R12),cbind(R21,R22))
  Zi1=matrix(c(rep(1,m),s),ncol=2)
  for(subj in 1:n)
  {
    treat=rbinom(1,1,trprob)
    gender=rbinom(1,1,genprob)
    age=sample(c(0,1,2),1,prob=ageprob)
    Xi1=matrix(c(c(Zi1),rep(treat,m)*s,rep(gender,m),rep(as.numeric(age==1),m)
                 ,rep(as.numeric(age==2),m)),ncol=6)
    Xi2=matrix(c(1,treat,gender,as.numeric(age==1),as.numeric(age==2)),ncol=5)
    ceni=rexp(1,rate)
    Yi=Xi1%*%beta1
    scale=exp(-Xi2%*%beta2/r)
    param=list(list(shape=r,scale=scale))
    for(i in 2:(m+1))
    {
      param[[i]]=list(mean=Yi[i-1],sd=sigma)
    }         
    Joint=rMvdc(1,mvdc(copula = ellipCopula(family = "normal", dim=m+1,
                                            dispstr = "un",param = c(R[lower.tri(R,diag=F)])),margins = c("weibull",rep("norm",m)),paramMargins=param))
    Ti=Joint[1]
    Yi=Joint[2:(m+1)]
    Tiobs=min(Ti,s[m],ceni)
    indi=as.numeric(Ti<s[m]&Ti<ceni)
    si=s[s<=Tiobs]
    dimyi=length(si)
    Datai=matrix(0,nrow=m,ncol=10);Dataobsi=matrix(0,nrow=dimyi,ncol=10)
    Datai[,1]=subj;Datai[,2]=Yi;Datai[,3]=s;Datai[,4]=Ti;Datai[,5]=Tiobs;Datai[,6]=ceni;Datai[,7]=indi;Datai[,8]=treat;Datai[,9]=gender;Datai[,10]=age
    Dataobsi[,1]=subj;Dataobsi[,2]=Yi[1:dimyi];Dataobsi[,3]=si;Dataobsi[,4]=Ti;Dataobsi[,5]=Tiobs;Dataobsi[,6]=ceni;Dataobsi[,7]=indi;Dataobsi[,8]=treat
    Dataobsi[,9]=gender;Dataobsi[,10]=age
    Dataobs=rbind(Dataobs,Dataobsi)
    Data=rbind(Data,Datai)
  }
  Data=data.frame(Data)
  Dataobs=data.frame(Dataobs)
  colnames(Data)=c('subj','resp','s','surti','obsti','ceni','indi','treat','gender','age')
  colnames(Dataobs)=c('subj','resp','s','surti','obsti','ceni','indi','treat','gender','age')
  results=list(Dataobs,Data)
  return(results)
}


Y.lmmsurtcopulaNOPLOT=function(beta1,beta2,s,n,r,sigma,rhoty,rhoy,rate,trprob,genprob,ageprob,lmmcor,surlmmcor,fred)
{
  m=length(s)
  Data=NULL
  Dataobs=NULL
  Rt=1
  if(lmmcor=="ex")  Ryt=rep(rhoty,m)
  if(lmmcor=="ar")  Ryt=rhoty^(m:1)
  if(surlmmcor=="ex") Ry=matrix(c(rep(c(1,rep(rhoy,m)),(m-1)),1),ncol=m)
  if(surlmmcor=="ar")  Ry=AR1(m,rhoy)
  Rty=t(Ryt)
  R=rbind(cbind(Rt,Rty),cbind(Ryt,Ry))
  Zi1=matrix(c(rep(1,m),s),ncol=2)
  for(subj in 1:n)
  {
    treat=rbinom(1,1,trprob)
    gender=rbinom(1,1,genprob)
    age=sample(c(0,1,2),1,prob=ageprob)
    Xi1=matrix(c(c(Zi1),rep(treat,m)*s,rep(gender,m),rep(as.numeric(age==1),m)
                 ,rep(as.numeric(age==2),m)),ncol=6)
    Xi2=matrix(c(1,treat,gender,as.numeric(age==1),as.numeric(age==2)),ncol=5)
    ceni=rexp(1,rate)
    Yi=Xi1%*%beta1
    scale=exp(-Xi2%*%beta2/r)
    param=list(list(shape=r,scale=scale))
    for(i in 2:(m+1))
    {
      param[[i]]=list(mean=Yi[i-1],sd=sigma)
    }         
    Joint=rMvdc(1,mvdc(copula = ellipCopula(family = "t", dim=m+1,df=fred, dispstr = "un",
                                            param = c(R[lower.tri(R,diag=F)])),margins = c("weibull",rep("norm",m)),paramMargins=param))
    Ti=Joint[1]
    Yi=Joint[2:(m+1)]
    Tiobs=min(Ti,s[m],ceni)
    indi=as.numeric(Ti<s[m]&Ti<ceni)
    si=s[s<=Tiobs]
    dimyi=length(si)
    Datai=matrix(0,nrow=m,ncol=10);Dataobsi=matrix(0,nrow=dimyi,ncol=10)
    Datai[,1]=subj;Datai[,2]=Yi;Datai[,3]=s;Datai[,4]=Ti;Datai[,5]=Tiobs;Datai[,6]=ceni;Datai[,7]=indi;Datai[,8]=treat;Datai[,9]=gender;Datai[,10]=age
    Dataobsi[,1]=subj;Dataobsi[,2]=Yi[1:dimyi];Dataobsi[,3]=si;Dataobsi[,4]=Ti;Dataobsi[,5]=Tiobs;Dataobsi[,6]=ceni;Dataobsi[,7]=indi;Dataobsi[,8]=treat
    Dataobsi[,9]=gender;Dataobsi[,10]=age
    Dataobs=rbind(Dataobs,Dataobsi)
    Data=rbind(Data,Datai)
  }
  Data=data.frame(Data)
  Dataobs=data.frame(Dataobs)
  colnames(Data)=c('subj','resp','s','surti','obsti','ceni','indi','treat','gender','age')
  colnames(Dataobs)=c('subj','resp','s','surti','obsti','ceni','indi','treat','gender','age')
  results=list(Dataobs,Data)
  return(results)
}


logliklmmsurGaucop=function(theta,data,lmmcor,surlmmcor,m)                #log-likelihood by assuming correlation among marginals introduced 
{                                                                         #by the multivariate Gaussian copula
  beta1=theta[1:6]
  beta2=theta[7:11]
  rhoty=theta[12]
  rhoy=theta[13]
  r=theta[14]
  sigma=theta[15]
  Rt=1
  if(lmmcor=="ex")  Ryt=rep(rhoty,m)
  if(lmmcor=="ar")  Ryt=rhoty^(m:1)
  if(surlmmcor=="ex") Ry=matrix(c(rep(c(1,rep(rhoy,m)),(m-1)),1),ncol=m)
  if(surlmmcor=="ar")  Ry=AR1(m,rhoy)
  Rty=t(Ryt)
  R=rbind(cbind(Rt,Rty),cbind(Ryt,Ry))
  n=length(unique(data$subj))
  ll=-1.0e40
  if(is.positive.definite(R)>0&sigma>0&r>0)
  {
    ll=0
    for(i in 1:n)
    {
      Xi1=matrix(c(rep(1,sum(data$subj==i)),data$s[data$subj==i],data$s[data$subj==i]*data$treat[data$subj==i],
                   data$gender[data$subj==i],as.numeric(data$age[data$subj==i]==1),as.numeric(data$age[data$subj==i]==2)),ncol=6)
      Xi2=c(1,unique(data$treat[data$subj==i]),unique(data$gender[data$subj==i]),unique(as.numeric(data$age[data$subj==i]==1)),
            unique(as.numeric(data$age[data$subj==i]==2)))
      obsti=unique(data$obsti[data$subj==i])
      indi=unique(data$indi[data$subj==i])
      si=data$s[data$subj==i]
      dimi=length(si)+1
      dimyi=dimi-1
      Ri=as.matrix(R[1:dimi,1:dimi])
      Ztobsi=qnorm(pweibull(obsti,shape=r,scale=exp(-t(Xi2)%*%beta2/r)))
      Zyobsi=c(data$resp[data$subj==i]-Xi1%*%beta1)/sigma
      if(indi==1)
      {
        Zobsi=c(Ztobsi,Zyobsi)
        ll=ll+(-dimyi*log(sigma)+log(dmvnorm(Zobsi,rep(0,dimi),Ri))+log(dweibull(obsti,shape=r,scale=exp(-t(Xi2)%*%beta2/r)))
               -log(dnorm(Ztobsi)))                                             #dmvnorm by default log=T
      }
      if(indi==0)
      {
        Ryi=as.matrix(Ri[-1,-1])
        Rtyi=matrix(Ri[1,-1],nrow=1)
        Ryti=matrix(Ri[-1,1],ncol=1)
        mutconi=Rtyi%*%solve(Ryi)%*%Zyobsi
        sigmatconi=sqrt(1-Rtyi%*%solve(Ryi)%*%Ryti)
        ll=ll+(-dimyi*log(sigma)+log(dmvnorm(Zyobsi,rep(0,dimyi),Ryi))+log(1-pnorm((Ztobsi-mutconi)/sigmatconi)))
      }
    }
  }
  if(is.na(ll)) ll=-1.0e40
  if(abs(ll)==Inf) ll=-1.0e40
  return(-ll)
}


logliklmmsurtcop=function(theta,data,lmmcor,surlmmcor,m,fred)                     #log-likelihood by assuming correlation among marginals introduced 
{                                                                                 #by the multivariate t copula
  beta1=theta[1:6]
  beta2=theta[7:11]
  rhoty=theta[12]
  rhoy=theta[13]
  r=theta[14]
  sigma=theta[15]
  Rt=1
  if(lmmcor=="ex")  Ryt=rep(rhoty,m)
  if(lmmcor=="ar")  Ryt=rhoty^(m:1)
  if(surlmmcor=="ex") Ry=matrix(c(rep(c(1,rep(rhoy,m)),(m-1)),1),ncol=m)
  if(surlmmcor=="ar")  Ry=AR1(m,rhoy)
  Rty=t(Ryt)
  R=rbind(cbind(Rt,Rty),cbind(Ryt,Ry))
  n=length(unique(data$subj))
  ll=-1.0e40
  if(is.positive.definite(R)>0&sigma>0&r>0)
  {
    ll=0
    for(i in 1:n)
    {
      Xi1=matrix(c(rep(1,sum(data$subj==i)),data$s[data$subj==i],data$s[data$subj==i]*data$treat[data$subj==i],
                   data$gender[data$subj==i],as.numeric(data$age[data$subj==i]==1),as.numeric(data$age[data$subj==i]==2)),ncol=6)
      Xi2=c(1,unique(data$treat[data$subj==i]),unique(data$gender[data$subj==i]),unique(as.numeric(data$age[data$subj==i]==1)),
            unique(as.numeric(data$age[data$subj==i]==2)))
      obsti=unique(data$obsti[data$subj==i])
      indi=unique(data$indi[data$subj==i])
      si=data$s[data$subj==i]
      dimi=length(si)+1
      dimyi=dimi-1
      Ri=R[1:dimi,1:dimi]
      Wtobsi=qt(pweibull(obsti,shape=r,scale=exp(-t(Xi2)%*%beta2/r)),df=fred)
      Wyobsi=c(qt(pnorm((data$resp[data$subj==i]-Xi1%*%beta1)/sigma),df=fred))
      if(indi==1)
      {
        Wobsi=c(Wtobsi,Wyobsi)
        ll=ll+(-dimyi*log(sigma)+log(dmvt(Wobsi,rep(0,dimi),Ri,fred,log=F))+sum(log(dnorm((data$resp[data$subj==i]-Xi1%*%beta1)/sigma)))
               -sum(log(dt(Wyobsi,fred)))+log(dweibull(obsti,shape=r,scale=exp(-t(Xi2)%*%beta2/r)))-log(dt(Wtobsi,df=fred)))   #dmvt by default log=T
      }                                                                                                                         #reset to log=F
      if(indi==0)
      {
        Ryi=as.matrix(Ri[-1,-1])
        Rtyi=matrix(Ri[1,-1],nrow=1)
        Ryti=matrix(Ri[-1,1],ncol=1)
        fredtconi=fred+dimyi
        mutconi=Rtyi%*%solve(Ryi)%*%Wyobsi
        sigmatconi=sqrt((fred+t(Wyobsi)%*%solve(Ryi)%*%Wyobsi)*(1-Rtyi%*%solve(Ryi)%*%Ryti)/(fred+dimyi))
        ll=ll+(-dimyi*log(sigma)+log(dmvt(Wyobsi,rep(0,dimyi),Ryi,fred,log=F))+sum(log(dnorm((data$resp[data$subj==i]-Xi1%*%beta1)/sigma)))
               -sum(log(dt(Wyobsi,fred)))+log(1-pt((Wtobsi-mutconi)/sigmatconi,fredtconi)))                      #dmvt by default log=T
      }                                                                                                           #reset to log=F
    }
  }
  if(is.na(ll)) ll=-1.0e40
  if(abs(ll)==Inf) ll=-1.0e40
  return(-ll)
}


Y.lmmsurtcopest.nlm=function(theta,m,Data,lmmcor,surlmmcor,fred,N)       #function for estimation by assuming correlation among the marginals  
{                                                                       #is introduced by the multivariate t copula
  estmatrix=matrix(0,nrow=N,ncol=15)
  gradientmatrix=matrix(0,nrow=N,ncol=15)
  sdmatrix=matrix(0,nrow=N,ncol=15)
  samplestmatrix=matrix(0,nrow=1,ncol=15)
  rebiasmatrix=matrix(0,nrow=N,ncol=15)
  RMSEmatrix=matrix(0,nrow=N,ncol=15)
  CPmatrix=matrix(0,nrow=N,ncol=15)
  ECPmatrix=matrix(0,nrow=N,ncol=15)
  loglik=0
  colnames(estmatrix)=c('beta01','beta11','beta21','beta31','beta41','beta51','beta02','beta12','beta22','beta32',
                        'beta42','rho1','rho2','r','sigma')
  colnames(gradientmatrix)=c('gbeta01','gbeta11','gbeta21','gbeta31','gbeta41','gbeta51','gbeta02','gbeta12','gbeta22','gbeta32',
                             'gbeta42','grho1','grho2','gr','gsigma')
  colnames(sdmatrix)=c('sdbeta01','sdbeta11','sdbeta21','sdbeta31','sdbeta41','sdbeta51','sdbeta02','sdbeta12',
                       'sdbeta22','sdbeta32','sdbeta42','sdrho1','sdrho2','sdr','sdsigma')
  colnames(samplestmatrix)=c('ssdbeta01','ssdbeta11','ssdbeta21','ssdbeta31','ssdbeta41','ssdbeta51','ssdbeta02','ssdbeta12',
                             'ssdbeta22','ssdbeta32','ssdbeta42','ssdrho1','ssdrho2','ssdr','ssdsigma')
  colnames(rebiasmatrix)=c('rbbeta01','rbbeta11','rbbeta21','rbbeta31','rbbeta41','rbbeta51','rbbeta02','rbbeta12',
                           'rbbeta22','rbbeta32','rbbeta42','rbrho1','rbrho2','rbr','rbsigma')
  colnames(RMSEmatrix)=c('rmbeta01','rmbeta11','rmbeta21','rmbeta31','rmbeta41','rmbeta51','rmbeta02','rmbeta12',
                         'rmbeta22','rmbeta32','rmbeta42','rmrho1','rmrho2','rmr','rmsigma')
  colnames(CPmatrix)=c('cpbeta01','cpbeta11','cpbeta21','cpbeta31','cpbeta41','cpbeta51','cpbeta02','cpbeta12',
                       'cpbeta22','cpbeta32','cpbeta42','cprho1','cprho2','cpr','cpsigma')
  colnames(ECPmatrix)=c('ecpbeta01','ecpbeta11','ecpbeta21','ecpbeta31','ecpbeta41','ecpbeta51','ecpbeta02','ecpbeta12',
                        'ecpbeta22','ecpbeta32','ecpbeta42','ecprho1','ecprho2','ecpr','ecpsigma')
  i=1
  j=0
  datasetid=NULL
  while(i<=N)
  {
    j=j+1
    esttcop=nlm(f=logliklmmsurtcop,p=theta,m=m,fred=fred,data=Data[[j]],lmmcor=lmmcor,surlmmcor=surlmmcor,hessian=T,iterlim=1000)
    if(esttcop$code==1&is.positive.definite(esttcop$hessian)>0)
    {
      estmatrix[i,]=esttcop$estimate
      gradientmatrix[i,]=esttcop$gradient
      sdmatrix[i,]=sqrt(diag(solve(esttcop$hessian)))
      rebiasmatrix[i,]=esttcop$estimate/theta-1
      RMSEmatrix[i,]=(esttcop$estimate-theta)^2
      loglik[i]=-esttcop$minimum
      print(c(j,i))
      datasetid=c(datasetid,j)
      i=i+1
    }
  }
  samplestmatrix[1,]=apply(estmatrix, 2, sd)
  eupmatrix=estmatrix+matrix(rep(qnorm(0.975)*samplestmatrix,N),nrow=N,byrow=T)
  elowmatrix=estmatrix-matrix(rep(qnorm(0.975)*samplestmatrix,N),nrow=N,byrow=T)
  upmatrix=estmatrix+qnorm(0.975)*sdmatrix
  lowmatrix=estmatrix-qnorm(0.975)*sdmatrix
  for(k in 1:N)
  {
    ECPmatrix[k,]=as.numeric(elowmatrix[k,]<theta&theta<eupmatrix[k,])
    CPmatrix[k,]=as.numeric(lowmatrix[k,]<theta&theta<upmatrix[k,])
  }
  results=list(estmatrix,gradientmatrix,sdmatrix,colMeans(estmatrix),colMeans(sdmatrix),colMeans(rebiasmatrix),colMeans(RMSEmatrix)^(0.5),
               colMeans(ECPmatrix),colMeans(CPmatrix),loglik,j,datasetid)
  return(results)
}


Y.lmmsurGaucopest.nlm=function(theta,m,Data,lmmcor,surlmmcor,N)            #function for estimation by assuming correlation among the marginals  
{                                                                          #is introduced by the multivariate Gaussian copula
  estmatrix=matrix(0,nrow=N,ncol=15)
  gradientmatrix=matrix(0,nrow=N,ncol=15)
  sdmatrix=matrix(0,nrow=N,ncol=15)
  samplestmatrix=matrix(0,nrow=1,ncol=15)
  rebiasmatrix=matrix(0,nrow=N,ncol=15)
  RMSEmatrix=matrix(0,nrow=N,ncol=15)
  CPmatrix=matrix(0,nrow=N,ncol=15)
  ECPmatrix=matrix(0,nrow=N,ncol=15)
  loglik=0
  colnames(estmatrix)=c('beta01','beta11','beta21','beta31','beta41','beta51','beta02','beta12','beta22','beta32',
                        'beta42','rho1','rho2','r','sigma')
  colnames(gradientmatrix)=c('gbeta01','gbeta11','gbeta21','gbeta31','gbeta41','gbeta51','gbeta02','gbeta12','gbeta22','gbeta32',
                             'gbeta42','grho1','grho2','gr','gsigma')
  colnames(sdmatrix)=c('sdbeta01','sdbeta11','sdbeta21','sdbeta31','sdbeta41','sdbeta51','sdbeta02','sdbeta12',
                       'sdbeta22','sdbeta32','sdbeta42','sdrho1','sdrho2','sdr','sdsigma')
  colnames(samplestmatrix)=c('ssdbeta01','ssdbeta11','ssdbeta21','ssdbeta31','ssdbeta41','ssdbeta51','ssdbeta02','ssdbeta12',
                             'ssdbeta22','ssdbeta32','ssdbeta42','ssdrho1','ssdrho2','ssdr','ssdsigma')
  colnames(rebiasmatrix)=c('rbbeta01','rbbeta11','rbbeta21','rbbeta31','rbbeta41','rbbeta51','rbbeta02','rbbeta12',
                           'rbbeta22','rbbeta32','rbbeta42','rbrho1','rbrho2','rbr','rbsigma')
  colnames(RMSEmatrix)=c('rmbeta01','rmbeta11','rmbeta21','rmbeta31','rmbeta41','rmbeta51','rmbeta02','rmbeta12',
                         'rmbeta22','rmbeta32','rmbeta42','rmrho1','rmrho2','rmr','rmsigma')
  colnames(CPmatrix)=c('cpbeta01','cpbeta11','cpbeta21','cpbeta31','cpbeta41','cpbeta51','cpbeta02','cpbeta12',
                       'cpbeta22','cpbeta32','cpbeta42','cprho1','cprho2','cpr','cpsigma')
  colnames(ECPmatrix)=c('ecpbeta01','ecpbeta11','ecpbeta21','ecpbeta31','ecpbeta41','ecpbeta51','ecpbeta02','ecpbeta12',
                        'ecpbeta22','ecpbeta32','ecpbeta42','ecprho1','ecprho2','ecpr','ecpsigma')
  i=1
  j=0
  datasetid=NULL
  while(i<=N)
  {
    j=j+1
    estGaucop=nlm(f=logliklmmsurGaucop,p=theta,m=m,data=Data[[j]],lmmcor=lmmcor,surlmmcor=surlmmcor,hessian=T,iterlim=1000)
    if(estGaucop$code==1&is.positive.definite(estGaucop$hessian)>0)
    {
      estmatrix[i,]=estGaucop$estimate
      gradientmatrix[i,]=estGaucop$gradient
      sdmatrix[i,]=sqrt(diag(solve(estGaucop$hessian)))
      rebiasmatrix[i,]=estGaucop$estimate/theta-1
      RMSEmatrix[i,]=(estGaucop$estimate-theta)^2
      loglik[i]=-estGaucop$minimum
      print(c(j,i))
      datasetid=c(datasetid,j)
      i=i+1
    }
  }
  samplestmatrix[1,]=apply(estmatrix, 2, sd)
  eupmatrix=estmatrix+matrix(rep(qnorm(0.975)*samplestmatrix,N),nrow=N,byrow=T)
  elowmatrix=estmatrix-matrix(rep(qnorm(0.975)*samplestmatrix,N),nrow=N,byrow=T)
  upmatrix=estmatrix+qnorm(0.975)*sdmatrix
  lowmatrix=estmatrix-qnorm(0.975)*sdmatrix
  for(k in 1:N)
  {
    ECPmatrix[k,]=as.numeric(elowmatrix[k,]<theta&theta<eupmatrix[k,])
    CPmatrix[k,]=as.numeric(lowmatrix[k,]<theta&theta<upmatrix[k,])
  }
  results=list(estmatrix,gradientmatrix,sdmatrix,colMeans(estmatrix),colMeans(sdmatrix),colMeans(rebiasmatrix),colMeans(RMSEmatrix)^(0.5),
               colMeans(ECPmatrix),colMeans(CPmatrix),loglik,j,datasetid)
  return(results)
}


#Dynamic prediction for survival process on simulated data
simpredSurGaucop=function(beta1,beta2,rho1,rho2,r,sigma,data,dynati,predinterv,acc)      #Prediction based on multivariate Gaussian copula joint model
{
  timepoint=data$s
  ni=length(timepoint)
  R21=rep(rho1,ni)
  R12=t(R21)
  R11=1
  R22=matrix(c(rep(c(1,rep(rho2,ni)),(ni-1)),1),ncol=ni)
  R=rbind(cbind(R11,R12),cbind(R21,R22))
  Xi1=matrix(c(rep(1,ni),timepoint,timepoint*data$treat,
               data$gender,as.numeric(data$age==1),as.numeric(data$age==2)),ncol=6)
  Xi2=c(1,unique(data$treat),unique(data$gender),unique(as.numeric(data$age==1)),unique(as.numeric(data$age==2)))
  if(dynati>=unique(data$obsti))  
  {
    t=seq(unique(data$obsti),predinterv,by=acc)
    dimyi=sum(timepoint<=unique(data$obsti))
  }
  if(dynati<unique(data$obsti))  
  {
    t=seq(dynati,predinterv,by=acc)
    dimyi=sum(timepoint<=dynati)
  }
  Ryi=R[-1,-1][1:dimyi,1:dimyi]
  Si=1-pweibull(t,shape=r,scale=exp(-t(Xi2)%*%beta2/r))
  Zyi=c(data$resp-Xi1%*%beta1)[1:dimyi]/sigma
  mutconi=as.numeric(t(R[1,-1][1:dimyi])%*%solve(Ryi)%*%Zyi)
  sigmatconi=as.numeric(sqrt(1-t(R[1,-1][1:dimyi])%*%solve(Ryi)%*%R[1,-1][1:dimyi]))
  predi=pnorm((qnorm(Si)+mutconi)/sigmatconi)/pnorm((qnorm(Si[1])+mutconi)/sigmatconi)
  results=list(t,Si/Si[1],predi,dimyi)
  return(results)
}



simpredSurtcop=function(beta1,beta2,rho1,rho2,r,sigma,data,dynati,predinterv,fred,acc)       #Prediction based on multivariate t copula joint model
{
  timepoint=data$s
  ni=length(timepoint)
  R21=rep(rho1,ni)
  R12=t(R21)
  R11=1
  R22=matrix(c(rep(c(1,rep(rho2,ni)),(ni-1)),1),ncol=ni)
  R=rbind(cbind(R11,R12),cbind(R21,R22))
  Xi1=matrix(c(rep(1,ni),timepoint,timepoint*data$treat,
               data$gender,as.numeric(data$age==1),as.numeric(data$age==2)),ncol=6)
  Xi2=c(1,unique(data$treat),unique(data$gender),unique(as.numeric(data$age==1)),unique(as.numeric(data$age==2)))
  if(dynati>=unique(data$obsti))  
  {
    t=seq(unique(data$obsti),predinterv,by=acc)
    dimyi=sum(timepoint<=unique(data$obsti))
  }
  if(dynati<unique(data$obsti))  
  {
    t=seq(dynati,predinterv,by=acc)
    dimyi=sum(timepoint<=dynati)
  }
  Ryi=R[-1,-1][1:dimyi,1:dimyi]
  Rtyi=matrix(R[1,-1][1:dimyi],nrow=1)
  Ryti=matrix(R[-1,1][1:dimyi],ncol=1)
  Si=1-pweibull(t,shape=r,scale=exp(-t(Xi2)%*%beta2/r))
  Wyi=c(qt(pnorm((data$resp-Xi1%*%beta1)[1:dimyi]/sigma),df=fred))
  fredtconi=fred+dimyi
  mutconi=as.numeric(Rtyi%*%solve(Ryi)%*%Wyi)
  sigmatconi=as.numeric(sqrt((fred+t(Wyi)%*%solve(Ryi)%*%Wyi)*(1-Rtyi%*%solve(Ryi)%*%Ryti)/(fred+dimyi)))
  predi=pt((qt(Si,df=fred)+mutconi)/sigmatconi,df=fredtconi)/pt((qt(Si[1],df=fred)+mutconi)/sigmatconi,df=fredtconi)
  results=list(t,predi,dimyi)
  return(results)
}


#Claculate the absolute different between the true resudial life time and the estimated residual life time
meanresditi=function(Data,beta1,beta2,rho1,rho2,r,sigma,dynati,predinterv,acc,n,cop,fred)
{
  estleft=0
  truleft=0
  j=0
  i=1
  while(i<=n)
  { 
    Datai=Data[Data$subj==i,]
    s=Datai$s
    if(max(s)>=dynati)
    { 
      j=j+1
    if(cop=="nocop")   estleft[j]=sum(simpredSurGaucop(beta1=beta1,beta2=beta2,rho1=rho1,rho2=rho2,r=r,sigma=sigma,      #for survival sub-model alone
                                                         data=Datai,dynati=dynati,predinterv=predinterv,acc=acc)[[2]]*acc)
    if(cop=="Gaucop")   estleft[j]=sum(simpredSurGaucop(beta1=beta1,beta2=beta2,rho1=rho1,rho2=rho2,r=r,sigma=sigma,      #for Gaussian joint model
                                              data=Datai,dynati=dynati,predinterv=predinterv,acc=acc)[[3]]*acc)
    if(cop=="tcop")   estleft[j]=sum(simpredSurtcop(beta1=beta1,beta2=beta2,rho1=rho1,rho2=rho2,r=r,sigma=sigma,
                                              data=Datai,dynati=dynati,predinterv=predinterv,fred=fred,acc=acc)[[2]]*acc) ##for t joint model
    truleft[j]=unique(Datai$surti)-dynati
    }
    i=i+1
  }
  return(abs(estleft-truleft))
}
