################################################################
### Simulation code under violation of M4' (c): 
### additional M-W interaction in outcome model E[Y|A,M,C,U,W]
### Xu Shi (xushi@hsph.harvard.edu)
### Isabel Fulcher (isabelfulcher@g.harvard.edu)
################################################################
rm(list=ls())
library(foreign)
library(sandwich)
library(lmtest)
library(boot)
n.rep=100
my.filepath=""
source(paste0(my.filepath,"GENIUS_func.R"))
expit = function(x) {exp(x)/(1+exp(x))}

## Coefficients
cbeta = c(c1=.5,c2=.5)
wbeta = c(intercept=0,c1=.2,c2=.2)
ubeta = c(intercept=0,c1=.2,c2=.2)
abeta = c(intercept=0,W=1,c1=.2,c2=.2)  # A=(1,W,C1,C2)%*%abeta
mbeta = c(intercept=0,A=1,U=1,c1=.2,c2=.2); # M=(1,A,U,C1,C2)%*%meta
ybeta = c(intercept=0,M=1,A=1,W=-1,U=-1,c1=.2,c2=.2,MW=0.5) # Y=(1,M,A,W,U,C1,C2,MW)%*%ybeta
myparam = list(abeta=abeta,mbeta=mbeta,ybeta=ybeta,cbeta=cbeta,wbeta=wbeta,ubeta=ubeta)
## true value
beta.MonY = myparam$ybeta["M"]
beta.AonM = myparam$mbeta["A"]
mysamplesize=1000
myR=2*1000

gen.data = function(N,param,DAG,addcov){
  ### add two covariates: one continues one binary
  if(addcov==T){
    C1 = rnorm(N,mean=param$cbeta["c1"],sd=2)
    C2 = rbinom(N,1,plogis(param$cbeta["c2"]))
  }else{
    C1 = C2 = rep(0,N)
  }
  ### generate data under true DGM
  if(DAG=="a"){
    U = E = 0
    W = cbind(1,C1,C2)%*%param$wbeta + rnorm(N,mean=0,sd=1)
  }else if(DAG=="b"){
    W = E = 0
    U = cbind(1,C1,C2)%*%param$ubeta + rnorm(N,mean=0,sd=1)
  }else if(DAG=="c"){
    E = 0
    W = cbind(1,C1,C2)%*%param$wbeta + rnorm(N,mean=0,sd=1)
    U = cbind(1,C1,C2)%*%param$ubeta + rnorm(N,mean=0,sd=1)
  }else{
    E = rnorm(N,mean=0,sd=1)
    W = cbind(1,C1,C2)%*%param$wbeta + rnorm(N,mean=0,sd=1)
    U = cbind(1,C1,C2)%*%param$ubeta + rnorm(N,mean=0,sd=1)
  }
  A = rbinom(n=N,size=1,p=expit(cbind(1,W,C1,C2)%*%param$abeta))
  ### intentionally make var(M|A)=f(A)
  lambda0=1; lambda1=0.5;
  M = cbind(1,A,U,C1,C2)%*%param$mbeta + 
    rnorm(N,mean=0,sd=abs(lambda0+as.vector(lambda1%*%t(A))))
  Y = cbind(1,M,A,W,U,C1,C2,M*W)%*%param$ybeta + rnorm(N,mean=0,sd=1)
  ### add in an indep meanzero measurement error
  Mstar = M+E
  return(data.frame(C1=C1,C2=C2,M=M,Mstar=Mstar,A=A,Y=Y,W=W,U=U))
}

run.one = function(N,param,DAG,addcov){ 
  data = gen.data(N,param,DAG,addcov)
  A = data$A; M = data$M; Mstar = data$Mstar; Y = data$Y
  C1 = data$C1; C2 = data$C2; C=cbind(data$C1,data$C2); W = data$W; U = data$U
  ##### FIRST get effect of A on M w/ robust se
  ## M~A
  m_Ma.OLS = lm(Mstar~A+C1+C2)
  (beta.Ma = m_Ma.OLS$coef["A"])
  (var.beta.Ma.robust = coeftest(m_Ma.OLS, vcov = sandwich)["A",2]^2)
  (var.beta.Ma.nonrobust = (summary(m_Ma.OLS)$coef["A",2])^2)
  
  ##### SECOND get effect of M on Y
  ###(1) OLS to obtain beta_m 
  m_OLS = lm(Y~Mstar+A+C1+C2)
  (ols.beta.Ym = m_OLS$coef["Mstar"])
  (ols.NIE = beta.Ma*ols.beta.Ym)
  ## CI of naive OLS
  var.OLS.beta.Ym = (summary(m_OLS)$coef["Mstar",2])^2
  var.OLS.NIE.robust = var.beta.Ma.robust*ols.beta.Ym^2+var.OLS.beta.Ym*beta.Ma^2
  ci.ols.robust = ols.NIE + c(-1,1)*qnorm(1-0.05/2)*sqrt(var.OLS.NIE.robust)
  var.OLS.NIE.nonrobust = var.beta.Ma.nonrobust*ols.beta.Ym^2+var.OLS.beta.Ym*beta.Ma^2
  ci.ols.nonrobust = ols.NIE + c(-1,1)*qnorm(1-0.05/2)*sqrt(var.OLS.NIE.nonrobust)
  ###(2) GENIUS to obtain beta_m 
  if(addcov==T){
    m_genius = genius_addY_withC(Y,A=Mstar,G=A,C=C,formula=A~G+C,alpha=0.05,lower=-10,upper=10)
  }else{
    m_genius = genius_addY(Y,A=Mstar,G=A,formula=A~G,alpha=0.05,lower=-10,upper=10)
  }
  (genius.beta.Ym = m_genius$beta.est)
  (genius.NIE = beta.Ma*genius.beta.Ym)
  ## CI of GENIUS NIE
  var.genius.beta.Ym = m_genius$beta.var
  var.genius.NIE.robust = var.beta.Ma.robust*genius.beta.Ym^2+var.genius.beta.Ym*beta.Ma^2
  ci.genius.robust = genius.NIE + c(-1,1)*qnorm(1-0.05/2)*sqrt(var.genius.NIE.robust)
  var.genius.NIE.nonrobust = var.beta.Ma.nonrobust*genius.beta.Ym^2+var.genius.beta.Ym*beta.Ma^2
  ci.genius.nonrobust = genius.NIE + c(-1,1)*qnorm(1-0.05/2)*sqrt(var.genius.NIE.nonrobust)
  ###(3) ORACLE -- observe U W M
  if(DAG=="a"){
    m_OLS_oracle = lm(Y~M*W+A+C1+C2)
    m_Ma.OLS.oracle = lm(M~A+C1+C2)
  } else if(DAG=="b"){
    m_OLS_oracle = lm(Y~M+A+U+C1+C2)
    m_Ma.OLS.oracle = lm(M~A+U+C1+C2)
  } else if(DAG=="c"){
    m_OLS_oracle = lm(Y~M*W+U+A+C1+C2)
    m_Ma.OLS.oracle = lm(M~A+U+C1+C2)  # m_Ma.OLS.oracle = lm(M~A+U+W+C1+C2)
  } else if(DAG=="d"){
    m_OLS_oracle = lm(Y~M*W+U+A+C1+C2)
    m_Ma.OLS.oracle = lm(M~A+U+C1+C2)  # m_Ma.OLS.oracle = lm(M~A+U+W+C1+C2)
  }
  (beta.Ma.oracle = m_Ma.OLS.oracle$coef["A"])
  if(DAG=="b"){
    (beta.Ym.oracle = m_OLS_oracle$coef["M"])
  }else{
    (beta.Ym.oracle = m_OLS_oracle$coef["M"]+m_OLS_oracle$coef["M:W"]*mean(W))##this is diff than no violation
  }
  (param$ybeta["M"]+param$ybeta["MW"]*mean(W)) ## true beta.Ym
  (NIE.oracle = beta.Ym.oracle*beta.Ma.oracle)
  (NIE.true = (param$ybeta["M"]+param$ybeta["MW"]*mean(W))*param$mbeta["A"])## true NIE
  ## CI of ORACLE
  (var.beta.Ma.oracle.robust = coeftest(m_Ma.OLS.oracle, vcov = sandwich)["A",2]^2)
  (var.beta.Ma.oracle.nonrobust = (summary(m_Ma.OLS.oracle)$coef["A",2])^2)
  ## diff from the no violation setting, here we need to get var of (\theta_M+\theta_{MW}*E[W])
  p=length(which(!is.na(coef(m_OLS_oracle))))
  Ym.vec = c(0,1,rep(0,p-3),mean(W))
  var.beta.Ym.oracle=Ym.vec%*%sandwich(m_OLS_oracle)%*%Ym.vec
  #beta.Ym.oracle and var.beta.Ym.oracle are diff from no vio 
  var.oracle.NIE.robust = as.numeric(var.beta.Ma.oracle.robust*beta.Ym.oracle^2+
                                       var.beta.Ym.oracle*beta.Ma.oracle^2)
  var.oracle.NIE.nonrobust = as.numeric(var.beta.Ma.oracle.nonrobust*beta.Ym.oracle^2+
                                          var.beta.Ym.oracle*beta.Ma.oracle^2)
  ci.oracle.robust = NIE.oracle + c(-1,1)*qnorm(1-0.05/2)*sqrt(var.oracle.NIE.robust)#no need /sqrt(n)
  ci.oracle.nonrobust = NIE.oracle + c(-1,1)*qnorm(1-0.05/2)*sqrt(var.oracle.NIE.nonrobust)
  
  #####BOOT: now compute bootstrapped CI for OLS, genius, oracle
  boot.type="perc"
  ###(1) bootstrap ci for OLS
  boot.func = function(data, indices,DAG){
    d = data[indices,] # allows boot to select sample 
    fit.theta_m = lm(Y~Mstar+A+C1+C2, data=d)
    fit.beta_a = lm(Mstar~A+C1+C2, data=d)
    NIE = fit.theta_m$coef["Mstar"]*fit.beta_a$coef["A"]
    return(NIE)
  }
  results = boot(data=data.frame(Y=Y,Mstar=Mstar,A=A,C1=C1,C2=C2), statistic=boot.func,DAG=DAG,R=myR)
  ci.ols.boot = boot.ci(results, type=boot.type)$percent[1,4:5]
  ###(2) bootstrap ci for GENIUS
  boot.func = function(data,indices,DAG){
    d = data[indices,] # allows boot to select sample 
    if(addcov==T){
      fit.theta_m = genius_addY_withC(d$Y,A=d$Mstar,G=d$A,C=cbind(d$C1,d$C2),formula=A~G+C,alpha=0.05,lower=-10,upper=10)
    }else{
      fit.theta_m = genius_addY(d$Y,A=d$Mstar,G=d$A,formula=A~G,alpha=0.05,lower=-10,upper=10) 
    }
    fit.beta_a = lm(Mstar~A+C1+C2, data=d)
    NIE = fit.theta_m$beta.est*fit.beta_a$coef["A"]
    return(NIE)
  }
  results = boot(data=data.frame(Y=Y,Mstar=Mstar,A=A,C1=C1,C2=C2), statistic=boot.func,DAG=DAG,R=myR)
  ci.genius.boot = boot.ci(results, type=boot.type)$percent[1,4:5]
  ###(3) bootstrap ci for oracle
  boot.func = function(data,indices,DAG){
    d = data[indices,] # allows boot to select sample 
    if(DAG=="a"){
      fit.theta_m = lm(Y~M*W+A+C1+C2, data=d)
      fit.beta_a = lm(M~A+C1+C2, data=d)
    } else if(DAG=="b"){
      fit.theta_m = lm(Y~M+A+U+C1+C2, data=d)
      fit.beta_a = lm(M~A+U+C1+C2, data=d)
    } else if(DAG=="c"){
      fit.theta_m = lm(Y~M*W+A+U+C1+C2, data=d)
      fit.beta_a = lm(M~A+U+C1+C2, data=d)     # fit.beta_a = lm(M~A+U+W+C1+C2, data=d)
    } else if(DAG=="d"){
      fit.theta_m = lm(Y~M*W+A+U+C1+C2, data=d)
      fit.beta_a = lm(M~A+U+C1+C2, data=d)     # fit.beta_a = lm(M~A+U+W+C1+C2, data=d)
    }
    if(DAG=="b"){
      NIE = fit.theta_m$coef["M"]*fit.beta_a$coef["A"]
    }else{
      (NIE = (fit.theta_m$coef["M"]+fit.theta_m$coef["M:W"]*mean(W))*
         fit.beta_a$coef["A"])
    }
    return(NIE)  
  }
  results = boot(data=data.frame(Y=Y,M=M,A=A,W=W,U=U,C1=C1,C2=C2), statistic=boot.func, DAG=DAG,R=myR)
  ci.oracle.boot = boot.ci(results, type=boot.type)$percent[1,4:5]
  
  
  ### collect rslts
  rslt = c(ols.beta.Ym=ols.beta.Ym,ols.NIE=ols.NIE,  
           genius.beta.Ym=genius.beta.Ym,genius.NIE=genius.NIE,
           beta.Ma=beta.Ma,NIE.true=NIE.true,NIE.oracle=NIE.oracle,
           ci.ols.robust=ci.ols.robust,ci.genius.robust=ci.genius.robust,ci.oracle.robust=ci.oracle.robust,
           ci.ols.nonrobust=ci.ols.nonrobust,ci.genius.nonrobust=ci.genius.nonrobust,ci.oracle.nonrobust=ci.oracle.nonrobust,
           ci.ols.boot=ci.ols.boot,ci.genius.boot=ci.genius.boot,ci.oracle.boot=ci.oracle.boot,
           var.OLS.NIE.robust=var.OLS.NIE.robust,var.genius.NIE.robust=var.genius.NIE.robust,var.oracle.NIE.robust=var.oracle.NIE.robust,
           var.OLS.NIE.nonrobust=var.OLS.NIE.nonrobust,var.genius.NIE.nonrobust=var.genius.NIE.nonrobust,var.oracle.NIE.nonrobust=var.oracle.NIE.nonrobust)
  return(rslt)
}



simu_rslt.a = run.one(N=mysamplesize,param=myparam,DAG="a",addcov=T)
simu_rslt.b = run.one(N=mysamplesize,param=myparam,DAG="b",addcov=T)
simu_rslt.c = run.one(N=mysamplesize,param=myparam,DAG="c",addcov=T)
simu_rslt.d = run.one(N=mysamplesize,param=myparam,DAG="d",addcov=T)
