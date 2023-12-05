
# EM algorithm using coxph(), glm(), glm()

expit.f <- function(x){u=exp(x); u/(1+u)}

curehaz <- function(Time,Status,wcoxexp){    
  death_point <- sort((subset(Time, Status==1)))#sort(unique(subset(Time, Status==1)))
  #if(model=='ph') coxexp <- exp((beta)%*%t(X))  
  lambda <- numeric()
  event <- numeric()
  for(i in 1: length(death_point)){
    event[i] <- sum(Status*as.numeric(Time==death_point[i]))
    #if(model=='ph')  temp <- sum(as.numeric(Time>=death_point[i])*w*drop(coxexp))
    #if(model=='aft')  temp <- sum(as.numeric(Time>=death_point[i])*w)
    temp <- sum(as.numeric(Time>=death_point[i])*wcoxexp)
    temp1 <- event[i]
    lambda[i] <- temp1/temp
  }
  
  list(hazmass=lambda)
}

p_y_g_x_coxph.f <- function(del, t, X, Z, beta, eta, hf0, survf0){
  covZ = drop(eta%*%t(Z))#eta[1] + eta[2]*x1 + eta[3]*x2;
  piX = expit.f(covZ)
  covT = drop(beta%*%t(X))#beta[1]*x1 +beta[2]*x2 #beta[1]*x1 + beta[2]*x2 
  hf<-hf0*exp(covT)
  survf<-survf0^(exp(covT))
  densf <- hf*survf*piX
  res0 <- 1-piX*(1-survf)
  res <- ifelse(del==1, densf, res0) #del*densf+(1-del)*res0  
  return(res)
}

p_x1_g_x2.f <- function(x1, V, gam){
  
  cov = drop(gam%*%t(V))
  p = expit.f(cov);
  res <- x1*p+(1-x1)*(1-p) 
  
  return(res)
}

wgt1_zeta_coxph.f <- function(scen, r, x1, del, t, X, Z, V, beta, eta, gam, hf0, survf0){
  
  #if(model=="binary"){
  
  X0=X1=X; Z0=Z1=Z
  if(scen=="sc"){
    X1[,ncol(X)]<-1; X0[,ncol(X)]<-0
    Z1[,ncol(Z)]<-1; Z0[,ncol(Z)]<-0
  }else if(scen=="sb"){
    X1[,ncol(X)]<-1; X0[,ncol(X)]<-0
    #Z1[,ncol(Z)]<-1; Z0[,ncol(Z)]<-0
  }else{
    #X1[,ncol(X)]<-1; X0[,ncol(X)]<-0
    Z1[,ncol(Z)]<-1; Z0[,ncol(Z)]<-0
  }
  
  num01 <- p_y_g_x_coxph.f(del, t, X0, Z0, beta, eta, hf0, survf0)
  num02 <- p_x1_g_x2.f(0, V, gam)
  
  num11 <- p_y_g_x_coxph.f(del, t, X1, Z1, beta, eta, hf0, survf0)
  num12 <- p_x1_g_x2.f(1, V, gam)
  num0 = num01*num02; num1=num11*num12
  denom <- num0 + num1
  res0 <- ((1-x1)*num0+x1*num1)/denom
  
  res <- ifelse(r==1, r, res0)
  
  return(res)
}

# estimation in phase I
em_coxph_H0 <-function(dat.orig, covXname, covZname, eta,beta,emmax,eps){     
  
  Status=dat.orig$del
  Time=dat.orig$time
  event_point=sort((subset(Time, Status==1)))#sort(unique(subset(Time, Status==1)))
  
  orderid=order(subset(Time, Status==1))
  ID=unique(dat.orig$id)
  useID=c((ID[Status==1])[orderid], ID[Status==0])
  
  dat.x <- as.data.frame(dat.orig)
  
  # expand the dat.x for susceptibility model 
  pz <- c(rep(1, dim(dat.x)[1]), rep(0, sum(1-dat.x$del)))
  dat.z <- cbind(rbind(subset(dat.x, del==1),
                       subset(dat.x, del==0),
                       subset(dat.x, del==0)),pz)
  dat.z <- dat.z[order(dat.z$id),]
  
  # weights for latent susceptibility
  w2<-dat.z$del
  dat.z$wgt2<-ifelse(dat.z$pz==1,w2,1-w2)
  
  dat.z$wgt<-dat.z$wgt2
  
  if(!is.null(beta)){
    X=as.matrix(dat.z[,covXname,drop=FALSE])
    covT=drop(exp(beta%*%t(X)))
  }else{
    covT=1.0
  }
  Z=as.matrix(cbind(rep(1,nrow(dat.z)), dat.z[,covZname,drop=FALSE]))
  
  wcoxexp0 <- dat.z$wgt*covT
  test0=data.frame(id=subset(dat.z, pz==1)$id, val=wcoxexp0[dat.z$pz==1])
  wcoxexp=aggregate(val~id, data=test0, sum)$val #test0$val#
  hazmass=curehaz(Time,Status,wcoxexp)$hazmass
  Hazf<-stepfun(event_point,c(0,cumsum(hazmass)))
  Hazv<-Hazf(dat.z$time)
  Hazv<-ifelse(dat.z$time>max(event_point), Inf, Hazv)
  Hazv<-ifelse(dat.z$time<min(event_point), 0, Hazv)
  s0=exp(-Hazv)
  
  df1h<-data.frame(id=useID, time=c(event_point,Time[Status==0]), lam0=c(hazmass, rep(1, sum(1-Status))))
  df2h<-dat.z
  dfh<- right_join(df1h,df2h)
  dfh=dfh[order(dfh$id),]
  
  
  
  
  convergence<- 1000;i <-1
  while (convergence > eps & i < emmax){  
    
    
    uncureprob <- matrix(exp((eta)%*%t(Z))/(1+exp((eta)%*%t(Z))),ncol=1)
    if(is.null(beta)){
      coxexp=1.0
    }else{
      coxexp=drop(exp((beta)%*%t(X)))
    }
    survival<-s0^(coxexp)
    #haz<-dfh$lam0*coxexp
    
    ## E step 
    # w1 <- wgt1_zeta_coxph.f(dat.z$r, dat.z$px, dat.z$del, dat.z$time, X, Z, V, beta, eta, gam, dfh$lam0, s0)
    w2 <- dat.z$del+(1-dat.z$del)*(uncureprob*survival)/((1-uncureprob)+uncureprob*survival)
    
    #dat.z$wgt1 <- w1
    dat.z$wgt2 <- ifelse(dat.z$pz==1, w2, 1-w2)
    #dat.z$wgt <- dat.z$wgt1*dat.z$wgt2
    
    ## M step
    
    #update_gam <- glm(px~1+v, family=quasibinomial(link="logit"), weights=wgt1, dat.z)$coef
    
    #logistfit<- eval(parse(text = paste("glm", "(", "w~Z[,-1]",",family = quasibinomial(link='", "logit", "'",")",")",sep = "")))
    #logistfit<-glm(pz~1+v, family=quasibinomial(link="logit"), weights=wgt2, dat.z)
    logistfit<- eval(parse(text = paste("glm", "(", "dat.z$pz~Z[,-1]",",family = quasibinomial(link='", "logit", "'","), weights=dat.z$wgt2",")",sep = "")))
    update_cureb <- logistfit$coef
    
    
    if(!is.null(beta)){
      update_beta<-coxph(Surv(dat.z$time, dat.z$del)~X, weights=dat.z$wgt2, subset=(dat.z$pz==1 & dat.z$wgt2!=0), method="breslow")$coef
      update_wcoxexp0=drop(exp(update_beta%*%t(X[dat.z$pz==1,,drop=FALSE]))*subset(dat.z, pz==1)$wgt2)
    }else{
      update_wcoxexp0=subset(dat.z, pz==1)$wgt2
    }
    
    test=data.frame(id=subset(dat.z, pz==1)$id, val=update_wcoxexp0)
    update_wcoxexp=aggregate(val~id, data=test, sum)$val#test$val#
    update_hazmass=curehaz(Time,Status,update_wcoxexp)$hazmass
    update_Hazf<-stepfun(event_point,c(0,cumsum(update_hazmass)))
    update_Hazv<-update_Hazf(dat.z$time)
    update_Hazv<-ifelse(dat.z$time>max(event_point), Inf, update_Hazv)
    update_Hazv<-ifelse(dat.z$time<min(event_point), 0, update_Hazv)
    update_s0=exp(-update_Hazv)
    
    df1h<-data.frame(id=useID, time=c(event_point,Time[Status==0]), lam0=c(update_hazmass, rep(1, sum(1-Status))))
    dfh<- right_join(df1h,df2h)
    dfh=dfh[order(dfh$id),]
    
    if(is.null(beta)){
      convergence<-sum(c(update_cureb-eta)^2)+sum((s0-update_s0)^2)
      
      eta <- update_cureb
      s0<-update_s0
      uncureprob <- matrix(exp((eta)%*%t(Z))/(1+exp((eta)%*%t(Z))),ncol=1)
      i <- i+1
      print(c(beta,eta))
      print(convergence)
      length(update_hazmass)
      
      em <- list(eta=eta, 
                 beta= beta,
                 lam0=update_hazmass, 
                 tk=event_point,
                 tau=convergence,
                 pdata=dat.z,
                 Z=Z)
      
    }else{
      convergence<-sum(c(update_cureb-eta,update_beta-beta)^2)+sum((s0-update_s0)^2)
      beta <- update_beta 
      eta <- update_cureb
      s0<-update_s0
      uncureprob <- matrix(exp((eta)%*%t(Z))/(1+exp((eta)%*%t(Z))),ncol=1)
      i <- i+1
      print(c(beta,eta))
      print(convergence)
      length(update_hazmass)
      
      em <- list(eta=eta, 
                 beta= beta,
                 lam0=update_hazmass, 
                 tk=event_point,
                 tau=convergence,
                 pdata=dat.z,
                 X=X,
                 Z=Z)
    }
    
    
  }
  
  return(em)
  
}

# estimation after phase II
em_coxph <-function(Status, Time, scen, dat.z, X, Z, V, est.H0,emmax,eps){     
  
  if(scen=="sc"){
    eta=c(est.H0$eta,0.01); gam=rep(0.01,ncol(V)); beta=c(est.H0$beta,0.01)
  }else if(scen=="sb"){
    eta=c(est.H0$eta); gam=rep(0.01,ncol(V)); beta=c(est.H0$beta,0.01)
  }else{
    eta=c(est.H0$eta,0.01); gam=rep(0.01,ncol(V)); beta=c(est.H0$beta)
  }
  
  event_point=sort(unique(subset(Time, Status==1)))
  
  dat.z$wgt<-dat.z$wgt1*dat.z$wgt2
  wcoxexp0 <- dat.z$wgt*drop(exp(beta%*%t(X)))
  test0=data.frame(id=subset(dat.z, pz==1)$id, val=wcoxexp0[dat.z$pz==1])
  wcoxexp=aggregate(val~id, data=test0, sum)$val
  hazmass=curehaz(Time,Status,wcoxexp)$hazmass
  Hazf<-stepfun(event_point,c(0,cumsum(hazmass)))
  Hazv<-Hazf(dat.z$time)
  Hazv<-ifelse(dat.z$time>max(event_point), Inf, Hazv)
  Hazv<-ifelse(dat.z$time<min(event_point), 0, Hazv)
  s0=exp(-Hazv)
  
  df1h<-data.frame(time=c(event_point,Time[Status==0]), lam0=c(hazmass, rep(1, sum(1-Status))))
  df2h<-dat.z
  dfh<- right_join(df1h,df2h)
  dfh=dfh[order(dfh$id),]
  
  convergence<- 1000;i <-1
  while (convergence > eps & i < emmax){  
    
    
    uncureprob <- matrix(exp((eta)%*%t(Z))/(1+exp((eta)%*%t(Z))),ncol=1)
    coxexp= drop(exp((beta)%*%t(X)))
    survival<-s0^(coxexp)
    
    ## E step 
    w1 <- wgt1_zeta_coxph.f(scen, dat.z$r, dat.z$px, dat.z$del, dat.z$time, X, Z, V, beta, eta, gam, dfh$lam0, s0)
    w2 <- dat.z$del+(1-dat.z$del)*(uncureprob*survival)/((1-uncureprob)+uncureprob*survival)
    
    dat.z$wgt1 <- w1
    dat.z$wgt2 <- ifelse(dat.z$pz==1, w2, 1-w2)
    dat.z$wgt <- dat.z$wgt1*dat.z$wgt2
    
    ## M step
    
    covXfit <- eval(parse(text = paste("glm", "(", "dat.z$px~V[,-1]",",family = quasibinomial(link='", "logit", "'","), weights=dat.z$wgt",")",sep = "")))
    update_gam=covXfit$coef
    
    logistfit<- eval(parse(text = paste("glm", "(", "dat.z$pz~Z[,-1]",",family = quasibinomial(link='", "logit", "'","), weights=dat.z$wgt",")",sep = "")))
    update_cureb <- logistfit$coef
    
    ## use dat.z for coxph() and baseline survival dist. estimation 
    update_beta<-coxph(Surv(dat.z$time, dat.z$del)~X, weights=dat.z$wgt, subset=(dat.z$pz==1 & dat.z$wgt!=0), method="breslow")$coef
    update_wcoxexp0=drop(exp(update_beta%*%t(X[dat.z$pz==1,,drop=FALSE]))*subset(dat.z, pz==1)$wgt)
    test=data.frame(id=subset(dat.z, pz==1)$id, val=update_wcoxexp0)
    update_wcoxexp=aggregate(val~id, data=test, sum)$val
    update_hazmass=curehaz(Time,Status,update_wcoxexp)$hazmass
    update_Hazf<-stepfun(event_point,c(0,cumsum(update_hazmass)))
    update_Hazv<-update_Hazf(dat.z$time)
    update_Hazv<-ifelse(dat.z$time>max(event_point), Inf, update_Hazv)
    update_Hazv<-ifelse(dat.z$time<min(event_point), 0, update_Hazv)
    update_s0=exp(-update_Hazv)
    
    df1h<-data.frame(time=c(event_point,Time[Status==0]), lam0=c(update_hazmass, rep(1, sum(1-Status))))
    dfh<- right_join(df1h,df2h)
    dfh=dfh[order(dfh$id),]
    
    convergence<-sum(c(update_cureb-eta,update_beta-beta)^2)+sum((s0-update_s0)^2)
    eta <- update_cureb
    beta <- update_beta 
    gam <- update_gam
    s0<-update_s0
    uncureprob <- matrix(exp((eta)%*%t(Z))/(1+exp((eta)%*%t(Z))),ncol=1)
    i <- i+1
    print(c(beta,eta,gam))
    print(convergence)
    length(update_hazmass)
  }
  em <- list(eta=eta, 
             beta= beta,
             gam=gam,
             lam0=update_hazmass, 
             tk=event_point,
             tau=convergence,
             pdata=dat.z)
}

