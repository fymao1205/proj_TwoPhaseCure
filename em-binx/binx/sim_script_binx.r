library(MASS)
library(parallel)
library(dplyr)
#library(Rcpp)
library(survival)
library(cubature)
source("/Users/maof3/Documents/research/2022/2phase_cure/withV/adj_TZ/clean/designfct/design.r")
source("/Users/maof3/Documents/research/2022/2phase_cure/withV/adj_TZ/clean/designfct/design_withV.r")
source("/Users/maof3/Documents/research/2022/2phase_cure/withV/adj_TZ/clean/binx/em_coxph.r")
source("/Users/maof3/Documents/research/2022/2phase_cure/withV/adj_TZ/clean/binx/avar_louis.r")
source("/Users/maof3/Documents/research/2022/2phase_cure/withV/adj_TZ/coxph/opt_str_design_gau.r") # not feasible in practice 


script_coxph_bin_withV.f <- function(nsim, n2_samp, scen, emmax=200,eps=1e-06,
                                     covXname.H0=c("v"), covZname.H0=c("v"), covVname=c("v"),
                                     phI.design.factor, phII.design, w, pcuts, n, 
                                   ga10, ga20, be10, be20, et10, et20, et30, sh0, sc0){
  
  phI.strat.style = paste0(phI.design.factor, collapse = "_")
  
  param.set = c(ga10, ga20, be10, be20, et10, et20, et30, sh0, sc0)
  name.param.set = c("gam0", "gam1",  "beta1",  "beta2", "eta1", "eta2", "eta3", "wei.shape", "wei.scale")
  
  names(param.set) = name.param.set
  
  if(phI.strat.style=="cbar"){
    ref.group <- 1:(length(pcuts)+1)
    
  }
  
  if(phI.strat.style=="del"){
    
    ref.group <- 0:1
    
  }
  
  if(phI.strat.style=="cbar_del"){
    
    group.mat <- expand.grid(1:3, 0:1)
    ref.group <- paste0(group.mat[,1], group.mat[,2])
  }
  
  if(scen=="sc"){
    
    covXname=c(covXname.H0,"px")
    covZname=c(covZname.H0, "px")
    
  }else if(scen=="sb"){
    
    covXname=c(covXname.H0,"px")
    covZname=c(covZname.H0)
    
  }else{
    
    covXname=c(covXname.H0)
    covZname=c(covZname.H0, "px")
    
  }
  
  if(!is.null(covXname)){
    est.para <- c(paste0("eta_", c(0, covZname)), 
                  paste0("gam_", c(0, covVname)),
                  paste0("beta_", c(covXname)))
  }else{
    est.para <- c(paste0("eta_", c(0, covZname)), 
                  paste0("gam_", c(0, covVname)))
  }
  
  
  file1.out <- paste(n2_samp, scen,"coxph.", phI.strat.style, phII.design, paste(append(paste(name.param.set, param.set, sep="_"), ".dat"), collapse=""), sep="")
  
  label.tx <- c("nsim", est.para, 
    paste0("ASE.", est.para
    ),# paste0("score.", est.para),
    paste0("p.", ref.group), paste0("g", ref.group))
  
  if ( file.exists( file1.out ) ) { unlink( file1.out ) }
  
  cat(as.character(label.tx), sep=" ", "\n", append=T, file=file1.out)
  
  iter = 1
  for (k in 1:(2*nsim)){
    #while(iter <= nsim){
    if(iter>nsim){break}
    
    ###########################
    ##### Data Generation #####
    ###########################
    set.seed(k)
    
    # inexpensive covariate v
    v <- rbinom(n, 1, 0.44)
    
    # expensive covariate x
    prob.x <- expit.f(ga10+ga20*v)
    x <- rbinom(n, 1, prob.x)
    
    # event time 
    u.sim <- runif(n, 0, 1)
    t <- ((-log(u.sim))*exp(-v*be20-x*be10))^(1/1.1)*0.5
    
    # susceptibility indicator 
    prob.z <- expit.f(et10+et20*x+et30*v)
    z <- rbinom(n, 1, prob.z)
    
    # generate right censoring time
    c <- pmin(rweibull(n, sh0, sc0),10)
    
    # event status: yes-1; no-0
    del <- ifelse(t<c & z==1, 1, 0)
    
    # observed time 
    time <- ifelse(del==1, t, c)
    
    # complete data 
    dat <- cbind(del, time, x,v)
    dat <- as.data.frame(dat)
    dat$id <- 1:n
    
    ###########################
    ##### Two-phase Design #####
    ###########################
    
    # Phase I stratification
    phI_strat = phaesI_strat.f(pcuts=c(1/3,2/3), dat, design.factor=phI.design.factor)
    
    # Phase II selection 
    if(is.null(covXname.H0)){
      beta.ini=NULL
    }else{
      beta.ini=rep(0.01, length(covXname.H0))
    }
    
    # fit a reduced model under H0: beta1=0 in phase I
    est.H0=em_coxph_H0(dat, covXname.H0, covZname.H0, 
                       rep(0.01, length(covZname.H0)+1), beta.ini, emmax=100, eps=1e-03)
    
    if(phII.design %in% c("RSD-Smu1", "RSD-Smu2", paste0("BRSD-eff-seq",w))){
      
      # residual-dependent sampling 
      phII_sel = designPhII_resid_coxph_withV(phaseI_strat=phI_strat, n2samp=n2_samp, 
                                              design=phII.design, w, 
                                              covXname.H0, covZname.H0, 
                                              est.H0$beta, est.H0$eta, est.H0$lam0, est.H0$tk)$sel
      
    }else if(phII.design=="opt1"){
      
      phII_sel <- try(opt1_ODS_asym_AY.f(scen, grad=1e-06,n2_samp, phI_strat, ini=rep(log((n2_samp/n)/(1-n2_samp/n)), 6), 
                                        c(log(2), 1.1, be10, be20, et10, et20, et30, ga10, ga20)
                                        , sc0, sh0, 0.44))
      
      if(isTRUE(class(phII_sel)=="try-error")) { next }
    }else if(phII.design=="opt2"){
      
      phII_sel <- try(opt2_ODS_asym_AY.f(scen, grad=1e-06,n2_samp, phI_strat, ini=rep(log((n2_samp/n)/(1-n2_samp/n)), 6), 
                                         c(log(2), 1.1, be10, be20, et10, et20, et30, ga10, ga20)
                                         , sc0, sh0, 0.44))
      
      if(isTRUE(class(phII_sel)=="try-error")) { next }
    }else if(phII.design=="optA"){
      
      phII_sel <- try(optA_ODS_asym_AY.f(grad=1e-06,n2_samp, phI_strat, ini=rep(log((n2_samp/n)/(1-n2_samp/n)), 6), 
                                         c(log(2), 1.1, be10, be20, et10, et20, et30, ga10, ga20)
                                         , sc0, sh0, 0.44))
      
      if(isTRUE(class(phII_sel)=="try-error")) { next }
    }else if(phII.design=="optD"){
      
      phII_sel <- try(optD_ODS_asym_AY.f(grad=1e-06,n2_samp, phI_strat, ini=rep(log((n2_samp/n)/(1-n2_samp/n)), 6), 
                                         c(log(2), 1.1, be10, be20, et10, et20, et30, ga10, ga20)
                                         , sc0, sh0, 0.44))
      
      if(isTRUE(class(phII_sel)=="try-error")) { next }
    }else{
      # stratified-dependent sampling
      phII_sel = design_strat.f(n2_samp, phII.design, phI_strat)
    }    
    
    # incomplete data set with partially missing covariate x
    dt.phII = dat
    dt.phII$x[-phII_sel$s_id] <- NA
    dt.phII$r <- 1
    dt.phII$r[-phII_sel$s_id] <- 0
    
    
    ###########################
    ##### Estimation and Inference #####
    ###########################
    
    #--- creating pseudo data set ---#
    dat.orig=dt.phII
    Status=dat.orig$del
    Time=dat.orig$time
    dat.orig$id=1:nrow(dat.orig)
    
    # expand the dat.orig for missing covariate x
    px <- c( subset(dat.orig, r==1)$x, rep(c(1,0), each=sum(1-dat.orig$r)) )
    xdata <- cbind(rbind(subset(dat.orig, r==1),
                         subset(dat.orig, r==0),
                         subset(dat.orig, r==0)), px)
    dat.x <- as.data.frame(xdata)
    
    # expand the dat.x for susceptibility model 
    pz <- c(rep(1, dim(dat.x)[1]), rep(0, sum(1-dat.x$del)))
    dat.z <- cbind(rbind(subset(dat.x, del==1),
                         subset(dat.x, del==0),
                         subset(dat.x, del==0)),pz)
    dat.z <- dat.z[order(dat.z$id),]
    
    # weights for missing covariates
    p1=mean(dat.orig$x,na.rm=T)
    w1<-ifelse(dat.z$px==1, p1, 1-p1)
    dat.z$wgt1<-ifelse(dat.z$r==1,1,w1)
    
    # weights for latent susceptibility
    w2<-dat.z$del
    dat.z$wgt2<-ifelse(dat.z$pz==1,w2,1-w2)
    
    # prepare the covariate matrix: X for the failure time model; Z for the susceptibility model; V for the covariate model
    X=as.matrix(dat.z[,covXname])#c("v","px")
    Z=as.matrix(cbind(rep(1,nrow(dat.z)), dat.z[,covZname]))
    V=as.matrix(cbind(rep(1,nrow(dat.z)), dat.z[,covVname]))
    
    # estimation after design 
    est <- em_coxph(Status, Time, scen, dat.z, X, Z, V, est.H0,
                    emmax,eps)
    
    # asymptotic variance computation  
    asvar <- CureCox.Louis.varMatrix(est$lam0, est$tk, est$pdata$id, est$pdata$pz, est$pdata$del, est$pdata$time, est$pdata$px, 
                                     X, Z, V, est$beta, est$eta, est$gam, est$pdata$wgt1, est$pdata$wgt2)
    ase <- sqrt(diag(asvar))
    
    
    ###########################
    ##### record the results #####
    ###########################
    
    cat(k, est$eta, est$gam, est$beta, ase, 
        phII_sel$s_prob, 
        phI_strat$strata.n,
        seq=" ", "\n", append = T, file=file1.out)
    
    print(iter)
    iter = iter + 1
    
  }
}


n=5000
be10=0.0; 
be20=0.22; 
ga10 = log(0.3/0.7); ga20=0
sh0 = 1.2; sc0=0.5 # nonrare
et10 = log(0.3/0.7)#log(0.6/0.4); #
et20=0.0; et30=0.11
phI.design.factor=c("cbar","del")
#phII.design="bal"
n2_samp=1000; 
nsim=100

mclapply(seq(0,1,by=0.2), function(w){
  script_coxph_bin_withV.f(nsim=100, n2_samp=n2_samp, scen="sc", emmax=1000,eps=1e-06,
                           covXname.H0=c("v"), covZname.H0=c("v"), covVname=c("v"),
                           phI.design.factor=phI.design.factor, phII.design=paste0("BRSD-eff-seq",w), w, pcuts=c(1/3,2/3),
                           n=n, ga10=ga10, ga20=ga20, be10=be10, be20=be20, et10=et10, et20=et20, et30=et30, sh0=sh0, sc0=sc0)
  
}, mc.cores=6)

script_coxph_bin_withV.f(nsim=100, n2_samp=n2_samp, scen="sc", emmax=1000,eps=1e-06,
                         covXname.H0=c("v"), covZname.H0=c("v"), covVname=c("v"),
                         phI.design.factor=phI.design.factor, phII.design=paste0("BRSD-mtp"), w=NULL, pcuts=c(1/3,2/3),
                         n=n, ga10=ga10, ga20=ga20, be10=be10, be20=be20, et10=et10, et20=et20, et30=et30, sh0=sh0, sc0=sc0)

script_coxph_bin_withV.f(nsim=100, n2_samp=n2_samp, scen="sc", emmax=1000,eps=1e-06,
                         covXname.H0=c("v"), covZname.H0=c("v"), covVname=c("v"),
                         phI.design.factor=phI.design.factor, phII.design=paste0("BRSD-add",1), w=1, pcuts=c(1/3,2/3),
                         n=n, ga10=ga10, ga20=ga20, be10=be10, be20=be20, et10=et10, et20=et20, et30=et30, sh0=sh0, sc0=sc0)


script_coxph_bin_withV.f(nsim=100, n2_samp=n2_samp, scen="sc", emmax=1000,eps=1e-06,
                         covXname.H0=c("v"), covZname.H0=c("v"), covVname=c("v"),
                         phI.design.factor=phI.design.factor, phII.design="opt2", w=NULL, pcuts=c(1/3,2/3),
                         n=n, ga10=ga10, ga20=ga20, be10=be10, be20=be20, et10=et10, et20=et20, et30=et30, sh0=sh0, sc0=sc0)

script_coxph_bin_withV.f(nsim=100, n2_samp=n2_samp, scen="sc", emmax=1000,eps=1e-06,
                         covXname.H0=c("v"), covZname.H0=c("v"), covVname=c("v"),
                         phI.design.factor=phI.design.factor, phII.design="optA", w=NULL, pcuts=c(1/3,2/3),
                         n=n, ga10=ga10, ga20=ga20, be10=be10, be20=be20, et10=et10, et20=et20, et30=et30, sh0=sh0, sc0=sc0)
script_coxph_bin_withV.f(nsim=100, n2_samp=n2_samp, scen="sc", emmax=1000,eps=1e-06,
                         covXname.H0=c("v"), covZname.H0=c("v"), covVname=c("v"),
                         phI.design.factor=phI.design.factor, phII.design="optD", w=NULL, pcuts=c(1/3,2/3),
                         n=n, ga10=ga10, ga20=ga20, be10=be10, be20=be20, et10=et10, et20=et20, et30=et30, sh0=sh0, sc0=sc0)

script_coxph_bin_withV.f(nsim=100, n2_samp=n2_samp, scen="sa", emmax=200,eps=1e-06,
                         covXname.H0=c("v"), covZname.H0=c("v"), covVname=c("v"),
                         phI.design.factor=phI.design.factor, phII.design="RSD-Smu1", pcuts=c(1/3,2/3),
                         n=n, ga10=ga10, ga20=ga20, be10=be10, be20=be20, et10=et10, et20=et20, et30=et30, sh0=sh0, sc0=sc0)

test <- read.table(file=file.choose(), header=1)

mclapply(c("srs","optn-Smu4"), function(x){
  script_coxph_bin_withV.f(nsim=100, phI.design.factor=phI.design.factor, phII.design=x, pcuts=c(1/3,2/3),
                         n=n, n2_samp=n2_samp, 
                         ga10=ga10, ga20=ga20, be10=be10, be20=be20, et10=et10, et20=et20, et30=et30, sh0=sh0, sc0=sc0)
}, mc.cores=2)

test=read.table(file=file.choose(), header=1)
colMeans(test)[2:8]
c(et10,et30,et20,ga10,ga20,be20,be10)
coverage(test[,2:8], test[,8+(1:7)], c(et10,et30,et20,ga10,ga20,be20,be10), 0.95)
