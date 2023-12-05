
# ------------------------------------------------------------------------------------------
# This is a script for one simulation study. 
# A for-loop is used for repeated simulations. 
# In each simulation run, 3 steps are included: 
#  1) data generation
#  2) two-phase design implementation
#  3) estimation and inference
# ------------------------------------------------------------------------------------------

library(MASS)
library(parallel)
library(dplyr)
library(survival)
library(cubature)
source("utility/design.r") # functions for phase II designs
source("utility/design_withV.r") # functions for phase II designs
source("utility/em_coxph.r") # functions for esitmation via EM
source("utility/avar_louis.r") # functions for asymptotic variance computation for EM
source("utility/opt_str_design_gau.r") # optimal stratified designs, not feasible in practice 

####### Test Data configuration #######

# phase I sample size 
n=5000 
# true regression coefs. for T|X,Z=1
be10=0.0; be20=0.22;
# true parameters for censoring time
sh0 = 1.2; sc0=0.5 
# true parameters for X1|X2
ga10 = log(0.3/0.7); ga20=0 
# true regression coefs. for Z|X
et10 = log(0.3/0.7); et20=0.0; et30=0.11 

# phase I stratifying factor: discretized observed times and event status
phI.design.factor=c("cbar","del"); 
# size of the phase II sample
n2_samp=1000; 

# balanced design used in phase II
phII.design="bal"; 

# scen: regression scenarios; 
#  scen="sa" if the expensive biomarker only adjusted in the susceptibility model;
#  scen="sb" if the expensive biomarker only adjusted in the failure time model;
#  scen="sc" if the expensive biomarker adjusted in both models.
scen="sc"

# covXname.H0: labels of inexpensive covariate adjusted in the failure time model;
# covZname.H0: labels of inexpensive covariate adjusted in the susceptibility model;
# covVname: labels of inexpensive covariate adjusted in the covariate model;
covXname.H0=c("v"); covZname.H0=c("v"); covVname=c("v")

if(scen=="sc"){
  
  covXname=c(covXname.H0,"px") # complete labels of covariates adjusted in the failure time model;
  covZname=c(covZname.H0, "px") # complete labels of covariates adjusted in the susceptibility model;
  
}else if(scen=="sb"){
  
  covXname=c(covXname.H0,"px")
  covZname=c(covZname.H0)
  
}else{
  
  covXname=c(covXname.H0)
  covZname=c(covZname.H0, "px")
  
}


####### prepare a .dat file that record the simulation results ######

phI.strat.style = paste0(phI.design.factor, collapse = "_")
param.set = c(ga10, ga20, be10, be20, et10, et20, et30, sh0, sc0)
name.param.set = c("gam0", "gam1",  "beta1",  "beta2", "eta1", "eta2", "eta3", "wei.shape", "wei.scale")
names(param.set) = name.param.set

file1.out <- paste(n2_samp, scen,"coxph.", phI.strat.style, phII.design, 
                   paste(append(paste(name.param.set, param.set, sep="_"), ".dat"), collapse=""), sep="")

est.para <- c(paste0("eta_", c(0, covZname)), 
              paste0("gam_", c(0, covVname)),
              paste0("beta_", c(covXname)))

group.mat <- expand.grid(1:3, 0:1)
ref.group <- paste0(group.mat[,1], group.mat[,2]) # labels of the phase I strata: "10" "20" "30" "11" "21" "31"

label.tx <- c("nsim", est.para, paste0("ASE.", est.para), paste0("p.", ref.group), paste0("g", ref.group))

if ( file.exists( file1.out ) ) { unlink( file1.out ) }

cat(as.character(label.tx), sep=" ", "\n", append=T, file=file1.out)


###### start the simulation ######

nsim=100 # number of simulation runs
iter = 1
for (k in 1:(2*nsim)){

  if(iter>nsim){break}
  
  # --- 1) Data Generation --- #
  
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
  
  # --- 2) Two-phase Design --- #
  
  # Phase I stratification
  phI_strat = phaesI_strat.f(pcuts=c(1/3,2/3), dat, design.factor=phI.design.factor)
  
  # fit a reduced model under H0: beta1=0 in phase I
  #' covXname: labels of inexpensive covariate adjusted in the failure time model;
  #' covZname: labels of inexpensive covariate adjusted in the susceptibility model;
  #' emmax: max iteration runs for the EM algorithm
  #' eps: convergence level for the EM algorithm
  #' eta, beta: initial values
  est.H0=em_coxph_H0(dat, covXname=c("v"), covZname=c("v"), 
                     eta=rep(0.01, 2), beta=rep(0.01, 1), 
                     emmax=100, eps=1e-03)
  
  # Phase II selection, multiple choices are available
  
  if(phII.design %in% c("RSD-Smu1", "RSD-Smu2", paste0("BRSD-eff-seq", seq(0,1, by=0.1)))){
    
    # RSD-Smu1: univeriate residual-dependent sampling using Smu1
    # RSD-Smu2: univeriate residual-dependent sampling using Smu2
    # BRSD-eff-seq: bivariate residual-dependent sampling using standardized Smu1 and Smu2
    # w: a value in [0,1], the proportion of a phase II sample was selected in wave 1;
    # 1-w: a value in [0,1], the proportion of a phase II sample was selected in wave 2;
    
    # residual-dependent sampling 
    phII_sel = designPhII_resid_coxph_withV(phaseI_strat=phI_strat, n2samp=n2_samp, 
                                            design=phII.design, w=0.5, 
                                            covXname=c("v"), covZname=c("v"), 
                                            est.H0$beta, est.H0$eta, est.H0$lam0, est.H0$tk)$sel
    
  }else if(phII.design=="opt1"){
    
    # opt1: optimal stratified sampling for beta1, biomarker effect on the failure time process
    
    phII_sel <- try(opt1_ODS_asym_AY.f(scen, grad=1e-06,n2_samp, phI_strat, ini=rep(log((n2_samp/n)/(1-n2_samp/n)), 6), 
                                       c(log(2), 1.1, be10, be20, et10, et20, et30, ga10, ga20)
                                       , sc0, sh0, 0.44))
    
    if(isTRUE(class(phII_sel)=="try-error")) { next }
    
  }else if(phII.design=="opt2"){
    
    # opt2: optimal stratified sampling for eta1, biomarker effect on the susceptibility 
    
    phII_sel <- try(opt2_ODS_asym_AY.f(scen, grad=1e-06,n2_samp, phI_strat, ini=rep(log((n2_samp/n)/(1-n2_samp/n)), 6), 
                                       c(log(2), 1.1, be10, be20, et10, et20, et30, ga10, ga20)
                                       , sc0, sh0, 0.44))
    
    if(isTRUE(class(phII_sel)=="try-error")) { next }
  }else if(phII.design=="optA"){
    
    # optA: optimal stratified sampling based on A-optimality
    
    phII_sel <- try(optA_ODS_asym_AY.f(grad=1e-06,n2_samp, phI_strat, ini=rep(log((n2_samp/n)/(1-n2_samp/n)), 6), 
                                       c(log(2), 1.1, be10, be20, et10, et20, et30, ga10, ga20)
                                       , sc0, sh0, 0.44))
    
    if(isTRUE(class(phII_sel)=="try-error")) { next }
    
  }else if(phII.design=="optD"){
    
    # optA: optimal stratified sampling based on D-optimality
    
    phII_sel <- try(optD_ODS_asym_AY.f(grad=1e-06,n2_samp, phI_strat, ini=rep(log((n2_samp/n)/(1-n2_samp/n)), 6), 
                                       c(log(2), 1.1, be10, be20, et10, et20, et30, ga10, ga20)
                                       , sc0, sh0, 0.44))
    
    if(isTRUE(class(phII_sel)=="try-error")) { next }
  }else{
    
    # stratified-dependent sampling
    # phII.design = "bal" or "srs"  
    phII_sel = design_strat.f(n2_samp, phII.design, phI_strat)
  }    
  
  # incomplete data set with partially missing covariate x
  dt.phII = dat
  dt.phII$x[-phII_sel$s_id] <- NA
  dt.phII$r <- 1
  dt.phII$r[-phII_sel$s_id] <- 0
  
  
  ###### 3) Estimation and Inference ######
  
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
  
  # prepare the covariate matrix: 
  # X for the failure time model; 
  # Z for the susceptibility model; 
  # V for the covariate model
  X=as.matrix(dat.z[,covXname])
  Z=as.matrix(cbind(rep(1,nrow(dat.z)), dat.z[,covZname]))
  V=as.matrix(cbind(rep(1,nrow(dat.z)), dat.z[,covVname]))
  
  # estimation after design: 
  #' scen: regression scenarios; 
  #'  scen="sa" if the expensive biomarker only adjusted in the susceptibility model;
  #'  scen="sb" if the expensive biomarker only adjusted in the failure time model;
  #'  scen="sc" if the expensive biomarker adjusted in both models.
  #' emmax: max iteration runs for the EM algorithm
  #' eps: convergence level for the EM algorithm
  est <- em_coxph(Status, Time, scen, dat.z, X, Z, V, est.H0,
                  emmax=200,eps=1e-06)
  
  # asymptotic variance matrix  
  asvar <- CureCox.Louis.varMatrix(est$lam0, est$tk, est$pdata$id, est$pdata$pz, est$pdata$del, est$pdata$time, est$pdata$px, 
                                   X, Z, V, est$beta, est$eta, est$gam, est$pdata$wgt1, est$pdata$wgt2)
  # estimated asymptotic standard error
  ase <- sqrt(diag(asvar))
  
  
# --- record the results --- #
  
  cat(k, est$eta, est$gam, est$beta, ase, 
      phII_sel$s_prob, 
      phI_strat$strata.n,
      seq=" ", "\n", append = T, file=file1.out)
  
  print(iter)
  iter = iter + 1
  
}

