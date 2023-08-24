
library(Rcpp)
library(MASS)
library(dplyr)
library(survival)
sourceCpp("commonf.cpp")
sourceCpp("loglik_pwc.cpp")
source("rldt_est_pwc.r")
source("design_pwc.r")
source("design_withV_pwc.r")

#######################################
###
### load and prepare the .RData
###
#######################################
load("~/Documents/research/2022/2phase_cure/real_data/PLCO/plco_imputed_06112022.RData")
cov.names=c("ID", "incidence.age","incidence.days", "age","case", "female", "family.history.lung.cancer", "bron", "bmi", "race", "qtyears", "pkyr.cat", "EUR_PRS", "Meta_PRS")
dat0 <- plco[,cov.names, drop=FALSE] #subset(plco, !is.na(Meta_PRS))
na.ind=apply(dat0, 1, function(x)(any(is.na(x))))
sub_dat0 <- subset(dat0, !na.ind)

sub_dat0$yr_fup=sub_dat0$incidence.days/365.25

plco_dt=data.frame(#id=1:nrow(sub_dat0), 
  time=sub_dat0$yr_fup,
  age=sub_dat0$age, 
  logage=log(sub_dat0$age),
  #time=sub_dat0$incidence.age,
  del=sub_dat0$case, 
  gender=sub_dat0$female, 
  bmi=(sub_dat0$bmi>30),
  race=(sub_dat0$race==0),
  bron=sub_dat0$bron,
  fmh=sub_dat0$family.history.lung.cancer,
  qtyr=(sub_dat0$qtyears),
  qtyr0=ifelse(sub_dat0$qtyears==0, 1, 0), 
  qtyr1=ifelse(sub_dat0$qtyears>0 & sub_dat0$qtyears<=15, 1, 0),
  qtyr2=ifelse(sub_dat0$qtyears>15, 1, 0),
  pkyr=(sub_dat0$pkyr.cat),
  logqtyr=ifelse(sub_dat0$qtyears==0, 0, log(sub_dat0$qtyears)), 
  logpkyr=ifelse(sub_dat0$pkyr.cat==0, 0, log(sub_dat0$pkyr.cat)),
  pkyr0=ifelse(sub_dat0$pkyr.cat==0, 1, 0), 
  pkyr1=ifelse(sub_dat0$pkyr.cat>0 & sub_dat0$pkyr.cat<20, 1, 0),
  pkyr2=ifelse(sub_dat0$pkyr.cat>=20 & sub_dat0$pkyr.cat<40, 1, 0),
  #pkyr3=ifelse(sub_dat0$pkyr.cat>=40 & sub_dat0$pkyr.cat<50, 1, 0),
  pkyr4=ifelse(sub_dat0$pkyr.cat>=40, 1, 0),
  x=sub_dat0$Meta_PRS,#(sub_dat0$Meta_PRS-min(sub_dat0$Meta_PRS)), #ifelse(sub_dat0$qtyears==0, 1, 0), 
  v=sub_dat0$female)

dt0=subset(plco_dt, pkyr0 ==0 )#& pkyr1==0 & qtyr2==0) #& qtyr0==1)#pkyr0==0 & qtyr2==0) #pkyr2==1 & qtyr0==1) #
indices = rbinom(nrow(dt0), 1, 1)
dt=dt0[indices==1,,drop=FALSE]
dt$id=1:nrow(dt)

### Phase II selection 
covXname.H0=c("gender", "bron", "fmh") 
covZname.H0=#c("gender", "bron", "fmh", "qtyr0", "pkyr2","pkyr4");
covVname=NULL


#' nsim: number of replications
#' n2samp: user-specified size of the phase II sub-sample
#' phII.design={"srs", "bal", "RSD-Smu1", "RSD-Smu2", paste0("BRSD-eff-seq",w)}:
#' "srs": simple random sampling
#' "bal": balanced sampling
#' "RSD-Smu1": univariate residual dependent sampling using Smu1
#' "RSD-Smu2": univariate residual dependent sampling using Smu2
#' paste0("BRSD-eff-seq",w): our proposed bivariate residual dependent sampling using both Smu1 and Smu2 (transformed)
#' w: an user-specified weight, a value between 0 and 1. select (w x n2samp) subjects for eta1 estimation and ((1-w) x n2samp) subjects for neta1
#' scen={"sa", "sb", "sc"}: 
#' sa=scenario a: adjust expensive covariate x only in the susceptibility Z model
#' sb=scenario b: adjust expensive covariate x only in the failure time T model
#' sc=scenario c: adjust expensive covariate x in both T and Z models
#' pwc.cuts: percentiles used to specify the cutpoints for the PWC model, e.g. pwc.cuts=c(1/3,2/3) the tertiles
#' covXname.H0: user-specfied inexpensive covariates adjusted in the T model
#' covZname.H0: user-specfied inexpensive covariates adjusted in the Z model
#' covVname=NULL: user-specfied inexpensive covariates adjusted in the covariate model X1|X2; 
#' here X1 is set to be normally distributed and X2=NULL
rldt_script_pwc.f <- function(nsim=1, n2samp, phII.design, w, scen="sc", pwc.cuts, 
                              covXname.H0, covZname.H0, covVname=NULL){
  
  if(is.null(covXname.H0)){
    X0=NULL
  }else{
    X0=as.matrix(dt[,covXname.H0])
  }
  
  if(is.null(covZname.H0)){
    Z0=as.matrix(cbind(rep(1,nrow(dt))))
  }else{
    Z0=as.matrix(cbind(rep(1,nrow(dt)), dt[,covZname.H0])) 
  }
  
  V=as.matrix(cbind(rep(1,nrow(dt)), dt[,covVname]))
  
  #######################################
  ###
  ### Phase I stratification
  ###
  #######################################
  phI.design.factor=c("cbar","del")
  phI.strat.style = paste0(phI.design.factor, collapse = "_")
  phI_strat = phaesI_strat.f(pcuts=c(1/3,2/3), dt, design.factor=phI.design.factor)
  
  #######################################
  ### estimation using a simplified null model for 2 purposes: 
  ### 1): to compute the score-type residuals
  ### 2): to be used as an initial for the main estimation 
  #######################################
  brks=quantile(dt$time[dt$del==1], pwc.cuts) # specify the cutpoints (using the empirical tertiles of the observed failure time) for the PWC constand hazard model
  est.H0 <- rldt_est_pwc_H0.f(dt$del, dt$time, X0, Z0, brks)
  
  if(phI.strat.style=="cbar_del"){
    group.mat <- expand.grid(1:3, 0:1)
    ref.group <- paste0(group.mat[,1], group.mat[,2])
  }
  
  if(scen=="sc"){# scenario c: adjust expensive covariate x in both models
    
    # add the label for the missing expensive covariate 
    covXname=c(covXname.H0,"px") # the complete covariate list adjusted in the failure time T model
    covZname=c(covZname.H0, "px") # the complete covariate list adjusted in the susceptibility Z model
    
    # form the initial value of parameters
    ini=c(est.H0$est[1:(length(brks)+1)], 0.01, 0.01, 
          est.H0$est[-(1:(length(brks)+1))], 0.01, 0.01)
    
    est.para <- c(paste0("logalp_", 1:(length(brks)+1)), 
                  c("eta1", "beta1"), 
                  paste0("eta_", c(0, covZname.H0)), 
                  paste0("gam_", c(0, 1)),
                  paste0("beta_", c(covXname.H0)))
    
  }else if(scen=="sb"){
    
    covXname=c(covXname.H0,"px")
    covZname=c(covZname.H0)
    
    ini=c(est.H0$est[1:(length(brks)+1)], 0.01, 
          est.H0$est[-(1:(length(brks)+1))], 0.01, 0.01)
    
    est.para <- c(paste0("logalp_", 1:(length(brks)+1)), 
                  c("beta1"), 
                  paste0("eta_", c(0, covZname.H0)), 
                  paste0("gam_", c(0, 1)),
                  paste0("beta_", c(covXname.H0)))
    
  }else{# scenario a: adjust expensive covariate x only in the Z model
    
    covXname=c(covXname.H0)
    covZname=c(covZname.H0, "px")
    
    ini=c(est.H0$est[1:(length(brks)+1)], 0.01, 
          est.H0$est[-(1:(length(brks)+1))], 0.01, 0.01)
    
    est.para <- c(paste0("logalp_", 1:(length(brks)+1)), 
                  c("eta1"), 
                  paste0("eta_", c(0, covZname.H0)), 
                  paste0("gam_", c(0, 1)),
                  paste0("beta_", c(covXname.H0)))
  }
  
  file1.out <- paste(n2samp, scen,"pwc.", paste0(round(brks),collapse = "_"), phI.strat.style, phII.design, 
                     paste(append(paste(c(covXname.H0, covZname.H0, covVname), collapse="_"), ".dat"), collapse=""), 
                     sep="")
  
  label.tx <- c("nsim", est.para, 
                paste0("ASE.", est.para), paste0("hessASE.", est.para), 
                paste0("score.", est.para),
                paste0("p.", ref.group), paste0("g", ref.group))
  
  if ( file.exists( file1.out ) ) { unlink( file1.out ) }
  cat(as.character(label.tx), sep=" ", "\n", append=T, file=file1.out)
  
  iter=1
  for(k in 1:nsim){
    
    if(iter>nsim){break}
    set.seed(k)
    
    if(phII.design %in% c("RSD-Smu1", "RSD-Smu2", paste0("BRSD-eff-seq",w))){
      
      phII_sel = designPhII_resid_pwc_noV(phaseI_strat=phI_strat, n2samp, 
                                          design=phII.design, w, 
                                          covXname.H0, covZname.H0, 
                                          est.H0$est,brks)#$sel
      
    }else{
      phII_sel = design_strat.f(n2samp, phII.design, phI_strat)
    }
    
    dt.phII = dt
    dt.phII$x[-phII_sel$s_id] <- NA
    dt.phII$r <- 1
    dt.phII$r[-phII_sel$s_id] <- 0
    dt.phII$pz<-dt$del
    
    
    est <-  rldt_est_pwc.f(ini=ini, 
                           scen,
                              dt.phII$del, dt.phII$time, dt.phII$x, X0, Z0, brks)
    print(sqrt(diag(ginv(est$hess))))
    ase <- ase_pwc_norm.f(grad=1e-08, scen, est, dt.phII$del, dt.phII$time, dt.phII$x, X0, Z0, brks)
    
    print(ase$ase)
    
    
    cat(k, est$par, 
        ase$ase, sqrt(diag(ginv(est$hess))), 
        ase$score,
        phII_sel$s_prob, 
        phI_strat$strata.n,
        seq=" ", "\n", append = T, file=file1.out)
    
    print(iter)
    iter = iter + 1
  }
}

# rldt_script_pwc.f(1, n2samp=3000, phII.design="RSD-Smu1", w=NULL, scen="sc", est.H0, brks, covXname.H0, covZname.H0, covVname)
# rldt_script_pwc.f(1, n2samp=5000, phII.design="BRSD-eff-seq0.5", w=0.5, scen="sc", est.H0, brks, covXname.H0, covZname.H0, covVname)
# rldt_script_pwc.f(1, n2samp=5000, phII.design="bal", w=NULL, scen="sc", est.H0, brks, covXname.H0, covZname.H0, covVname)
rldt_script_pwc.f(1, n2samp=nrow(dt), phII.design="srs", w=NULL, scen="sc", c(1/3,2/3), covXname.H0, covZname.H0, covVname)


