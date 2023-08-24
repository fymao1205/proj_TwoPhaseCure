
# for residual-dependent sampling 

wgt2_p_coxph_withV.f <- function(z, del, t, X, Z, beta.H0, eta.H0, lam_vec, tk_vec){
  
  covZ = drop(eta.H0 %*% t(Z));
  piX = expit.f(covZ)
  Lam0.f <- stepfun(tk_vec, c(0,cumsum(lam_vec)))
  
  if(is.null(beta.H0)){
    covX=0.0
  }else{
    covX=drop(beta.H0 %*% t(X));
  }
  
  Hf=Lam0.f(t)*exp(covX)
  survf <- exp(-Hf)
  res0 <- piX*survf
  res1 <- del + (1-del)*(res0/(1-piX+res0))
  res <- z*res1+(1-z)*(1-res1)
  
  return(res)
}

Smu_vec_coxph_withV.f <- function(del, t, X, Z, beta.H0, eta.H0, lam_vec, tk_vec){
  
  covZ = drop(eta.H0 %*% t(Z));
  piX = expit.f(covZ)
  Lam0.f <- stepfun(tk_vec, c(0,cumsum(lam_vec)))
  
  if(is.null(beta.H0)){
    covX=0.0
  }else{
    covX=drop(beta.H0 %*% t(X));
  }
  
  Hf=Lam0.f(t)*exp(covX)
  wgt <- wgt2_p_coxph_withV.f(1, del, t, X, Z, beta.H0, eta.H0, lam_vec, tk_vec)
  res2 <- del - wgt*Hf
  res1 <- wgt - piX
  #res <- list(Smu1=res1, Smu2=res2)
  
  res3 <- ifelse(del==0, res2, Hf*(1-piX)) #(1-exp(-Hf))*(2*del-1)*(1-piX)
  
  #res <- list(Smu1=(res1-mean(res1))/sd(res1), Smu2=(res2-mean(res2))/sd(res2), Smu3=res3)
  #res <- list(Smu1=res1/sd(res1), Smu2=res2/sd(res2), Smu3=res3)
  res <- list(Smu1=res1, Smu2=res2, Smu3=res3)
  
  return(res)
  
}


designPhII_resid_coxph_withV <- function(phaseI_strat, n2samp, design="RSD-Smu1", w=0.5,
                                         covXname.H0, covZname.H0, 
                                         beta.H0, eta.H0, lam_vec, tk_vec){
  if(is.null(beta)){
    X.H0=NULL
  }else{
    X.H0=phaseI_strat$dt_ext[,covXname.H0,drop=F]
  }
  
  Z.H0=as.matrix(cbind(rep(1, nrow(phaseI_strat$dt_ext)), phaseI_strat$dt_ext[,covZname.H0,drop=F]))
  
  qzi <- Smu_vec_coxph_withV.f(phaseI_strat$dt_ext$del, phaseI_strat$dt_ext$time, X.H0, Z.H0, beta.H0, eta.H0, lam_vec, tk_vec)
  
  k=round(n2samp/2*0.15)
  
  
  if(design %in% c("RSD-Smu1")){
    res0 <- rkd_ext_opt.f(k, n2samp, qzi[[1]], phaseI_strat) 
  }else if(design %in% c("RSD-Smu2")){
    res0 <- rkd_ext_opt.f(k, n2samp, qzi[[2]], phaseI_strat)
  }else if(design %in% paste0("BRSD-eff-seq",w)){
    mat=as.matrix(cbind(qzi[[1]], qzi[[2]]))

    covS=cov(mat)
    Seff1=mat[,1]-covS[2,1]/covS[2,2]*mat[,2]
    Seff2=mat[,2]-covS[1,2]/covS[1,1]*mat[,1]
    
    if(w==0 | w==1){
      qzi_use=w*Seff1+(1-w)*Seff2
      res0 <- rkd_ext_opt.f(k, n2samp, qzi_use, phaseI_strat)
    }else{
      
      n2a=round(n2samp*w)
      ka=round(n2a/2*0.15)
      res01 <- rkd_ext_opt.f(ka, n2a, Seff1, phaseI_strat)
      
      dt_a_use <- phaseI_strat$dt_ext #dt_c
      dt_a_use$R<-1
      dt_a_use[!(dt_a_use$id %in% res01$s_id), "R"] <- 0
      dt_b=subset(dt_a_use, R==0)
      
      n2b=n2samp-n2a
      kb=round(n2b/2*0.15)
      phaseI_strat.b=phaseI_strat
      phaseI_strat.b$dt_ext=dt_b
      phaseI_strat.b$strata.n=phaseI_strat$strata.n-res01$s_m
      res0 <- rkd_ext_opt_2st.f(kb, n2b, Seff2, phaseI_strat, res01$s_id)
      
    }
    
  }
  
  
  res = list(sel=res0, qzi=qzi)
  
  return(res)
}