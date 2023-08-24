
expit.f<-function(x){exp(x)/(1+exp(x))}

wgt2_p_withV.f <- function(z, del, t, X0, Z0, eta2, beta2, alp, brks){
  
  muZ=drop(eta2 %*% t(Z0))
  piX = expit.f(muZ)
  
  if(is.null(X0)){
    muT=rep(0, dim(Z0)[1])
  }else{
    muT=drop(beta2 %*% t(X0))
  }
  survf <- 1-vppwc_H0(t, muT, exp(alp), brks)
  res0 <- piX*survf
  res1 <- del + (1-del)*(res0/(1-piX+res0))
  res <- z*res1+(1-z)*(1-res1)#ifelse(z==1,res1, 1-res1)
  
  return(res)
}

Smu_vec_pwc_withV.f <- function(del, t, X0, Z0, para.H0, brks){
  
  len <- length(brks)+1
  
  if(is.null(X0)){
    len.be2=0
  }else{
    len.be2=(ncol(X0))
  }
  len.et2=ncol(Z0)
  
  alp=(para.H0[1:len]); eta2=para.H0[len+(1:len.et2)]; 
  if(!is.null(X0)){
    beta2=para.H0[len+len.et2+(1:len.be2)]
    muT=drop(beta2 %*% t(X0))
  }else{
      muT=rep(0, dim(Z0)[1])
    }
  muZ=drop(eta2 %*% t(Z0))
  piX = expit.f(muZ)
  
  Hf <- vHpwc_H0(t, muT,exp(alp), brks)
  wgt <- wgt2_p_withV.f(1, del, t, X0, Z0, eta2, beta2, alp,brks)
  res1 <- del - wgt*Hf
  res2 <- wgt - piX
  
  #res3 <- ifelse(del==0, res2, Hf*(1-piX)) #(1-exp(-Hf))*(2*del-1)*(1-piX)
  
  res <- list(Smu1=res1, Smu2=res2)
  
  return(res)
  
}


designPhII_resid_pwc_withV <- function(phaseI_strat, n2samp, design="RSD-Smu1", w=0.5,
                                         covXname.H0, covZname.H0, para.H0,brks){
  
  if(is.null(covXname.H0)){
    X.H0=NULL
  }else{
    X.H0=phaseI_strat$dt_ext[,covXname.H0,drop=F]
  }
  
  Z.H0=as.matrix(cbind(rep(1, nrow(phaseI_strat$dt_ext)), phaseI_strat$dt_ext[,covZname.H0,drop=F]))
  
  qzi <- Smu_vec_pwc_withV.f(phaseI_strat$dt_ext$del, phaseI_strat$dt_ext$time, X.H0, Z.H0, para.H0, brks)
  
  #print(k)
  
  #k=min(round(n2samp/2), sum(phaseI_strat$dt_ext$del))
  k=round(n2samp/2*0.15)
  
  
  if(design %in% c("RSD-Smu1")){
    #k = min((round(n2samp/2)), sum(phaseI_strat$dt_ext$del))
    res0 <- rkd_ext_opt.f(k, n2samp, qzi[[1]], phaseI_strat)
  }else if(design %in% c("RSD-Smu2")){
    res0 <- rkd_ext_opt.f(k, n2samp, qzi[[2]], phaseI_strat)
  }else if(design %in% c("optn-Smu3")){
    res0 <- rkd_ext.f(k, n2samp, qzi[[3]], phaseI_strat)
  # }else if(design %in% c("BRSD")){
  #   qzi_use=qzi[[1]]-qzi[[2]]
  #   #print(head(qzi_use))
  #   res0 <- rkd_ext.f(round(n2samp/2), n2samp, -qzi_use, phaseI_strat) #rkd_max.f(n2samp, qzi_use, phaseI_strat) #rkd_ext.f(k, n2samp, qzi_use, phaseI_strat)
  #   
  # }else if(design %in% c("BRSD2")){
  #   qzi_use=qzi[[1]]+qzi[[2]]
  #   #print(head(qzi_use))
  #   res0 <- rkd_ext.f(round(n2samp/2), n2samp, -qzi_use, phaseI_strat) #rkd_max.f(n2samp, qzi_use, phaseI_strat) #rkd_ext.f(k, n2samp, qzi_use, phaseI_strat)
     
  }else if(design %in% paste0("BRSD",w)){
    #sign_use=sign(qzi[[1]]*qzi[[2]])
    #qzi_use=qzi[[1]]*w-sign_use*qzi[[2]]*(1-w)
    mat=as.matrix(cbind(qzi[[1]], qzi[[2]]))
    qzi0=mat %*% (ginv(cov(mat)))
    qzi_use=qzi0[,1]*w-qzi0[,2]*(1-w)
    #print(head(qzi_use))
    res0 <- rkd_ext_opt.f(k, n2samp, qzi_use, phaseI_strat) #rkd_max.f(n2samp, qzi_use, phaseI_strat) #rkd_ext.f(k, n2samp, qzi_use, phaseI_strat)
    
  }else if(design %in% "BRSD-dc"){
    
    X=cbind(qzi[[1]], qzi[[2]])
    Zm=matrix(1,nrow(X), nrow(X))
    D=X-Zm%*%X/nrow(X)
    T=D%*%sqrtm(ginv(t(D)%*%D))
    qzi_use=T[,1]*T[,2]
    res0 <- rkd_ext_opt.f(k, n2samp, qzi_use, phaseI_strat)
  }else if(design %in% paste0("BRSD-mtp")){
    
    #qzi_use=qzi[[1]]*qzi[[2]]#(qzi[[1]]-mean(qzi[[1]]))*(qzi[[2]]-mean(qzi[[2]]))
    #qzi0=cbind(qzi[[1]], qzi[[2]]) %*% ginv(cov(cbind(qzi[[1]], qzi[[2]])))
    
    mat=as.matrix(cbind(qzi[[1]], qzi[[2]]))
    covS=cov(mat)
    Seff1=mat[,1]-covS[1,2]/covS[2,2]*mat[,2]
    Seff2=mat[,2]-covS[2,1]/covS[1,1]*mat[,1]
    qzi_use=Seff1*Seff2
    res0 <- rkd_ext_opt.f(k, n2samp, qzi_use, phaseI_strat)
    
  }else if(design %in% paste0("BRSD-add",w)){
    #sign_use=sign(qzi[[1]]*qzi[[2]])
    #qzi_use=qzi[[1]]*w-sign_use*qzi[[2]]*(1-w)
    
    mat=as.matrix(cbind(qzi[[1]], qzi[[2]]))
    qzi0= mat %*% (ginv(cov(mat)))
    qzi_use=qzi0[,1]*w+qzi0[,2]*(1-w)
    #print(head(qzi_use))
    res0 <- rkd_ext_opt.f(k, n2samp, qzi_use, phaseI_strat) #rkd_max.f(n2samp, qzi_use, phaseI_strat) #rkd_ext.f(k, n2samp, qzi_use, phaseI_strat)
    
  }else if(design %in% paste0("BRSD-quad",w)){
    mat=as.matrix(cbind(qzi[[1]], qzi[[2]]))
    
    #invB=ginv(cov(mat))
    #qzi0=apply(mat, 1, function(x){
    #  x%*% invB %*% x
    #})
    
    qzi0= mat %*% chol(ginv(cov(mat)))
    qzi_use=w*qzi0[,1]^2+(1-w)*qzi0[,2]^2
    res0 <- rkd_max.f(n2samp, qzi_use, phaseI_strat)
  }else if(design %in% paste0("BRSD-eff",w)){
    mat=as.matrix(cbind(qzi[[1]], qzi[[2]]))
    
    #qzi0= mat %*% (ginv(cov(mat)))
    #qzi_use=mat[,1]- qzi0[,2] #w*qzi0[,1]^2+(1-w)*qzi0[,2]^2
    covS=cov(mat)
    Seff1=mat[,1]-covS[1,2]/covS[2,2]*mat[,2]
    Seff2=mat[,2]-covS[2,1]/covS[1,1]*mat[,1]
    qzi_use=Seff1-Seff2#w*Seff1+(1-w)*Seff2#Seff2^2
    res0 <- rkd_ext_opt.f(k, n2samp, qzi_use, phaseI_strat)
    #res0 <- rkd_max.f(n2samp, qzi_use, phaseI_strat)
  }else if(design %in% paste0("BRSD-eff-seq",w)){
    mat=as.matrix(cbind(qzi[[1]], qzi[[2]]))
    
    #qzi0= mat %*% (ginv(cov(mat)))
    #qzi_use=mat[,1]- qzi0[,2] #w*qzi0[,1]^2+(1-w)*qzi0[,2]^2
    covS=cov(mat)
    Seff1=mat[,1]-covS[2,1]/covS[2,2]*mat[,2]
    Seff2=mat[,2]-covS[1,2]/covS[1,1]*mat[,1]
    #qzi_use=Seff1-Seff2#w*Seff1+(1-w)*Seff2#Seff2^2
    
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
      #res02_1=rkd_ext_opt.f(kb, n2b, Seff2[dt_a_use$R==0], phaseI_strat.b)
      res0 <- rkd_ext_opt_2st.f(kb, n2b, Seff2, phaseI_strat, res01$s_id)
      
      #s_m=res01$s_m+res02$s_m
      #s_prob=s_m/phaseI_strat$strata.n
      
      #res0=res02#list(s_id=c(res01$s_id, res02$s_id), s_m=s_m, s_prob=s_prob)
    }
    
    #res0 <- rkd_max.f(n2samp, qzi_use, phaseI_strat)
  }
  
  
  res = list(sel=res0, qzi=qzi)
  
  return(res)
}


designPhII_resid_pwc_noV <- function(phaseI_strat, n2samp, design="RSD-Smu1", w=0.5,
                                         covXname.H0, covZname.H0, 
                                         para.H0,brks){
  
  if(is.null(covXname.H0)){
    X.H0=NULL
  }else{
    X.H0=phaseI_strat$dt_ext[,covXname.H0,drop=F]
  }
  
  Z.H0=as.matrix(cbind(rep(1, nrow(phaseI_strat$dt_ext)), phaseI_strat$dt_ext[,covZname.H0,drop=F]))
  
  qzi <- Smu_vec_pwc_withV.f(phaseI_strat$dt_ext$del, phaseI_strat$dt_ext$time, X.H0, Z.H0, para.H0, brks)
  
  #print(k)
  
  #k=min(round(n2samp/2), sum(phaseI_strat$dt_ext$del))
  k=round(n2samp/2*0.15)
  
  
  if(design %in% c("RSD-Smu1")){
    #k = min((round(n2samp/2)), sum(phaseI_strat$dt_ext$del))
    #res0 <- rkd_ext_opt_ignoreV.f(k, n2samp, qzi[[1]], phaseI_strat)
    res0 <- rkd_ext.f(n2samp/2, n2samp,  qzi[[1]], phaseI_strat)
  }else if(design %in% c("RSD-Smu2")){
    res0 <- rkd_ext.f(n2samp/2, n2samp, qzi[[2]], phaseI_strat)
  }else if(design %in% c("optn-Smu3")){
    res0 <- rkd_ext_opt_ignoreV.f(k, n2samp, qzi[[3]], phaseI_strat)

  }else if(design %in% paste0("BRSD",w)){
    #sign_use=sign(qzi[[1]]*qzi[[2]])
    #qzi_use=qzi[[1]]*w-sign_use*qzi[[2]]*(1-w)
    mat=as.matrix(cbind(qzi[[1]], qzi[[2]]))
    qzi0=mat %*% (ginv(cov(mat)))
    qzi_use=qzi0[,1]*w-qzi0[,2]*(1-w)
    #print(head(qzi_use))
    res0 <- rkd_ext_opt_ignoreV.f(k, n2samp, qzi_use, phaseI_strat) #rkd_max.f(n2samp, qzi_use, phaseI_strat) #rkd_ext.f(k, n2samp, qzi_use, phaseI_strat)
    
  }else if(design %in% "BRSD-dc"){
    
    X=cbind(qzi[[1]], qzi[[2]])
    Zm=matrix(1,nrow(X), nrow(X))
    D=X-Zm%*%X/nrow(X)
    T=D%*%sqrtm(ginv(t(D)%*%D))
    qzi_use=T[,1]*T[,2]
    res0 <- rkd_ext_opt_ignoreV.f(k, n2samp, qzi_use, phaseI_strat)
  }else if(design %in% paste0("BRSD-mtp")){
    
    #qzi_use=qzi[[1]]*qzi[[2]]#(qzi[[1]]-mean(qzi[[1]]))*(qzi[[2]]-mean(qzi[[2]]))
    #qzi0=cbind(qzi[[1]], qzi[[2]]) %*% ginv(cov(cbind(qzi[[1]], qzi[[2]])))
    
    mat=as.matrix(cbind(qzi[[1]], qzi[[2]]))
    covS=cov(mat)
    Seff1=mat[,1]-covS[1,2]/covS[2,2]*mat[,2]
    Seff2=mat[,2]-covS[2,1]/covS[1,1]*mat[,1]
    qzi_use=Seff1*Seff2
    res0 <- rkd_ext_opt_ignoreV.f(k, n2samp, qzi_use, phaseI_strat)
    
  }else if(design %in% paste0("BRSD-add",w)){
    #sign_use=sign(qzi[[1]]*qzi[[2]])
    #qzi_use=qzi[[1]]*w-sign_use*qzi[[2]]*(1-w)
    
    mat=as.matrix(cbind(qzi[[1]], qzi[[2]]))
    qzi0= mat %*% (ginv(cov(mat)))
    qzi_use=qzi0[,1]*w+qzi0[,2]*(1-w)
    #print(head(qzi_use))
    res0 <- rkd_ext_opt_ignoreV.f(k, n2samp, qzi_use, phaseI_strat) #rkd_max.f(n2samp, qzi_use, phaseI_strat) #rkd_ext.f(k, n2samp, qzi_use, phaseI_strat)
    
  }else if(design %in% paste0("BRSD-quad",w)){
    mat=as.matrix(cbind(qzi[[1]], qzi[[2]]))
    
    #invB=ginv(cov(mat))
    #qzi0=apply(mat, 1, function(x){
    #  x%*% invB %*% x
    #})
    
    qzi0= mat %*% chol(ginv(cov(mat)))
    qzi_use=w*qzi0[,1]^2+(1-w)*qzi0[,2]^2
    res0 <- rkd_max.f(n2samp, qzi_use, phaseI_strat)
  }else if(design %in% paste0("BRSD-eff",w)){
    mat=as.matrix(cbind(qzi[[1]], qzi[[2]]))
    
    #qzi0= mat %*% (ginv(cov(mat)))
    #qzi_use=mat[,1]- qzi0[,2] #w*qzi0[,1]^2+(1-w)*qzi0[,2]^2
    covS=cov(mat)
    Seff1=mat[,1]-covS[1,2]/covS[2,2]*mat[,2]
    Seff2=mat[,2]-covS[2,1]/covS[1,1]*mat[,1]
    qzi_use=Seff1-Seff2#w*Seff1+(1-w)*Seff2#Seff2^2
    res0 <- rkd_ext.f(k, n2samp, qzi_use, phaseI_strat)
    #res0 <- rkd_max.f(n2samp, qzi_use, phaseI_strat)
  }else if(design %in% paste0("BRSD-eff-seq",w)){
    mat=as.matrix(cbind(qzi[[1]], qzi[[2]]))
    
    #qzi0= mat %*% (ginv(cov(mat)))
    #qzi_use=mat[,1]- qzi0[,2] #w*qzi0[,1]^2+(1-w)*qzi0[,2]^2
    covS=cov(mat)
    Seff1=mat[,1]-covS[2,1]/covS[2,2]*mat[,2]
    Seff2=mat[,2]-covS[1,2]/covS[1,1]*mat[,1]
    #qzi_use=Seff1-Seff2#w*Seff1+(1-w)*Seff2#Seff2^2
    
    if(w==0 | w==1){
      qzi_use=w*Seff1+(1-w)*Seff2
      #res0 <- rkd_ext_opt_noV.f(k, n2samp, qzi_use, phaseI_strat)
      res0 <- rkd_ext.f(n2samp/2, n2samp, qzi_use, phaseI_strat)
    }else{
      
      n2a=round(n2samp*w)
      ka=round(n2a/2*0.15)
      #res01 <- rkd_ext_opt_ignoreV.f(ka, n2a, Seff1, phaseI_strat)
      res01 <- rkd_ext.f(n2a/2, n2a, Seff1, phaseI_strat)
      
      dt_a_use <- phaseI_strat$dt_ext #dt_c
      dt_a_use$R<-1
      dt_a_use[!(dt_a_use$id %in% res01$s_id), "R"] <- 0
      dt_b=subset(dt_a_use, R==0)
      
      n2b=n2samp-n2a
      kb=round(n2b/2*0.15)
      phaseI_strat.b=phaseI_strat
      phaseI_strat.b$dt_ext=dt_b
      phaseI_strat.b$strata.n=phaseI_strat$strata.n-res01$s_m
      #res02_1=rkd_ext_opt.f(kb, n2b, Seff2[dt_a_use$R==0], phaseI_strat.b)
      #res0 <- rkd_ext_opt_ignoreV_2st.f(kb, n2b, Seff2, phaseI_strat, res01$s_id)
      #res02 <- rkd_ext.f(n2b/2, n2b, Seff2[!(dt_a_use$id %in% res01$s_id)], phaseI_strat.b)
      res02<-rkd_ext_opt_ignoreV_2st.f(kb, n2b, Seff2, phaseI_strat, res01$s_id)
      #res02_1<-rkd_ext_opt_2st.f(kb, n2b, Seff2, phaseI_strat, res01$s_id)
      
      #s_m=res01$s_m+res02$s_m
      s_m=res02$s_m
      s_prob=s_m/phaseI_strat$strata.n
      
      #res0=list(s_id=c(res01$s_id, res02$s_id), s_m=s_m, s_prob=s_prob)
      res0=list(s_id=c(res02$s_id), s_m=s_m, s_prob=s_prob)
    }
    
    #res0 <- rkd_max.f(n2samp, qzi_use, phaseI_strat)
  }
  
  
  res = list(sel=res0, qzi=qzi)
  
  return(res)
}




