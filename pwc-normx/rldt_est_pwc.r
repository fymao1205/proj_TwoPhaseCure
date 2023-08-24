

# max. obs loglik

rldt_est_pwc_H0.f <- function(del, time, X0, Z0, brks){
  
  len = length(brks)+1
  
  if(is.null(X0)){
    len.be2=0
  }else{
    len.be2=(ncol(X0))
  }
  
  len.et2=ncol(Z0)

  
  obj.f <- function(para){
    
    alp = exp(para[1:len]) 
    
    eta2 = para[len+(1:len.et2)]
    muZ=drop(eta2%*%t(Z0))
    
    if(!is.null(X0)){
      beta2 = para[len+len.et2+(1:len.be2)]
      muT=drop(beta2%*%t(X0))
    }else{
      muT=rep(0, length(time))
    }
    
    res0 <- logLn_obs_pwc_norm_H0( del, time, muZ, muT, alp, brks)
    
    res <- -sum(res0)
    
    print(para)
    print(res)
    
    return(res)
    
  }
  
  est <- nlm(f=obj.f, p=rep(0.01, len+len.et2+len.be2), gradtol = 1e-03)
  return(est)
}


rldt_est_pwc.f <- function(ini, scen, del, time, x, X0, Z0, brks){
  
  if(is.null(ini)){
    if(scen=="sc"){
      ini=rep(0.01,len+4+len.et2+len.be2)
    }else{
      ini=rep(0.01,len+3+len.et2+len.be2)
    }
    
  }
  
  r <- ifelse(is.na(x), 0, 1)
  len = length(brks)+1
  if(is.null(X0)){
    len.be2=0
  }else{
    len.be2=(ncol(X0))
  }
  
  len.et2=ncol(Z0)
  
  if(scen=="sc"){
    obj.f <- function(para){
      
      alp = exp(para[1:len]) 
      eta1 = para[len+1]
      beta1 = para[len+2]
      
      eta2 = para[len+2+(1:len.et2)]
      muX1=para[len+2+len.et2+len.be2+1]
      sigma2 = exp(para[len+2+len.et2+len.be2+2])
      
      muZ=drop(eta2%*%t(Z0)); 
      
      if(!is.null(X0)){
        beta2 = para[len+2+len.et2+(1:len.be2)]
        muT=drop(beta2%*%t(X0))
      }else{
        muT=rep(0, length(time))
      }
      
      
      res0 <- logLn_obs_pwc_norm(r, del, time, x, eta1, beta1, muZ, muT, muX1, sigma2, alp, brks)
      
      res <- -sum(res0)
      
      print(para)
      print(res)
      
      return(res)
      
    }
  }else if(scen=="sa"){
    
    obj.f <- function(para){
      
      alp = exp(para[1:len]) 
      eta1 = para[len+1]
      beta1 = 0
      
      eta2 = para[len+1+(1:len.et2)]
      #beta2 = para[len+1+len.et2+(1:len.be2)]
      muX1=para[len+1+len.et2+len.be2+1]
      sigma2 = exp(para[len+1+len.et2+len.be2+2])
      
      muZ=drop(eta2%*%t(Z0));
      
      if(!is.null(X0)){
        beta2 = para[len+1+len.et2+(1:len.be2)]
        muT=drop(beta2%*%t(X0))
      }else{
        muT=rep(0, length(time))
      }
    
      res0 <- logLn_obs_pwc_norm(r, del, time, x, eta1, beta1, muZ, muT, muX1, sigma2, alp, brks)
      
      
      res <- -sum(res0)
      
      print(para)
      print(res)
      
      return(res)
      
    }
  }else if(scen=="sb"){
    
    obj.f <- function(para){
      
      alp = exp(para[1:len]) 
      eta1 = 0
      beta1 = para[len+1]
      
      eta2 = para[len+1+(1:len.et2)]
      #beta2 = para[len+1+len.et2+(1:len.be2)]
      muX1=para[len+1+len.et2+len.be2+1]
      sigma2 = exp(para[len+1+len.et2+len.be2+2])
      
      muZ=drop(eta2%*%t(Z0)); #muT=drop(beta2%*%t(X0))
      
      
      if(!is.null(X0)){
        beta2 = para[len+1+len.et2+(1:len.be2)]
        muT=drop(beta2%*%t(X0))
      }else{
        muT=rep(0, length(time))
      }
      
      
      res0 <- logLn_obs_pwc_norm(r, del, time, x, eta1, beta1, muZ, muT, muX1, sigma2, alp, brks)
      #}
      
      
      res <- -sum(res0)
      #  res<- -mean(res0)
      
      print(para)
      print(res)
      
      return(res)
      
    }
  }
  
  
  est <- optim(par=ini, obj.f, method=("L-BFGS-B"), hessian=TRUE)
  #est <- nlm(f=obj.f, p=ini, hessian=TRUE, gradtol = 1e-05)
  return(est)
}


logL_pwc_norm.f <- function(scen, est, del, time, x, X0, Z0, brks){
  
  r <- ifelse(is.na(x), 0, 1)
  len <- length(brks)+1
  
  if(is.null(X0)){
    len.be2=0
  }else{
    len.be2=(ncol(X0))
  }
  len.et2=ncol(Z0)
  
  alp = exp(est[1:len]) 
  
  if(scen=="sc"){
    eta1 = est[len+1]
    beta1 = est[len+2]
    
    eta2 = est[len+2+(1:len.et2)]
    if(!is.null(X0)){beta2 =est[len+2+len.et2+(1:len.be2)]}
    muX1=est[len+2+len.et2+len.be2+1]
    sigma2 = exp(est[len+2+len.et2+len.be2+2])
    
  }else if(scen=="sa"){
    eta1 = est[len+1]
    beta1 = 0#est[len+2]
    
    eta2 = est[len+1+(1:len.et2)]
    if(!is.null(X0)){beta2 =est[len+1+len.et2+(1:len.be2)]}
    muX1=est[len+1+len.et2+len.be2+1]
    sigma2 = exp(est[len+1+len.et2+len.be2+2])
  }else if(scen=="sb"){
    eta1 = 0#est[len+1]
    beta1 = est[len+1]
    
    eta2 = est[len+1+(1:len.et2)]
    if(!is.null(X0)){beta2 =est[len+1+len.et2+(1:len.be2)]}
    muX1=est[len+1+len.et2+len.be2+1]
    sigma2 = exp(est[len+1+len.et2+len.be2+2])
  }
  
  
  muZ=drop(eta2%*%t(Z0)); 
  
  if(is.null(X0)){
    muT=rep(0, length(time))
  }else{
    muT=drop(beta2%*%t(X0))
  }
  
  res0 <- logLn_obs_pwc_norm(r, del, time, x, eta1, beta1, muZ, muT, muX1, sigma2, alp, brks)
  
  return(res0)
}


score_pwc_norm.f <- function(grad=1e-06, scen, para, del, time, x, X0, Z0, brks){
  
  logL0 <- logL_pwc_norm.f(scen, para, del, time, x, X0, Z0, brks)
  sc.list <- lapply(1:length(para), function(j){
    para.p <- para
    para.p[j] <- para[j]+grad
    
    logL.p <- logL_pwc_norm.f(scen, para.p, del, time, x, X0, Z0, brks)
    (logL0-logL.p)/grad 
  })
  
  res <- do.call("cbind", sc.list)
  return(res)
}

hess_pwc_norm.f <- function(grad=1e-06, scen, para, del, time, x, X0, Z0, brks){
  
  s0 <- score_pwc_norm.f(grad, scen, para, del, time, x, X0, Z0, brks)
  h.list <- lapply(1:length(para), function(j){
    para.p <- para
    para.p[j] <- para[j]+grad
    
    s.p <- score_pwc_norm.f(grad, scen, para.p, del, time, x, X0, Z0, brks)
    (colSums(s.p)-colSums(s0))/grad 
  })
  
  res <- do.call("cbind", h.list)
  return(res)
  
}

ase_pwc_norm.f <- function(grad=1e-06, scen, est, del, time, x, X0, Z0,brks){
  # score vector
  score.mat = score_pwc_norm.f(grad, scen, est$par, del, time, x, X0, Z0, brks)
  #score.mat = score_trueh.f(grad, est$par, del, time, v.mat, pos.cov.T, pos.cov.Z, pos.cov.X)
  
  #hess.mat = hess_pwc_norm.f(grad, scen, est$par, del, time, x, X0, Z0, brks)
  hess.mat = est$hess 
  #hessian(lik.f, est, x1=dt.cov[,1], v=dt.cov[,-1, drop=FALSE], del1, del2, b1, b2, cov_dist) # hessian matrix by numeircal approximation
  B = t(score.mat) %*% score.mat
  A = (hess.mat)
  #A = (hesmat)
  avar = ginv(A)%*% B%*% ginv(A)
  ase = sqrt(diag(avar))
  res = list(ase=ase, A=A, B=B, score=colSums(score.mat))
  return(res)
}


BootVar_pwc <- function(bootnum, ini, scen, dat, X0, Z0, brks){
  
  WRAPPERFUNC2 = function(data, indices){
    dat_b=data[indices,,drop=FALSE]
    fitTEMP = rldt_est_pwc.f(ini=ini, scen=scen, dat_b$del, dat_b$time, dat_b$x, 
                             X0[indices,,drop=FALSE], Z0[indices,,drop=FALSE], brks=brks)
    return(fitTEMP$par)
  }	
  
  WRAPPERFUNC = tryCatch(WRAPPERFUNC2,error = function(cond){return(rep(NA,P))},warning = function(cond){return(rep(NA,P))})
  
  BOOT = boot::boot(data = dat, statistic = WRAPPERFUNC, R = bootnum, parallel = "multicore", ncpus=8)
  ### Guarding against non-converging bootstrap samples
  #BOOT$t = ifelse(abs(BOOT$t)>5, matrix(NA,ncol = ncol(BOOT$t), nrow = nrow(BOOT$t)), BOOT$t)
  
  SD = apply(BOOT$t, 2, sd)
  
  return(list(Boot=BOOT, sd=SD))
}

