P_Y_g_X.f <- function(del, t, x1, x2, eta, alp, beta){
  
  covZ = eta[1] + eta[2]*x1+ eta[3]*x2
  piX = expit.f(covZ)
  hf <- alp[2]*alp[1]*(alp[1]*t)^(alp[2]-1)*exp(beta[1]*x1+beta[2]*x2) #hpwc(t, brks, exp(alp), 0)
  Hf <- (alp[1]*t)^(alp[2])*exp(beta[1]*x1+beta[2]*x2)
  survf <- exp(-Hf) #vppwc(t, brks, exp(alp), 1, 0)
  densf <- hf*survf
  res0 <- 1-piX*(1-survf)
  res <- del*densf*piX+(1-del)*res0  #ifelse(del==1, densf, res0)
  
  return(res)
  
}

# gaussian legendre 
# prepare gaussleg.f(): reference of Jooyoung's code
gaussleg.f <- function(n, x1, x2) {
  EPS <- 3e-14
  
  m <- (n + 1)/2
  xm <- 0.5 * (x2 + x1)
  xl <- 0.5 * (x2 - x1)
  
  x <- rep(0, n)
  w <- rep(0, n)
  for (i in 1:m) {
    z <- cos(pi * (i - 0.25)/(n + 0.5))
    
    tol <- 9999
    while (tol > EPS) {
      p1 <- 1
      p2 <- 0
      for (j in 1:n) {
        p3 <- p2
        p2 <- p1
        p1 <- ((2 * j - 1) * z * p2 - (j - 1) * p3)/j
      }
      
      pp <- n * (z * p1 - p2)/(z * z - 1)
      z1 <- z
      z <- z1 - p1/pp
      
      tol <- abs(z - z1)
      if (tol <= EPS) {
        break
      }
    }
    
    x[i] = xm - xl * z
    x[n + 1 - i] = xm + xl * z
    w[i] = (2 * xl)/((1 - z * z) * pp * pp)
    w[n + 1 - i] = w[i]
  }
  
  return(as.matrix(data.frame(location = x, weight = w)))
}

logL_wei.f <- function(est, r, del, time, x, v){
  
  #r <- ifelse(is.na(x), 0, 1)
  
  alp = c(exp(est[1]), est[2]) #exp(est[1]) #
  beta = est[3:4]
  eta = est[(5:7)]
  gam = est[(8:9)]
  
  #res0 <- logLn_obs_wei(r, del, time, x, v, c(eta), alp, gam, c(beta))
  px1=expit.f(gam[1]+gam[2]*v); px0=1-px1
  
  py1 = P_Y_g_X.f(del, time, 1, v, eta, alp, beta)*px1
  py0 = P_Y_g_X.f(del, time, 0, v, eta, alp, beta)*px0
  
  if(r==1){
    pyxr=x*py1+(1-x)*py0 
  }else{
    pyxr=py1+py0
  }
  #pyxr=r*pyx+(1-r)*(py1+py0)#ifelse(r==1,pyx,py1+py0)
  res0=log(pyxr)
  
  return(res0)
}

score_wei_ele.f <- function(grad=1e-06, pos, para, r, del, time, x, v){
  
  para.p = para
  para.p[pos] = para[pos]+grad
  
  logL0 <- logL_wei.f(para, r, del, time, x, v)
  logL.p <- logL_wei.f(para.p, r, del, time, x, v)
  
  res <- (logL.p-logL0)/grad 
  
  return(res)
}

SS_wei_ele.f <- function(grad=1e-06, pos, para, r, del, time, x, v){
  
  s0 <- score_wei_ele.f(grad, pos, para, r, del, time, x, v)
  
  s0^2
}

S1S2_wei_ele.f <- function(grad=1e-06, pos1, pos2, para, r, del, time, x, v){
  
  s1 <- score_wei_ele.f(grad, pos1, para, r, del, time, x, v)
  s2 <- score_wei_ele.f(grad, pos2, para, r, del, time, x, v)
  
  s1*s2
}

calS1S2_ele_approx2.f <- function(grad, pos1, pos2, cuts, sp11, sp12, sp13, sp01, sp02, sp03, para, sc0, sh0, pv){
  
  alp = c(exp(para[1]), para[2]) 
  beta = para[3:4]
  eta = para[(5:7)]
  gam = para[(8:9)]
  
  if(gam[2]==0){
    p1v1=p1v0 = expit.f(gam[1]); 
  }else{
    p1v1 = expit.f(gam[1]+gam[2]); p1v0 = expit.f(gam[1]); 
  }
  p0v1=1-p1v1;p0v0=1-p1v0
  #piZ_1v1 = expit.f(eta[1]+eta[2]+eta[3]); piZ_0v1 = expit.f(eta[1]+eta[3])
  #piZ_1v0 = expit.f(eta[1]+eta[2]); piZ_0v0 = expit.f(eta[1])
  
  sp.f<-function(del, t){
    
    #sp_vec=del*c(sp11, sp12, sp13)+(1-del)*c(sp01, sp02, sp03)
    
    ind1=ifelse(t<=cuts[1],1,0)
    ind2=ifelse(t>cuts[1]&t<=cuts[2],1,0)
    ind3=ifelse(t>cuts[2],1,0)
    
    res=ind1*(del*sp11+(1-del)*sp01)+ind2*(del*sp12+(1-del)*sp02)+ind3*(del*sp13+(1-del)*sp03) 
    return(res)
  }
  
  int1.f <- function(t){
    
    sp=sp.f(1, t)
    
    py1x1v1 = P_Y_g_X.f(1, t, 1, 1, eta, alp, beta) 
    py1x1v0 = P_Y_g_X.f(1, t, 1, 0, eta, alp, beta) 
    py1x0v1 = P_Y_g_X.f(1, t, 0, 1, eta, alp, beta) 
    py1x0v0 = P_Y_g_X.f(1, t, 0, 0, eta, alp, beta) 
    
    sr1x1v1=S1S2_wei_ele.f(grad, pos1, pos2, para, 1, 1, t, 1, 1) 
    sr1x0v1=S1S2_wei_ele.f(grad, pos1, pos2, para, 1, 1, t, 0, 1)
    sr1x1v0=S1S2_wei_ele.f(grad, pos1, pos2, para, 1, 1, t, 1, 0) 
    sr1x0v0=S1S2_wei_ele.f(grad, pos1, pos2, para, 1, 1, t, 0, 0)
    sr0v1=S1S2_wei_ele.f(grad, pos1, pos2, para, 0, 1, t, NA, 1) 
    sr0v0=S1S2_wei_ele.f(grad, pos1, pos2, para, 0, 1, t, NA, 0) 
    
    #val1v1=py1x1v1*piZ_1v1*p1v1; val0v1=py1x0v1*piZ_0v1*p0v1
    #val1v0=py1x1v0*piZ_1v0*p1v0; val0v0=py1x0v0*piZ_0v0*p0v0
    val1v1=py1x1v1*p1v1; val0v1=py1x0v1*p0v1
    val1v0=py1x1v0*p1v0; val0v0=py1x0v0*p0v0
    
    res0=(((sr1x1v1*val1v1+sr1x0v1*val0v1)*pv+(sr1x1v0*val1v0+sr1x0v0*val0v0)*(1-pv))*sp
          + ((sr0v1*(val1v1+val0v1))*pv+sr0v0*(val1v0+val0v0)*(1-pv))*(1-sp))*pweibull(t,sh0,sc0,lower.tail = F)
    
    res0 <- ifelse(is.na(res0), 0, res0)
    res0 <- ifelse(is.finite(res0), res0, 0)
    res0
  }
  
  int0.f <- function(t){
    
    sp=sp.f(0, t)
    
    py1x1v1 = P_Y_g_X.f(0, t, 1, 1, eta, alp, beta) 
    py1x1v0 = P_Y_g_X.f(0, t, 1, 0, eta, alp, beta) 
    py1x0v1 = P_Y_g_X.f(0, t, 0, 1, eta, alp, beta) 
    py1x0v0 = P_Y_g_X.f(0, t, 0, 0, eta, alp, beta) 
    
    sr1x1v1=S1S2_wei_ele.f(grad, pos1, pos2, para, 1, 0, t, 1, 1) 
    sr1x0v1=S1S2_wei_ele.f(grad, pos1, pos2, para, 1, 0, t, 0, 1)
    sr1x1v0=S1S2_wei_ele.f(grad, pos1, pos2, para, 1, 0, t, 1, 0) 
    sr1x0v0=S1S2_wei_ele.f(grad, pos1, pos2, para, 1, 0, t, 0, 0)
    sr0v1=S1S2_wei_ele.f(grad, pos1, pos2, para, 0, 0, t, NA, 1) 
    sr0v0=S1S2_wei_ele.f(grad, pos1, pos2, para, 0, 0, t, NA, 0) 
    
    val1v1=py1x1v1*p1v1; val0v1=py1x0v1*p0v1
    val1v0=py1x1v0*p1v0; val0v0=py1x0v0*p0v0
    
    res0=(((sr1x1v1*val1v1+sr1x0v1*val0v1)*pv+(sr1x1v0*val1v0+sr1x0v0*val0v0)*(1-pv))*sp
          + ((sr0v1*(val1v1+val0v1))*pv+sr0v0*(val1v0+val0v0)*(1-pv))*(1-sp))*dweibull(t,sh0,sc0)
    
    res0 <- ifelse(is.na(res0), 0, res0)
    res0 <- ifelse(is.finite(res0), res0, 0)
    res0
  }
  
  uu=gaussleg.f(20,0,1)
  
  res11=sum(int1.f(uu[,1])*uu[,2])
  res10=sum(int0.f(uu[,1])*uu[,2])
  
  res21=sum(int1.f(1/uu[,1])*(uu[,1])^(-2)*uu[,2])
  res20=sum(int0.f(1/uu[,1])*(uu[,1])^(-2)*uu[,2])
  
  res=res11+res10+res21+res20
  
  return(res)
}

calSS_ele_approx2.f <- function(grad, pos, cuts, sp11, sp12, sp13, sp01, sp02, sp03, para, sc0, sh0, pv){
  
  alp = c(exp(para[1]), para[2]) 
  beta = para[3:4]
  eta = para[(5:7)]
  gam = para[(8:9)]
  
  if(gam[2]==0){
    p1v1=p1v0 = expit.f(gam[1]); 
  }else{
    p1v1 = expit.f(gam[1]+gam[2]); p1v0 = expit.f(gam[1]); 
  }
  p0v1=1-p1v1;p0v0=1-p1v0
  #piZ_1v1 = expit.f(eta[1]+eta[2]+eta[3]); piZ_0v1 = expit.f(eta[1]+eta[3])
  #piZ_1v0 = expit.f(eta[1]+eta[2]); piZ_0v0 = expit.f(eta[1])
  
  sp.f<-function(del, t){
    
    #sp_vec=del*c(sp11, sp12, sp13)+(1-del)*c(sp01, sp02, sp03)
    
    ind1=ifelse(t<=cuts[1],1,0)
    ind2=ifelse(t>cuts[1]&t<=cuts[2],1,0)
    ind3=ifelse(t>cuts[2],1,0)
    
    res=ind1*(del*sp11+(1-del)*sp01)+ind2*(del*sp12+(1-del)*sp02)+ind3*(del*sp13+(1-del)*sp03) 
    
    return(res)
  }
  
  int1.f <- function(t){
    
    sp=sp.f(1,t)
    
    
    py1x1v1 = P_Y_g_X.f(1, t, 1, 1, eta, alp, beta) 
    py1x1v0 = P_Y_g_X.f(1, t, 1, 0, eta, alp, beta) 
    py1x0v1 = P_Y_g_X.f(1, t, 0, 1, eta, alp, beta) 
    py1x0v0 = P_Y_g_X.f(1, t, 0, 0, eta, alp, beta) 
    
    sr1x1v1 =  SS_wei_ele.f(grad, pos, para, 1, 1, t, 1, 1) 
    sr1x0v1 = SS_wei_ele.f(grad, pos, para, 1, 1, t, 0, 1)
    sr1x1v0 =  SS_wei_ele.f(grad, pos, para, 1, 1, t, 1, 0) 
    sr1x0v0 = SS_wei_ele.f(grad, pos, para, 1, 1, t, 0, 0)
    sr0v1 = SS_wei_ele.f(grad, pos, para, 0, 1, t, NA, 1)
    sr0v0 = SS_wei_ele.f(grad, pos, para, 0, 1, t, NA, 0)
    
    #val1v1=py1x1v1*piZ_1v1*p1v1; val0v1=py1x0v1*piZ_0v1*p0v1
    #val1v0=py1x1v0*piZ_1v0*p1v0; val0v0=py1x0v0*piZ_0v0*p0v0
    val1v1=py1x1v1*p1v1; val0v1=py1x0v1*p0v1
    val1v0=py1x1v0*p1v0; val0v0=py1x0v0*p0v0
    
    
    res0=(((sr1x1v1*val1v1+sr1x0v1*val0v1)*pv+(sr1x1v0*val1v0+sr1x0v0*val0v0)*(1-pv))*sp
          + ((sr0v1*(val1v1+val0v1))*pv+sr0v0*(val1v0+val0v0)*(1-pv))*(1-sp))*pweibull(t,sh0,sc0,lower.tail = F)
    
    res0 <- ifelse(is.na(res0), 0, res0)
    res0 <- ifelse(is.finite(res0), res0, 0)
    res0
  }
  
  
  int0.f <- function(t){
    
    sp=sp.f(0,t)
    
    py0x1v1 = P_Y_g_X.f(0, t, 1, 1, eta, alp, beta) 
    py0x1v0 = P_Y_g_X.f(0, t, 1, 0, eta, alp, beta) 
    py0x0v1 = P_Y_g_X.f(0, t, 0, 1, eta, alp, beta) 
    py0x0v0 = P_Y_g_X.f(0, t, 0, 0, eta, alp, beta) 
    
    sr1x1v1 =  SS_wei_ele.f(grad, pos, para, 1, 0, t, 1, 1) 
    sr1x0v1 = SS_wei_ele.f(grad, pos, para, 1, 0, t, 0, 1)
    sr1x1v0 =  SS_wei_ele.f(grad, pos, para, 1, 0, t, 1, 0) 
    sr1x0v0 = SS_wei_ele.f(grad, pos, para, 1, 0, t, 0, 0)
    sr0v1 = SS_wei_ele.f(grad, pos, para, 0, 0, t, NA, 1)
    sr0v0 = SS_wei_ele.f(grad, pos, para, 0, 0, t, NA, 0)
    
    val1v1=py0x1v1*p1v1; val0v1=py0x0v1*p0v1
    val1v0=py0x1v0*p1v0; val0v0=py0x0v0*p0v0
    
    res0=(((sr1x1v1*val1v1+sr1x0v1*val0v1)*pv+(sr1x1v0*val1v0+sr1x0v0*val0v0)*(1-pv))*sp
          + ((sr0v1*(val1v1+val0v1))*pv+sr0v0*(val1v0+val0v0)*(1-pv))*(1-sp))*dweibull(t,sh0,sc0)
    
    res0 <- ifelse(is.na(res0), 0, res0)
    res0 <- ifelse(is.finite(res0), res0, 0)
    res0
  }
  
  
  uu=gaussleg.f(20,0,1)
  
  res11=sum(int1.f(uu[,1])*uu[,2])
  res10=sum(int0.f(uu[,1])*uu[,2])
  
  res21=sum(int1.f(1/uu[,1])*(uu[,1])^(-2)*uu[,2])
  res20=sum(int0.f(1/uu[,1])*(uu[,1])^(-2)*uu[,2])
  
  res=res11+res10+res21+res20
  
  return(res)
}

info_mat.f <- function(pos_vec, grad, cuts, sp11, sp12, sp13, sp01, sp02, sp03, para, sc0, sh0, pv){
  
  
  p = length(pos_vec)
  res <- matrix(NA, p, p)
  
  run0=expand.grid(pos_vec, pos_vec)
  run1=run0[run0[,1]<run0[,2],]
  run2=run0[run0[,1]==run0[,2],]
  
  vec1 = mclapply(1:nrow(run1), function(i){

    x=as.vector(unlist(run1[i,]))
    calS1S2_ele_approx2.f(grad, x[1], x[2], cuts, sp11, sp12, sp13, sp01, sp02, sp03, para, sc0, sh0, pv)
    
  }, mc.cores=3)
  
  vec2 = mclapply(1:nrow(run2), function(i){

    x=as.vector(unlist(run2[i,]))
    calSS_ele_approx2.f(grad, x[1], cuts, sp11, sp12, sp13, sp01, sp02, sp03, para, sc0, sh0, pv)
    
  }, mc.cores=3)
  
  diag(res) <- as.vector(unlist(vec2))
  res[upper.tri(res)] <-as.vector(unlist(vec1))
  tres <- t(res)
  tres[upper.tri(tres)] <- as.vector(unlist(vec1))
  
  return(tres)
}


opt1_ODS_asym_AY.f <- function(scen, grad, n_sample, phaseI_res, ini, para, sc0, sh0, pv){
  
  if(scen=="sc"){
    pos_vec=1:length(para)
    target_pos=6
  }else if(scen=="sa"){
    pos_vec=(1:length(para))[-3]
    target_pos=5
  }
  
  dt_ext = phaseI_res$dt_ext
  ref.group = phaseI_res$ref.strata
  strata.n = phaseI_res$strata.n
  n <- dim(dt_ext)[1]
  dt_ext$id <- 1:n
  
  p_y_a = strata.n/n
  P_R = n_sample/n
  asvar_con <- function(logit.pi_mat){
    num <- exp(logit.pi_mat)
    inf.num <- is.infinite(num)
    
    pi_mat <- rep(1, length(num))
    pi_mat[!inf.num] <- num[!inf.num]/(1+num[!inf.num])
    
    constraint <- sum(p_y_a*pi_mat)-P_R
    return(list(ceq=constraint,c=NULL))
    #      return(constraint)
  }
  
  asvar_obj <- function(logit.pi_mat){
    
    num <- exp(logit.pi_mat)
    inf.num <- is.infinite(num)
    
    #print(num)
    
    pi_mat <- rep(1, length(num))
    pi_mat[!inf.num] <- num[!inf.num]/(1+num[!inf.num])
    
    #print(pi_mat)
    
    info = info_mat.f(pos_vec=pos_vec,grad=grad, cuts=phaseI_res$cuts, sp11=pi_mat[4], sp12=pi_mat[5], sp13=pi_mat[6],
                      sp01=pi_mat[1], sp02=pi_mat[2], sp03=pi_mat[3], para, sc0, sh0, pv)
    
    se_eta1 = sqrt(diag(ginv(info))/n)[target_pos]
    #se_eta1 = sqrt(1/var_eta1/n)
    #se_eta1 <- ifelse(is.nan(se_eta1), 100, se_eta1)
    
    #print(logit.pi_mat)
    print(se_eta1)
    return(se_eta1)
  }
  
  #alp=log(1); et0=et10; et1=et20; gam=ga10
  
  #lb = rep(logit(0.0), length(strata.n))
  opt_design <- solnl(X=ini, objfun = asvar_obj, confun=asvar_con, 
                      tolCon = 1e-05, tolFun = 1e-05)
  print(opt_design$fn)
  pi_opt <- as.vector(round(expit.f(opt_design$par),5))
  inf.num <- as.vector(exp(opt_design$par))
  pi_opt <- ifelse(is.finite(inf.num), pi_opt, 1)
  s_opt <- floor(pi_opt*strata.n)
  print(s_opt)
  
  s_id <-  unlist(sapply(which(s_opt!=0), function(x){as.numeric(sample(as.character(dt_ext$id[dt_ext$group %in% ref.group[x]]),s_opt[x]))}))
  if(sum(s_opt)<n_sample){
    s_id <- c(s_id, as.numeric(sample(as.character(dt_ext$id[!(dt_ext$id %in% s_id)]), n_sample-sum(s_opt))))
  }
  
  dt_use = dt_ext[(dt_ext$id %in% s_id),]
  s_num <- sapply(ref.group, function(g){
    ind = (dt_use$group %in% g)
    sum(ind)
  })
  s.prob <- s_num/strata.n
  s.prob <- ifelse(strata.n==0, 0, s.prob)
  
  res <- list(s_id=s_id,
              s_prob=s.prob)
  
  return(res)
  
}

opt2_ODS_asym_AY.f <- function(scen, grad, n_sample, phaseI_res, ini, para, sc0, sh0, pv){
  
  if(scen=="sc"){
    pos_vec=1:length(para)
  }else if(scen=="sb"){
    pos_vec=(1:length(para))[-6]
  }
  target_pos=3
  
  dt_ext = phaseI_res$dt_ext
  ref.group = phaseI_res$ref.strata
  strata.n = phaseI_res$strata.n
  n <- dim(dt_ext)[1]
  dt_ext$id <- 1:n
  
  p_y_a = strata.n/n
  P_R = n_sample/n
  asvar_con <- function(logit.pi_mat){
    num <- exp(logit.pi_mat)
    inf.num <- is.infinite(num)
    
    pi_mat <- rep(1, length(num))
    pi_mat[!inf.num] <- num[!inf.num]/(1+num[!inf.num])
    
    constraint <- sum(p_y_a*pi_mat)-P_R
    return(list(ceq=constraint,c=NULL))
    #      return(constraint)
  }
  
  asvar_obj <- function(logit.pi_mat){
    
    num <- exp(logit.pi_mat)
    inf.num <- is.infinite(num)
    
    #print(num)
    
    pi_mat <- rep(1, length(num))
    pi_mat[!inf.num] <- num[!inf.num]/(1+num[!inf.num])
    
    #print(pi_mat)
    
    info = info_mat.f(pos_vec=pos_vec,grad=grad, cuts=phaseI_res$cuts, sp11=pi_mat[4], sp12=pi_mat[5], sp13=pi_mat[6],
                      sp01=pi_mat[1], sp02=pi_mat[2], sp03=pi_mat[3], para, sc0, sh0, pv)
    se_eta1 = sqrt(diag(ginv(info))/n)[target_pos]
    
    print(se_eta1)
    return(se_eta1)
  }
  
  #alp=log(1); et0=et10; et1=et20; gam=ga10
  
  #lb = rep(logit(0.0), length(strata.n))
  opt_design <- solnl(X=ini, objfun = asvar_obj, confun=asvar_con, 
                      tolCon = 1e-05, tolFun = 1e-05)
  print(opt_design$fn)
  pi_opt <- as.vector(round(expit.f(opt_design$par),5))
  inf.num <- as.vector(exp(opt_design$par))
  pi_opt <- ifelse(is.finite(inf.num), pi_opt, 1)
  s_opt <- floor(pi_opt*strata.n)
  print(s_opt)
  
  s_id <-  unlist(sapply(which(s_opt!=0), function(x){as.numeric(sample(as.character(dt_ext$id[dt_ext$group %in% ref.group[x]]),s_opt[x]))}))
  if(sum(s_opt)<n_sample){
    s_id <- c(s_id, as.numeric(sample(as.character(dt_ext$id[!(dt_ext$id %in% s_id)]), n_sample-sum(s_opt))))
  }
  
  dt_use = dt_ext[(dt_ext$id %in% s_id),]
  s_num <- sapply(ref.group, function(g){
    ind = (dt_use$group %in% g)
    sum(ind)
  })
  s.prob <- s_num/strata.n
  s.prob <- ifelse(strata.n==0, 0, s.prob)
  
  res <- list(s_id=s_id,
              s_prob=s.prob)
  
  return(res)
  
}

optA_ODS_asym_AY.f <- function(grad, n_sample, phaseI_res, ini, para, sc0, sh0, pv){
  
  pos_vec=1:length(para)
  
  dt_ext = phaseI_res$dt_ext
  ref.group = phaseI_res$ref.strata
  strata.n = phaseI_res$strata.n
  n <- dim(dt_ext)[1]
  dt_ext$id <- 1:n
  
  p_y_a = strata.n/n
  P_R = n_sample/n
  asvar_con <- function(logit.pi_mat){
    num <- exp(logit.pi_mat)
    inf.num <- is.infinite(num)
    
    pi_mat <- rep(1, length(num))
    pi_mat[!inf.num] <- num[!inf.num]/(1+num[!inf.num])
    
    constraint <- sum(p_y_a*pi_mat)-P_R
    return(list(ceq=constraint,c=NULL))
    #      return(constraint)
  }
  
  asvar_obj <- function(logit.pi_mat){
    
    num <- exp(logit.pi_mat)
    inf.num <- is.infinite(num)
    
    #print(num)
    
    pi_mat <- rep(1, length(num))
    pi_mat[!inf.num] <- num[!inf.num]/(1+num[!inf.num])
    
    #print(pi_mat)
    
    #info = E_DR.Abar_AY(pi_vec=pi_mat, cuts=phaseI_res$Abar_cut, beta=beta, eta=eta, P2=P2, param.A=param.A)
    #var_beta1 = sqrt(diag(ginv(info))/n)[2]
    
    info = info_mat.f(pos_vec=pos_vec, grad=grad, cuts=phaseI_res$cuts, sp11=pi_mat[4], sp12=pi_mat[5], sp13=pi_mat[6],
                      sp01=pi_mat[1], sp02=pi_mat[2], sp03=pi_mat[3], para, sc0, sh0, pv)
    # var_eta1 = Veta1.f(cuts=phaseI_res$cuts, sp11=pi_mat[4], sp12=pi_mat[5], sp13=pi_mat[6],
    #                    sp01=pi_mat[1], sp02=pi_mat[2], sp03=pi_mat[3], alp, et0, et1, gam, sc0, sh0)
    
    # var_eta1 = calSS_et1et1.f(cuts=phaseI_res$cuts, sp11=pi_mat[4], sp12=pi_mat[5], sp13=pi_mat[6],
    #                sp01=pi_mat[1], sp02=pi_mat[2], sp03=pi_mat[3], alp, et0, et1, gam, sc0, sh0)
    se_eta1 = sum(diag(ginv(info)/n)[c(3,6)])
    #se_eta1 = sqrt(1/var_eta1/n)
    #se_eta1 <- ifelse(is.nan(se_eta1), 100, se_eta1)
    
    #print(logit.pi_mat)
    print(se_eta1)
    return(se_eta1)
  }
  
  #alp=log(1); et0=et10; et1=et20; gam=ga10
  
  #lb = rep(logit(0.0), length(strata.n))
  opt_design <- solnl(X=ini, objfun = asvar_obj, confun=asvar_con, 
                      tolCon = 1e-05, tolFun = 1e-05)
  print(opt_design$fn)
  pi_opt <- as.vector(round(expit.f(opt_design$par),5))
  inf.num <- as.vector(exp(opt_design$par))
  pi_opt <- ifelse(is.finite(inf.num), pi_opt, 1)
  s_opt <- floor(pi_opt*strata.n)
  print(s_opt)
  s_id <-  unlist(sapply(which(s_opt!=0), function(x){as.numeric(sample(as.character(dt_ext$id[dt_ext$group %in% ref.group[x]]),s_opt[x]))}))
  if(sum(s_opt)<n_sample){
    s_id <- c(s_id, as.numeric(sample(as.character(dt_ext$id[!(dt_ext$id %in% s_id)]), n_sample-sum(s_opt))))
  }
  
  dt_use = dt_ext[(dt_ext$id %in% s_id),]
  s_num <- sapply(ref.group, function(g){
    ind = (dt_use$group %in% g)
    sum(ind)
  })
  s.prob <- s_num/strata.n
  s.prob <- ifelse(strata.n==0, 0, s.prob)
  
  res <- list(s_id=s_id,
              s_prob=s.prob)
  
  return(res)
  
}

optD_ODS_asym_AY.f <- function(grad, n_sample, phaseI_res, ini, para, sc0, sh0, pv){
  
  pos_vec=1:length(para)
  
  dt_ext = phaseI_res$dt_ext
  ref.group = phaseI_res$ref.strata
  strata.n = phaseI_res$strata.n
  n <- dim(dt_ext)[1]
  dt_ext$id <- 1:n
  
  p_y_a = strata.n/n
  P_R = n_sample/n
  asvar_con <- function(logit.pi_mat){
    num <- exp(logit.pi_mat)
    inf.num <- is.infinite(num)
    
    pi_mat <- rep(1, length(num))
    pi_mat[!inf.num] <- num[!inf.num]/(1+num[!inf.num])
    
    constraint <- sum(p_y_a*pi_mat)-P_R
    return(list(ceq=constraint,c=NULL))
    #      return(constraint)
  }
  
  asvar_obj <- function(logit.pi_mat){
    
    num <- exp(logit.pi_mat)
    inf.num <- is.infinite(num)
    
    #print(num)
    
    pi_mat <- rep(1, length(num))
    pi_mat[!inf.num] <- num[!inf.num]/(1+num[!inf.num])
    
    info = info_mat.f(pos_vec=pos_vec,grad=grad, cuts=phaseI_res$cuts, sp11=pi_mat[4], sp12=pi_mat[5], sp13=pi_mat[6],
                      sp01=pi_mat[1], sp02=pi_mat[2], sp03=pi_mat[3], para, sc0, sh0, pv)

    se_eta1=det(ginv(info)[c(3,6),c(3,6)])/n
    
    #print(logit.pi_mat)
    print(se_eta1)
    return(se_eta1)
  }
  
  #alp=log(1); et0=et10; et1=et20; gam=ga10
  
  #lb = rep(logit(0.0), length(strata.n))
  opt_design <- solnl(X=ini, objfun = asvar_obj, confun=asvar_con, 
                      tolCon = 1e-05, tolFun = 1e-05)
  print(opt_design$fn)
  pi_opt <- as.vector(round(expit.f(opt_design$par),5))
  inf.num <- as.vector(exp(opt_design$par))
  pi_opt <- ifelse(is.finite(inf.num), pi_opt, 1)
  s_opt <- floor(pi_opt*strata.n)
  print(s_opt)
  s_id <-  unlist(sapply(which(s_opt!=0), function(x){as.numeric(sample(as.character(dt_ext$id[dt_ext$group %in% ref.group[x]]),s_opt[x]))}))
  if(sum(s_opt)<n_sample){
    s_id <- c(s_id, as.numeric(sample(as.character(dt_ext$id[!(dt_ext$id %in% s_id)]), n_sample-sum(s_opt))))
  }
  
  dt_use = dt_ext[(dt_ext$id %in% s_id),]
  s_num <- sapply(ref.group, function(g){
    ind = (dt_use$group %in% g)
    sum(ind)
  })
  s.prob <- s_num/strata.n
  s.prob <- ifelse(strata.n==0, 0, s.prob)
  
  res <- list(s_id=s_id,
              s_prob=s.prob)
  
  return(res)
  
}
