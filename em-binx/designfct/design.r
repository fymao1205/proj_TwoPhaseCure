
# design function 

# phase I stratification: stratify on del, discretized T, or (del, discretized T)
phaesI_strat.f <- function(pcuts, dt, design.factor){
  
  style = paste0(design.factor, collapse="_")
  
  if(is.null(pcuts)) pcuts <- c(1/3, 2/3)
  
  c_cut <- quantile(dt$time, p=pcuts)
  cbar <- cut(dt$time, c(-0.01,c_cut,Inf), labels = paste0(1:(length(c_cut)+1)))
  
  dt$cbar <- cbar
  
  if(style=="cbar"){
    
    z <- xtabs(~cbar, dt)
    group <- factor(paste0(dt$cbar))
    group.mat <- expand.grid(1:3)
    #ref.group <- paste0(group.mat[,1])
    ref.group=1:3
    strata.n <- as.data.frame(z)$Freq[!ind] # give size of each strata
    
  }
  
  if(style=="del"){
    
    z <- xtabs(~del, dt)
    group <- factor(paste0(dt$del))
    group.mat <- expand.grid(0:1)
    #ref.group <- paste0(group.mat)
    ref.group <- 0:1
    strata.n <- as.data.frame(z)$Freq # give size of each strata
    
  }
  
  if(style=="cbar_del"){
    
    z <- xtabs(~cbar+del, dt)
    group <- factor(paste0(dt$cbar, dt$del))
    group.mat <- expand.grid(1:3, 0:1)
    ref.group <- paste0(group.mat[,1], group.mat[,2])
    strata.n <- as.data.frame(z)$Freq # give size of each strata
    
  }
  
  dt$group <- group
  
  res <- list(strata.n=strata.n,
              dt_ext = dt,
              ref.strata = ref.group,
              cuts=c_cut)
  
  return(res)
  
}

design_strat.f <- function(n2samp, design="srs", phaseI_strat){
  
  dt_ext = phaseI_strat$dt_ext
  ref.group = phaseI_strat$ref.strata
  strata.n = phaseI_strat$strata.n
  n <- sum(strata.n)
  dt_ext$id <- 1:n
  
  if(design=="srs"){
    
    s_id <- as.numeric(sample(as.character(dt_ext$id),n2samp))
    
  }
  
  if(design=="bal"){
    
    num_strata <- length(strata.n)
    n_circ_vec <- pmin(n2samp/num_strata, strata.n)
    sum_n_circ <- sum(n_circ_vec)
    
    left <- n2samp - sum_n_circ
    
    while(left >0 ){
      left_bal <- strata.n - n_circ_vec
      add_bal <- numeric(num_strata)
      for(i in 1:num_strata){
        div <- sum(left_bal >0)
        add_bal[i] <- ifelse(left_bal[i] < left/div,left_bal[i],left/div)
      }
      n_circ_vec <- n_circ_vec + add_bal
      left <- n2samp - sum(n_circ_vec)
    }
    
    s_bal <- floor(n_circ_vec)
    
    s.prob <- s_bal/strata.n
    
    s_id <- unlist(sapply(which(s_bal!=0), function(x){as.numeric(sample(as.character(dt_ext$id[dt_ext$group%in% ref.group[x]]),s_bal[x]))}))

  }
  
  
  dt_use = dt_ext[(dt_ext$id %in% s_id),]
  s_num <- sapply(ref.group, function(g){
    ind = (dt_use$group %in% g)
    sum(ind)
  })
  
  if(design=="srs"){
    s.prob <- rep(n2samp/sum(strata.n), length(strata.n))
  }else{
    s.prob <- s_num/strata.n #
  }  
  
  res <- list(s_id=s_id,
              s_prob=s.prob,
              s_m=s_num)
  
} 

rkd_ext_opt.f <- function(k, n2samp, qzi, phaseI_strat){
  
  dt_ext = phaseI_strat$dt_ext
  ref.group = phaseI_strat$ref.strata
  strata.n = phaseI_strat$strata.n
  n <- sum(strata.n)
  #dt_ext$id <- 1:n
  
  indV0=which(dt_ext$v==0)
  nV0 = length(indV0)
  indV1 = which(dt_ext$v == 1)
  nV1 = length(indV1)
  
  #### optimal sampling ####################################################################################
  
  order_resi0 = order(qzi[indV0])
  order_resi1 = order(qzi[indV1])
  
  best_k = k#round(n2samp*0.15)
  phase2_id = c(indV0[order_resi0[1:best_k]], indV0[order_resi0[(nV0-best_k+1):nV0]], 
                indV1[order_resi1[1:(n2samp/2-best_k)]], indV1[order_resi1[(nV1-(n2samp/2-best_k)+1):nV1]])
  mart.opt = qzi[phase2_id]
  simZ.opt = dt_ext$v[phase2_id]
  best_var = (var(mart.opt[which(simZ.opt == 0)])*sum(simZ.opt==0)+var(mart.opt[which(simZ.opt == 1)])*sum(simZ.opt==1))/n2samp
  
  for (k in (best_k+1):(n2samp/2-best_k))
  {
    phase2_id = c(indV0[order_resi0[1:k]], indV0[order_resi0[(nV0-k+1):nV0]], 
                  indV1[order_resi1[1:(n2samp/2-k)]], indV1[order_resi1[(nV1-(n2samp/2-k)+1):nV1]])
    mart.opt = qzi[phase2_id]
    simZ.opt = dt_ext$v[phase2_id]
    tmp_var = (var(mart.opt[which(simZ.opt == 0)])*sum(simZ.opt==0)+var(mart.opt[which(simZ.opt == 1)])*sum(simZ.opt==1))/n2samp
    if (tmp_var > best_var) {
      best_k = k
      best_var = tmp_var
    }
  }
  
  s_id = c(indV0[order_resi0[1:best_k]], indV0[order_resi0[(nV0-best_k+1):nV0]], 
           indV1[order_resi1[1:(n2samp/2-best_k)]], indV1[order_resi1[(nV1-(n2samp/2-best_k)+1):nV1]])
  
  s_id = dt_ext$id[s_id]
  
  s_dt = subset(dt_ext, id %in% s_id)
  s_prob <- sapply(1:length(ref.group), function(g){
    
    if(any(s_dt$group %in% ref.group[g])){
      
      dim(subset(s_dt, group %in% ref.group[g]))[1]/strata.n[g]
      
    }else{
      
      0
    }
    
  })
  
  s_m <- sapply(1:length(ref.group), function(g){
    
    if(any(s_dt$group %in% ref.group[g])){
      
      nrow(subset(s_dt, group %in% ref.group[g]))
      
    }else{
      
      0
    }
    
  })
  
  res <- list(s_prob=s_prob,
              s_id=s_id,
              s_m=s_m)
  
  return(res)
  
  
}

rkd_ext_opt_2st.f <- function(k, n2samp, qzi, phaseI_strat, s_id.a=NULL){
  
  dt_ext = phaseI_strat$dt_ext
  ref.group = phaseI_strat$ref.strata
  strata.n = phaseI_strat$strata.n
  n <- sum(strata.n)
  #dt_ext$id <- 1:n
  
  indV0=which(dt_ext$v==0)
  nV0 = length(indV0)
  indV1 = which(dt_ext$v == 1)
  nV1 = length(indV1)
  
  order_resi0 = order(qzi[indV0])
  order_resi1 = order(qzi[indV1])
  
  if(!is.null(s_id.a)){
    
    dt_ext_b=subset(dt_ext, !(dt_ext$id %in% s_id.a))
    qzi_b=qzi[!(dt_ext$id %in% s_id.a)]
    
    indV0_b=which(dt_ext_b$v==0)
    nV0_b = length(indV0_b)
    indV1_b = which(dt_ext_b$v == 1)
    nV1_b = length(indV1_b)
    
    #### optimal sampling in phase IIb ####################################################################################
    
    order_resi0_b = order(qzi_b[indV0_b])
    order_resi1_b = order(qzi_b[indV1_b])
    
    best_k = k#round(n2samp*0.15)
    phase2_id_b_ind = c(indV0_b[order_resi0_b[1:best_k]], indV0_b[order_resi0_b[(nV0_b-best_k+1):nV0_b]], 
                  indV1_b[order_resi1_b[1:(n2samp/2-best_k)]], indV1_b[order_resi1_b[(nV1_b-(n2samp/2-best_k)+1):nV1_b]])
    phase2_id_b=dt_ext_b$id[phase2_id_b_ind]
    
    phase2_id=c(phase2_id_b, s_id.a)
    mart.opt = qzi[phase2_id]
    simZ.opt = dt_ext$v[phase2_id]
    best_var = (var(mart.opt[which(simZ.opt == 0)])*sum(simZ.opt==0)+var(mart.opt[which(simZ.opt == 1)])*sum(simZ.opt==1))/n2samp
    
    for (k in (best_k+1):(n2samp/2-best_k))
    {
      phase2_id_b_ind = c(indV0_b[order_resi0_b[1:k]], indV0_b[order_resi0_b[(nV0_b-k+1):nV0_b]], 
                    indV1_b[order_resi1_b[1:(n2samp/2-k)]], indV1_b[order_resi1_b[(nV1_b-(n2samp/2-k)+1):nV1_b]])
      phase2_id_b=dt_ext_b$id[phase2_id_b_ind]
      phase2_id=c(phase2_id_b, s_id.a)
      mart.opt = qzi[phase2_id]
      simZ.opt = dt_ext$v[phase2_id]
      tmp_var = (var(mart.opt[which(simZ.opt == 0)])*sum(simZ.opt==0)+var(mart.opt[which(simZ.opt == 1)])*sum(simZ.opt==1))/n2samp
      #print(tmp_var)
      if (tmp_var > best_var) {
        best_k = k
        best_var = tmp_var
      }
    }
    
    s_id_b_ind = c(indV0_b[order_resi0_b[1:best_k]], indV0_b[order_resi0_b[(nV0_b-best_k+1):nV0_b]], 
             indV1_b[order_resi1_b[1:(n2samp/2-best_k)]], indV1_b[order_resi1_b[(nV1_b-(n2samp/2-best_k)+1):nV1_b]])
    s_id_b=dt_ext_b$id[s_id_b_ind]
    s_id=c(s_id_b, s_id.a)
    
  }else{
    best_k = k#round(n2samp*0.15)
    phase2_id = c(indV0[order_resi0[1:best_k]], indV0[order_resi0[(nV0-best_k+1):nV0]], 
                  indV1[order_resi1[1:(n2samp/2-best_k)]], indV1[order_resi1[(nV1-(n2samp/2-best_k)+1):nV1]])
    mart.opt = qzi[phase2_id]
    simZ.opt = dt_ext$v[phase2_id]
    best_var = (var(mart.opt[which(simZ.opt == 0)])*sum(simZ.opt==0)+var(mart.opt[which(simZ.opt == 1)])*sum(simZ.opt==1))/n2samp
    
    for (k in (best_k+1):(n2samp/2-best_k))
    {
      phase2_id = c(indV0[order_resi0[1:k]], indV0[order_resi0[(nV0-k+1):nV0]], 
                    indV1[order_resi1[1:(n2samp/2-k)]], indV1[order_resi1[(nV1-(n2samp/2-k)+1):nV1]])
      mart.opt = qzi[phase2_id]
      simZ.opt = dt_ext$v[phase2_id]
      tmp_var = (var(mart.opt[which(simZ.opt == 0)])*sum(simZ.opt==0)+var(mart.opt[which(simZ.opt == 1)])*sum(simZ.opt==1))/n2samp
      if (tmp_var > best_var) {
        best_k = k
        best_var = tmp_var
      }
    }
    
    s_id = c(indV0[order_resi0[1:best_k]], indV0[order_resi0[(nV0-best_k+1):nV0]], 
             indV1[order_resi1[1:(n2samp/2-best_k)]], indV1[order_resi1[(nV1-(n2samp/2-best_k)+1):nV1]])
  }
  
  s_id = dt_ext$id[s_id]
  
  s_dt = subset(dt_ext, id %in% s_id)
  s_prob <- sapply(1:length(ref.group), function(g){
    
    if(any(s_dt$group %in% ref.group[g])){
      
      dim(subset(s_dt, group %in% ref.group[g]))[1]/strata.n[g]
      
    }else{
      
      0
    }
    
  })
  
  s_m <- sapply(1:length(ref.group), function(g){
    
    if(any(s_dt$group %in% ref.group[g])){
      
      nrow(subset(s_dt, group %in% ref.group[g]))
      
    }else{
      
      0
    }
    
  })
  
  res <- list(s_prob=s_prob,
              s_id=s_id,
              s_m=s_m)
  
  return(res)
  
  
}








