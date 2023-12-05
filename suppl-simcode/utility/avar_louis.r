

CureCox.Louis.varMatrix=function(lam_vec, tk_vec, pid, z, del, time, x1, X, Xp, V, beta, eta, gam, w1, w2)
{
  gradn=CureCox.Louis.Egradn(lam_vec, tk_vec, pid, z, del, time, x1, X, Xp, V, beta, eta, gam, w1, w2)
  
  Hess=CureCox.Louis.Ehess(lam_vec, tk_vec, pid, z, del, time, x1, X, Xp, V, beta, eta, gam, w1, w2)
  
  grad.sq=CureCox.Louis.Egradsq(lam_vec, tk_vec, pid, z, del, time, x1, X, Xp, V, beta, eta, gam, w1, w2)
  
  #  print(max(abs(apply(gradn,2,sum))))
  #  print(min(eigen(-Hess[1:5,1:5])$value))
  #  print(min(eigen(-Hess)$value))
  
  N=length(unique(pid)) # Sample size
  ker=rep(1,N)%*%t(rep(1,N))
  diag(ker)=0
  Louis=-Hess-grad.sq-t(gradn)%*%ker%*%gradn
  
  
  #out=t(gradn)%*%ker%*%gradn+grad.sq
  
  out=as.matrix(solve(Louis)[1:(length(eta)+length(beta)+length(gam)),1:(length(eta)+length(beta)+length(gam))])
  return(out)
}

# Gradient^2 of log-likelihood
CureCox.Louis.Egradsq=function(lam_vec, tk_vec, pid, z, del, time, x1, X, Xp, V, beta, eta, gam, w1, w2)
{
  K=length(tk_vec)        # number of unique event times
  N=length(time)
  
  # Logisitc Probability
  covZ=drop(eta%*% t(Xp))
  expitcovZ=drop(expit.f(covZ))
  
  # Hazard Proportion
  phaz=drop(exp(beta%*% t(X)))
  
  Lam0.f <- stepfun(tk_vec, c(0,cumsum(lam_vec)))
  Lam0=Lam0.f(time)
  
  # Hazard and Density
  Lnk=phaz%*%t(Lam0)
  Snk=exp(-Lnk)
  #lamk=diff(c(0,Lam0(tk)))
  col.lamk=(rep(1,N)%*%t(lam_vec))
  
  # Indicators
  I.Xi.tk= (time%*%t(rep(1,K))) >= (rep(1,N)%*%t(tk_vec))
  I.Xieqtk= (time%*%t(rep(1,K))) == (rep(1,N)%*%t(tk_vec))
  
  # Covariates 
  covX=drop(gam%*% t(V))
  expitcovX=expit.f(covX)
  
  # Lambda terms
  wll.I=(I.Xieqtk*del/col.lamk- I.Xi.tk*phaz)*z
  #wll.II=PTitk/col.lamk-PTijtk*phaz
  #wll.III=PTitk/col.lamk^2-2*PTitk/col.lamk*phaz+PTijtk*phaz^2
  
  # beta terms
  wb.I=(del-phaz*Lam0)*z
  
  # eta terms
  we.I=(z-expitcovZ)
  
  # gam terms
  wg.I=(x1-expitcovX)
  
  # Grad.lamlam
  wll.I=(I.Xieqtk/col.lamk-I.Xi.tk*phaz)*z
  wll.offdiag=apply(-w1*w2*wll.I*phaz,2,sum)
  maxid=pmax(rep(1:K,K),rep(1:K,rep(K,K)))
  wll.block=matrix(wll.offdiag[maxid],K,K)
  diag(wll.block)=apply(wll.I^2*w1*w2,2,sum)
  grad.lamlam= wll.block
  
  # Grad.elam
  welam=wll.I*we.I*w1*w2 
  grad.elam=t(Xp)%*%welam
  
  # Grad.glam
  wglam=wll.I*wg.I*w1*w2
  grad.glam=t(V)%*%wglam
  
  # Grad.ab
  #if(length(beta)>0)
  #{
    
    # Grad.blam
    wblam= wll.I * (wb.I)*w1*w2 
    grad.blam=t(X)%*%wblam 
    
    # Grad.bb
    wbb=(wb.I)^2*w1*w2
    grad.bb=t(X) %*% diag(wbb) %*% X
    
    # Grad.be
    web=we.I*wb.I*w1*w2
    grad.be=t(X) %*% diag(web)%*% Xp
    
     # Grad.bg
    wgb=wg.I*wb.I*w1*w2
    grad.bg=t(X) %*% diag(wgb)%*% V
    
    
  # }else{
  #   grad.ab = matrix(NA,length(a),0)
  #   grad.bb = matrix(NA,0,0)
  #   grad.blam = matrix(NA,0,K)
  # }
  
  
  # Grad.ee
  wee=we.I^2*w1*w2
  grad.ee=t(Xp) %*% diag(wee)%*% Xp
  
  # Grad.eg
  weg=we.I*wg.I*w1*w2
  grad.eg=t(Xp) %*% diag(weg)%*% V
  
  # Grad.gg
  wgg=wg.I^2*w1*w2
  grad.gg=t(V) %*% diag(wgg)%*% V
  
  # out=cbind(grad.bb,grad.be, grad.bg)
  # out=rbind(out,cbind(t(grad.be),grad.ee, grad.eg))
  # out=rbind(out,cbind(t(grad.bg),t(grad.eg), grad.gg))
  # out=cbind(out,rbind((grad.blam),(grad.elam), grad.glam))
  # out=rbind(out,cbind(t(grad.blam),t(grad.elam),t(grad.glam),grad.lamlam))
  out=cbind(grad.ee,grad.eg, t(grad.be))
  out=rbind(out,cbind(t(grad.eg),grad.gg, t(grad.bg)))
  out=rbind(out,cbind((grad.be),(grad.bg), grad.bb))
  out=cbind(out,rbind((grad.elam),(grad.glam), grad.blam))
  out=rbind(out,cbind(t(grad.elam),t(grad.glam),t(grad.blam),grad.lamlam))
  
  return(out)
}


# Hessian of log-likelihood
CureCox.Louis.Ehess=function(lam_vec, tk_vec, pid, z, del, time, x1, X, Xp, V, beta, eta, gam, w1, w2)
{
  K=length(tk_vec)        # number of unique event times
  N=length(time)
  nk=apply(outer(time[del==1], tk_vec, "=="),2,sum) # counting ties, usually all 1
  
  # Logisitc Probability
  covZ=drop(eta%*% t(Xp))
  expitcovZ=drop(expit.f(covZ))
  
  # Hazard Proportion
  phaz=drop(exp(beta%*% t(X)))
  
  Lam0.f <- stepfun(tk_vec, c(0,cumsum(lam_vec)))
  Lam0=Lam0.f(time)
  
  # Hazard and Density
  Lnk=phaz%*%t(Lam0)
  Snk=exp(-Lnk)
  #lamk=diff(c(0,Lam0(tk)))
  col.lamk=(rep(1,N)%*%t(lam_vec))
  
  # Indicators
  I.Xi.tk= (time%*%t(rep(1,K))) >= (rep(1,N)%*%t(tk_vec))
  I.Xieqtk= (time%*%t(rep(1,K))) == (rep(1,N)%*%t(tk_vec))
  
  # Covariates 
  covX=drop(gam%*% t(V))
  expitcovX=expit.f(covX)
  
  # Lambda terms
  #wll.I=(I.Xieqtk*del/col.lamk- I.Xi.tk*phaz)*z
  #wll.II=PTitk/col.lamk-PTijtk*phaz
  #wll.III=PTitk/col.lamk^2-2*PTitk/col.lamk*phaz+PTijtk*phaz^2
  
  # beta terms
  #wb.I=(del-phaz*Lam0)*z
  
  # eta terms
  #we.I=(z-expitcovZ)
  
  # gam terms
  #wg.I=(x1-expitcovX)
  
  # Hessian for bb
  wbb=(w1*w2*z*(-Lam0*phaz))#+Eyn*phaz*Lam0
  Hbb=t(X) %*% diag(wbb)%*% X
  
  # Hessian for ee
  wee=-expitcovZ*(1-expitcovZ)*w1*w2
  Hee=t(Xp)%*% diag(wee) %*% Xp
  
  # Hessian for gg
  wgg=-expitcovX*(1-expitcovX)*w1*w2
  Hgg=t(V)%*% diag(wgg) %*% V
  
  # Hessian for blamb
  wblam=(-I.Xi.tk)*phaz*z*w1*w2
  Hblam=t(X)%*%wblam
  
  # Hessian for lamlam
  #Hll=diag((-(nk)/lam_vec^2))
  col.lamk2=(rep(1,N)%*%t(lam_vec^2))
  Hll=diag(colSums(-I.Xieqtk/col.lamk2*z*w1*w2))
  
  temp=rbind(cbind(Hbb,Hblam),cbind(t(Hblam),Hll))
  
  return(Matrix::bdiag(Hee,Hgg,temp))
}

# Gradient of log-likelihood
CureCox.Louis.Egradn=function(lam_vec, tk_vec, pid, z, del, time, x1, X, Xp, V, beta, eta, gam, w1, w2)
{
  K=length(tk_vec)        # number of unique event times
  N=length(time)
  nk=apply(outer(time[del==1], tk_vec, "=="),2,sum) # counting ties, usually all 1
  
  # Logisitc Probability
  covZ=drop(eta%*% t(Xp))
  expitcovZ=drop(expit.f(covZ))
  
  # Hazard Proportion
  phaz=drop(exp(beta%*% t(X)))
  
  Lam0.f <- stepfun(tk_vec, c(0,cumsum(lam_vec)))
  Lam0=Lam0.f(time)
  
  # Hazard and Density
  Lnk=phaz%*%t(Lam0)
  Snk=exp(-Lnk)
  #lamk=diff(c(0,Lam0(tk)))
  col.lamk=(rep(1,N)%*%t(lam_vec))
  
  # Indicators
  I.Xi.tk= (time%*%t(rep(1,K))) >= (rep(1,N)%*%t(tk_vec))
  I.Xieqtk= (time%*%t(rep(1,K))) == (rep(1,N)%*%t(tk_vec))
  
  # Covariates 
  covX=drop(gam%*% t(V))
  expitcovX=expit.f(covX)
  
  # Lambda terms
  wll.I=(I.Xieqtk*del/col.lamk- I.Xi.tk*phaz)*z*w1*w2
  gradn.lam=wll.I
  
  # beta terms
  wb.I=(del-phaz*Lam0)*z*w1*w2
  gradn.b=X*wb.I
  
  # eta terms
  we.I=(z-expitcovZ)*w1*w2
  gradn.e=Xp*we.I
  
  # gam terms
  wg.I=(x1-expitcovX)*w1*w2
  gradn.g=V*wg.I
  
  res0=cbind(gradn.e, gradn.g, gradn.b,gradn.lam)
  dt=data.frame(cbind(pid=pid,res0))
  
  res=aggregate(.~pid, data=dt, sum)[,-1]
  
  return(as.matrix(res))
}



#test = survival::Surv(time=dat$time, event=dat$del, type="right")
