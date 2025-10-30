tgvis_fast=function(by,bXest,LD,Lvec=c(1:8),n_outcome,var_inf=0.15/1e6,eigen_thres=0.999,pip_thres_cred=0.95,susie_iter=500,estimate_residual_variance=F,scaled_prior_variance=0.5,residual_variance=1,prior_weight_theta=NULL,prior_weight_gamma=NULL,ebic_factor=1){
######################### basic infromation ####################################
m=length(by)
p=ncol(bXest)
bXest=Matrix(bXest,sparse=T)
fiteigen=matrixEigen(LD)
idx=find_cumvar_index(fiteigen$values,thres=eigen_thres)
if(idx<2){
stopifnot("LD region is too small. Try other methods like Mendelian randomization.")
}
U=fiteigen$vectors[,1:idx]
D=fiteigen$values[1:idx]
lambda=n_outcome*var_inf
Theta=matrixMultiply(U,t(U)*(1/(residual_variance*D+lambda*D^2)))
Theta_Pratt=matrixMultiply(U,t(U)*(1/D))
###################### XtX, Xty, and yty #######################################
XR=as.matrix(cbind(LD%*%bXest,LD))
XRinv=matrixMultiply(t(XR),Theta)
XtX=matrixMultiply(XRinv,XR)
dXtX <- diag(XtX)
dXtX[is.na(dXtX)] <- 1
dXtX[dXtX == 0] <- 1
XtX[is.na(XtX)] <- 0
diag(XtX) <- dXtX
XtXadjust=diag(XtX)
XtX <- cov2cor(XtX)
XtX=(t(XtX)+XtX)/2
Xty=matrixVectorMultiply(XRinv,by)
Xty=Xty/sqrt(XtXadjust)
############################## fit SuSiE #######################################
if(is.null(prior_weight_gamma)==T){
prior_weight_gamma=rep(1/m,m)
}
if(is.null(prior_weight_theta)==T){
prior_weight_theta=rep(1/p,p)
}
prior_weights=c(prior_weight_theta,prior_weight_gamma)
Bicvec=Lvec
for(i in 1:length(Lvec)){
fiti=susie_rss(z=Xty,R=XtX,n=n_outcome,L=Lvec[i],residual_variance=1,estimate_residual_variance=estimate_residual_variance,estimate_prior_method="EM",max_iter=susie_iter)
theta=coef(fiti)[-1]*sqrt(n_outcome)/sqrt(XtXadjust)
res=by-matrixVectorMultiply(XR,theta)
rss=sum(res*matrixVectorMultiply(Theta,res))
Bicvec[i]=log(rss)+(log(idx)+log(m+p)*ebic_factor)/idx*Lvec[i]
}
istar=which.min(Bicvec)
fit.causal=susie_rss(z=Xty,R=XtX,n=n_outcome,L=Lvec[istar],residual_variance=1,estimate_residual_variance=estimate_residual_variance,estimate_prior_method="EM",max_iter=susie_iter)
causal.cs=group.pip.filter(pip.summary=summary(fit.causal)$var,pip.thres.cred=pip_thres_cred)
pip.alive=causal.cs$ind.keep
if(length(pip.alive)>0){
pip.remove=setdiff(1:ncol(XtX),pip.alive)
gammatheta=coef.susie(fit.causal)[-1]*sqrt(n_outcome)/sqrt(XtXadjust)
gammatheta[pip.remove]=0
gamma=LD[,1]*0
gamma=gammatheta[-c(1:p)]
theta=gammatheta[1:p]
gamma.pip=gamma*0
gamma.pip=fit.causal$pip[-c(1:p)]
theta.pip=fit.causal$pip[1:p]
gammatheta.cs=causal.cs$cs
gammatheta.cs.pip=causal.cs$cs.pip
theta.cs=gammatheta.cs[1:p]
theta.cs.pip=gammatheta.cs.pip[1:p]
gamma.cs=gamma.cs.pip=gamma*0
gamma.cs=gammatheta.cs[-c(1:p)]
gamma.cs.pip=gammatheta.cs.pip[-c(1:p)]
gamma.pratt=prattestimation_gamma(by=by,LD=LD,Theta=Theta_Pratt,gamma=gamma)
theta.pratt=prattestimation_theta(by=by,bXest=bXest,LD=LD,Theta=Theta_Pratt,theta=theta)
names(theta)=names(theta.pip)=colnames(bXest)
names(gamma)=names(gamma.pip)=rownames(bXest)
}else{
theta=theta.pip=theta.pratt=theta.cs=theta.cs.pip=rep(0,p)
gamma=gamma.pip=gamma.pratt=gamma.cs=gamma.cs.pip=rep(0,m)
names(theta)=names(theta.pip)=colnames(bXest)
names(gamma)=names(gamma.pip)=rownames(bXest)
}
A=list(theta=theta,gamma=gamma,theta.pip=theta.pip,gamma.pip=gamma.pip,theta.pratt=theta.pratt,gamma.pratt=gamma.pratt,theta.cs=theta.cs,gamma.cs=gamma.cs,theta.cs.pip=theta.cs.pip,gamma.cs.pip=gamma.cs.pip,fit.causal=fit.causal,cs.summary=causal.cs,Bicvec=Bicvec)
Asummary=summary.tgvis(fit=A,bXest=bXest)
bb <- summary(fit.causal)$cs
Asummary$cs.logBF <- NA_real_
if (!is.null(bb) && "cs" %in% names(bb) && "cs_log10bf" %in% names(bb)) {
  for (i in seq_len(nrow(Asummary))) {
    cs_i <- Asummary$cs[i]
    if (!is.na(cs_i) && cs_i %in% bb$cs) {
      Asummary$cs.logBF[i] <- bb$cs_log10bf[bb$cs == cs_i]
    }
  }
}
need_cols <- c("variable","cs","cs.logBF","cs.pratt","cs.pip",
               "xqtl","type","estimate","pip","pratt")
Asummary=Asummary[,need_cols]
A$Summary=Asummary
return(A)
}
