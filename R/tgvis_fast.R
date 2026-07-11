tgvis_fast=function(by,bXest,LD,Lvec=c(1:8),n_outcome,var_inf=0.15/1e6,eigen_thres=0.999,pip_thres_cred=0.95,susie_iter=500,estimate_residual_variance=F,scaled_prior_variance=0.5,residual_variance=1,prior_weight_theta=NULL,prior_weight_gamma=NULL,ebic_factor=1,pleiotropy_rm=NULL){
######################### basic infromation ####################################
Lvec=sort(Lvec)
m=length(by)
p=ncol(bXest)
if(is.null(pleiotropy_rm)==T){
pleiotropy_rm=findUniqueNonZeroRows(bXest)
}
pleiotropy.keep=setdiff(1:m,pleiotropy_rm)
bXest=Matrix(bXest,sparse=T)
fiteigen=matrixEigen(LD)
idx=find_cumvar_index(fiteigen$values,thres=eigen_thres)
if(idx<2){
stopifnot("LD region is too small. Try other methods like Mendelian randomization.")
}
U=fiteigen$vectors[,1:idx]
D=fiteigen$values[1:idx]
lambda=n_outcome*var_inf
theta_weight=1/(residual_variance*D+lambda*D^2)
pratt_weight=1/D
###################### XtX, Xty, and yty #######################################
UXtheta=as.matrix(crossprod(U,bXest))*D
if(length(pleiotropy.keep)>0){
UXgamma=t(U[pleiotropy.keep,,drop=FALSE])*D
}else{
UXgamma=matrix(0,nrow=idx,ncol=0)
}
UXR=as.matrix(cbind(UXtheta,UXgamma))
Uy=as.numeric(crossprod(U,by))
XtX=matrixMultiply(t(UXR),UXR*theta_weight)
dXtX <- diag(XtX)
dXtX[is.na(dXtX)] <- 1
dXtX[dXtX == 0] <- 1
XtX[is.na(XtX)] <- 0
diag(XtX) <- dXtX
XtXadjust=diag(XtX)
XtX <- cov2cor(XtX)
XtX=(t(XtX)+XtX)/2
Xty=matrixVectorMultiply(t(UXR),Uy*theta_weight)
Xty=Xty/sqrt(XtXadjust)
pratt_vary=sum(Uy^2*pratt_weight)
############################## fit SuSiE #######################################
if(is.null(prior_weight_gamma)==T){
if(length(pleiotropy.keep)>0){
prior_weight_gamma=rep(1/length(pleiotropy.keep),length(pleiotropy.keep))
}else{
prior_weight_gamma=numeric(0)
}
}else{
prior_weight_gamma=prior_weight_gamma[pleiotropy.keep]
}
if(is.null(prior_weight_theta)==T){
prior_weight_theta=rep(1/p,p)
}
prior_weights=c(prior_weight_theta,prior_weight_gamma)
Bicvec=Lvec
fit.store=vector("list",length(Lvec))
fit.prev=NULL
for(i in 1:length(Lvec)){
fit.init=tgvis_susie_init_for_L(fit.prev,Lvec[i],ncol(XtX),
                                scaled_prior_variance=scaled_prior_variance,
                                estimate_residual_variance=estimate_residual_variance)
fiti=susie_rss(z=Xty,R=XtX,n=n_outcome,L=Lvec[i],residual_variance=residual_variance,estimate_residual_variance=estimate_residual_variance,estimate_prior_method="optim",max_iter=susie_iter,model_init=fit.init)
theta=coef(fiti)[-1]*sqrt(n_outcome)/sqrt(XtXadjust)
fitted=as.numeric(bXest%*%theta[1:p])
if(length(pleiotropy.keep)>0){
fitted[pleiotropy.keep]=fitted[pleiotropy.keep]+theta[-c(1:p)]
}
res=by-matrixVectorMultiply(LD,fitted)
Ures=as.numeric(crossprod(U,res))
rss=sum(Ures^2*theta_weight)
Bicvec[i]=log(rss)+(log(idx)+log(ncol(XtX))*ebic_factor)/idx*Lvec[i]
fit.store[[i]]=fiti
fit.prev=fiti
}
istar=which.min(Bicvec)
fit.init=tgvis_susie_init_for_L(fit.store[[istar]],Lvec[istar],ncol(XtX),
                                scaled_prior_variance=scaled_prior_variance,
                                estimate_residual_variance=estimate_residual_variance)
fit.causal=susie_rss(z=Xty,R=XtX,n=n_outcome,L=Lvec[istar],residual_variance=residual_variance,estimate_residual_variance=estimate_residual_variance,estimate_prior_method="optim",max_iter=susie_iter,model_init=fit.init)
causal.cs=group.pip.filter(pip.summary=summary(fit.causal)$var,pip.thres.cred=pip_thres_cred)
pip.alive=causal.cs$ind.keep
if(length(pip.alive)>0){
pip.remove=setdiff(1:ncol(XtX),pip.alive)
gammatheta=coef.susie(fit.causal)[-1]*sqrt(n_outcome)/sqrt(XtXadjust)
gammatheta[pip.remove]=0
gamma=LD[,1]*0
gamma[pleiotropy.keep]=gammatheta[-c(1:p)]
theta=gammatheta[1:p]
gamma.pip=gamma*0
gamma.pip[pleiotropy.keep]=fit.causal$pip[-c(1:p)]
theta.pip=fit.causal$pip[1:p]
gammatheta.cs=causal.cs$cs
gammatheta.cs.pip=causal.cs$cs.pip
theta.cs=gammatheta.cs[1:p]
theta.cs.pip=gammatheta.cs.pip[1:p]
gamma.cs=gamma.cs.pip=gamma*0
gamma.cs[pleiotropy.keep]=gammatheta.cs[-c(1:p)]
gamma.cs.pip[pleiotropy.keep]=gammatheta.cs.pip[-c(1:p)]
gamma.pratt=by*gamma/pratt_vary
theta.pratt=as.numeric(crossprod(bXest,by))*theta/pratt_vary
gamma.pratt[is.na(gamma.pratt)]=0
theta.pratt[is.na(theta.pratt)]=0
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
Asummary$cs.logBF <- rep(NA_real_, nrow(Asummary))
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
