#' tgvis: Function to estimate and select the optimal number of single effects using profile-likelihood and BIC
#'
#' This function estimates the number of single effects in a locus by combining profile-likelihood methods and Bayesian Information Criterion (BIC) to optimize the model. It includes resampling for estimating standard errors and performs score tests for infinitesimal effects.
#'
#' @param by A vector of Z-scores of the marginal effects from the outcome GWAS.
#' @param bXest A matrix of direct effect estimates based on the Z-scores of tissue-gene pairs.
#' @param LD The LD matrix of variants.
#' @param Noutcome The sample size of the outcome GWAS.
#' @param L_vec A vector of candidate numbers of single effects used in BIC. Default is `c(1:8)`.
#' @param estimate_inf An indicator of whether estimating the infinitesimal effect. Default is F.
#' @param var_inf When estimate_inf = F, the variance of infinitesimal effect (estimated by LDSC possibly). Default is 1e-7.
#' @param estimate_residual_variance An indicator of whether of not estimating the variance of residuals in SuSiE. Default is F.
#' @param residual_variance The residual variance.  Default is 1.
#' @param max_iter The maximum number of iterations for the profile-likelihood algorithm. Default is 50.
#' @param max_eps The convergence tolerance for the profile-likelihood algorithm. Default is 1e-3.
#' @param susie_iter The maximum number of iterations for `susie_rss` within the profile-likelihood algorithm. Default is 50.
#' @param scaled_prior_variance The prior variance of signals in SuSiE. Default is 0.5 which is slightly larger than 0.2 in SuSiE software.
#' @param pip_thres_cred The cumulative PIP threshold for variables in a credible set. Default is 0.95.
#' @param eigen_thres The threshold of eigenvalues for modelling the infinitesimal effect. Default is 1.
#' @param varinf_upper_boundary The upper boundary for the prior variance of infinitesimal effects, multiplied by var(y) to adapt to different locus variances. Default is 0.25.
#' @param varinf_lower_boundary The lower boundary for the prior variance of infinitesimal effects, not multiplied by var(y). Default is 0.001.
#' @param ebic_beta The extended BIC factor for causal effects of tissue-gene pairs and direct causal variants used in BIC computation. Default is 1.
#' @param ebic_upsilon The extended BIC factor for infinitesimal effects used in BIC computation. Default is 1.
#' @param pip_min The minimum PIP threshold for individual causal effects in the profile-likelihood. This is used to specify which tissue-gene pairs and direct causal variants to include in the score test of variance of infinitesimal effects. Default is 0.05.
#' @param pv_thres The p-value threshold for the score test. Default is 0.05.
#' @param pleiotropy_rm A vector of indices specifying which variants should not be considered as having direct causal effects.
#' @param prior_weight_theta A vector of prior weights of gene-tissue pairs, which will be used as input in SuSiE. Default is \code{NULL}.
#' @param prior_weight_gamma A vector of prior weights of direct causal variants, which will be used as input in SuSiE. Default is \code{NULL}.
#' @param standization A indicator of whether standardizing the input when performing SuSiE for fine-mapping causal gene-tissue pairs and direct causal variants. Default is \code{T}.
#'
#' @return A list containing:
#' \item{theta}{The estimated effects for tissue-gene pairs, scaled by the outcome GWAS sample size.}
#' \item{gamma}{The estimated effects for direct causal variants, scaled by the outcome GWAS sample size.}
#' \item{theta.pip}{Posterior inclusion probabilities (PIP) for tissue-gene pairs.}
#' \item{gamma.pip}{Posterior inclusion probabilities (PIP) for direct causal variants.}
#' \item{theta.pratt}{Pratt estimations for tissue-gene pairs.}
#' \item{gamma.pratt}{Pratt estimations for direct causal variants.}
#' \item{theta.cs}{Credible set indicators for tissue-gene pairs.}
#' \item{gamma.cs}{Credible set indicators for direct causal variants.}
#' \item{theta.cs.pip}{PIP within credible sets for tissue-gene pairs.}
#' \item{gamma.cs.pip}{PIP within credible sets for direct causal variants.}
#' \item{upsilon}{The estimated infinitesimal effects.}
#' \item{var.upsilon}{The estimated variance of infinitesimal effects.}
#' \item{fit.causal}{The SuSiE fit object for the causal analysis.}
#' \item{cs.summary}{A summary of the credible sets obtained from the analysis.}
#' \item{Bicvec}{A vector of BIC values for each candidate number of single effects.}
#'
#' @importFrom CppMatrix matrixInverse matrixMultiply matrixVectorMultiply matrixEigen matrixListProduct matrixGeneralizedInverse
#' @importFrom susieR susie_rss susie_get_cs coef.susie
#' @importFrom Matrix Matrix solve
#' @export
#'
tgvis=function(estimate_inf=F,by,bXest,LD,Noutcome,L_vec=c(1:8),
               var_inf=1e-7,estimate_residual_variance=F,
               scaled_prior_variance=0.5,residual_variance=1,
               max_iter=50,max_eps=1e-3,susie_iter=500,pip_thres_cred=0.95,
               eigen_thres=0.999,varinf_upper_boundary=0.25,varinf_lower_boundary=0.001,
               ebic_beta=1,ebic_upsilon=1,pip_min=0.05,pv_thres=0.05,pleiotropy_rm=NULL,
               prior_weight_theta=NULL,prior_weight_gamma=NULL,standization=T){
L_vec=sort(L_vec)
if(estimate_inf==T){
############################## Preparing the data ##############################
n=length(by);p=dim(bXest)[2]
if(is.null(pleiotropy_rm)==T){
pleiotropy_rm=findUniqueNonZeroRows(bXest)
}
pleiotropy.keep=setdiff(1:n,pleiotropy_rm)
bXest=Matrix(bXest,sparse=T)
fiteigen=matrixEigen(LD)
idx=find_cumvar_index(fiteigen$values,thres=eigen_thres)
if(idx<2){
stopifnot("LD region is too small. Try other methods like Mendelian randomization.")
}
Umat=fiteigen$vectors[,1:idx]
Dvec=fiteigen$values[1:idx]
Uy=as.numeric(crossprod(Umat,by))
pratt_vary=sum(Uy^2/Dvec)
varinf_upper_boundary=varinf_upper_boundary*pratt_vary/idx
UXtheta=as.matrix(crossprod(Umat,bXest))*Dvec
if(length(pleiotropy.keep)>0){
UXgamma=t(Umat[pleiotropy.keep,,drop=FALSE])*Dvec
}else{
UXgamma=matrix(0,nrow=idx,ncol=0)
}
UXR=as.matrix(cbind(UXtheta,UXgamma))
XtX=matrixMultiply(t(UXR),UXR*(1/Dvec))
dXtX <- diag(XtX)
dXtX[is.na(dXtX)] <- 1
dXtX[dXtX == 0] <- 1
XtX[is.na(XtX)] <- 0
diag(XtX) <- dXtX
XtXadjust=diag(XtX)
XtX <- cov2cor(XtX)
XtX=(t(XtX)+XtX)/2
if(is.null(prior_weight_theta)==T){
prior_weight_theta=rep(1/p,p)
}
if(is.null(prior_weight_gamma)==T){
prior_weight_gamma=rep(1/length(pleiotropy.keep),length(pleiotropy.keep))
}else{
prior_weight_gamma=prior_weight_gamma[pleiotropy.keep]
}
prior_weights=c(prior_weight_theta,prior_weight_gamma)
########################### Selecting the number of single effects ##############################
Bicvec=L_vec
fit.store=vector("list",length(L_vec))
beta.store=vector("list",length(L_vec))
upsilon.store=vector("list",length(L_vec))
varinf.store=rep(NA_real_,length(L_vec))
fit.prev=NULL
beta.prev=NULL
upsilon.prev=NULL
varinf.prev=NA_real_
for(i in 1:length(L_vec)){
L.current=max(1,L_vec[i])
upsilon=if(is.null(upsilon.prev)) 0*by else upsilon.prev
varinf=if(is.finite(varinf.prev)) varinf.prev else 1
iter=0
error=1
beta=if(is.null(beta.prev)) rep(0,ncol(XtX)) else beta.prev
fit.causal=tgvis_susie_init_for_L(fit.prev,L.current,ncol(XtX),
                                  scaled_prior_variance=scaled_prior_variance,
                                  estimate_residual_variance=estimate_residual_variance)
Hinv=1/(Dvec+1/varinf)
while(error>max_eps&iter<max_iter){
beta1=beta
res.beta=by-matrixVectorMultiply(LD,upsilon)
Xty=c(as.numeric(crossprod(bXest,res.beta)),
      if(length(pleiotropy.keep)>0) res.beta[pleiotropy.keep] else numeric(0))
XtyZ=Xty/sqrt(XtXadjust)
fit.causal=susie_rss(z=XtyZ,R=XtX,n=Noutcome,L=L.current,
                     estimate_prior_method="optim",max_iter=susie_iter,intercept=F,
                     standardize=standization,prior_weights=prior_weights,
                     model_init=fit.causal,estimate_residual_variance=estimate_residual_variance,
                     scaled_prior_variance=scaled_prior_variance,residual_variance=residual_variance)
beta=coef.susie(fit.causal)[-1]*sqrt(Noutcome)/sqrt(XtXadjust)
############# Score test needs to determine the fixed effect #######################
############# We remove the variants in the 95% credible sets with small PIP #######################
causal.cs=group.pip.filter(pip.summary=summary(fit.causal)$var,pip.thres.cred=pip_min)
pip.alive=causal.cs$ind.keep
if(length(pip.alive)>0){
beta[-pip.alive]=0
}else{
beta[]=0
}
eta.beta=tgvis_eta_from_beta(bXest,beta,p,pleiotropy.keep)
res.upsilon=by-matrixVectorMultiply(LD,eta.beta)
outcome=as.numeric(crossprod(Umat,res.upsilon))
########### Performing Score test to pre-remoing infinitesimal effect in each step ###################
active=which(fit.causal$pip>pip_min)
pv=tgvis_inf_test_eigen(res.inf=res.upsilon,D=Dvec,UA=UXR[,active,drop=FALSE])
######################## Performing REML ###############################
upsilon=by*0
if(pv<pv_thres|iter<5){
for(ii in 1:5){
Hinv=1/(Dvec+1/varinf)
upsilon=matrixVectorMultiply(Umat,outcome*Hinv)
for(jj in 1:5){
df=sum(Hinv)
varinf=min((sum(upsilon^2)+df)/idx,varinf_upper_boundary)
pv=ifelse(varinf<varinf_lower_boundary,0.5,pv)
}
}
}
if(iter>5){
error=norm(beta-beta1,"2")/sqrt(length(beta))
}
iter=iter+1
}
df=sum(Dvec*Hinv)*ifelse(sum(abs(upsilon))==0,0,1)
eta.beta=tgvis_eta_from_beta(bXest,beta,p,pleiotropy.keep)
res=by-matrixVectorMultiply(LD,eta.beta+upsilon)
Ures=as.numeric(crossprod(Umat,res))
rss=sum(Ures^2/Dvec)
Bicvec[i]=log(rss)+(log(idx)+ebic_beta*log(dim(XtX)[1]))/idx*L_vec[i]+(ebic_upsilon*log(idx)+log(idx))/idx*df
fit.store[[i]]=fit.causal
beta.store[[i]]=beta
upsilon.store[[i]]=upsilon
varinf.store[i]=varinf
fit.prev=fit.causal
beta.prev=beta
upsilon.prev=upsilon
varinf.prev=varinf
}
################### Reestimating using optimal number of single effect #####################
istar=which.min(Bicvec)
L.best=max(1,L_vec[istar])
upsilon=if(is.null(upsilon.store[[istar]])) 0*by else upsilon.store[[istar]]
varinf=if(is.finite(varinf.store[istar])) varinf.store[istar] else 1
iter=0
error=1
beta=if(is.null(beta.store[[istar]])) rep(0,ncol(XtX)) else beta.store[[istar]]
fit.causal=tgvis_susie_init_for_L(fit.store[[istar]],L.best,ncol(XtX),
                                  scaled_prior_variance=scaled_prior_variance,
                                  estimate_residual_variance=estimate_residual_variance)
Hinv=1/(Dvec+1/varinf)
while(error>max_eps&iter<max_iter){
beta1=beta
res.beta=by-matrixVectorMultiply(LD,upsilon)
Xty=c(as.numeric(crossprod(bXest,res.beta)),
      if(length(pleiotropy.keep)>0) res.beta[pleiotropy.keep] else numeric(0))
XtyZ=Xty/sqrt(XtXadjust)
fit.causal=susie_rss(z=XtyZ,R=XtX,n=Noutcome,L=L.best,
                     estimate_prior_method="optim",max_iter=susie_iter,intercept=F,
                     standardize=standization,prior_weights=prior_weights,
                     model_init=fit.causal,estimate_residual_variance=estimate_residual_variance,
                     scaled_prior_variance=scaled_prior_variance,residual_variance=residual_variance)
beta=coef.susie(fit.causal)[-1]*sqrt(Noutcome)/sqrt(XtXadjust)
causal.cs=group.pip.filter(pip.summary=summary(fit.causal)$var,pip.thres.cred=pip_min)
pip.alive=causal.cs$ind.keep
if(length(pip.alive)>0){
beta[-pip.alive]=0
}else{
beta[]=0
}
eta.beta=tgvis_eta_from_beta(bXest,beta,p,pleiotropy.keep)
res.upsilon=by-matrixVectorMultiply(LD,eta.beta)
outcome=as.numeric(crossprod(Umat,res.upsilon))
active=which(fit.causal$pip>pip_min)
pv=tgvis_inf_test_eigen(res.inf=res.upsilon,D=Dvec,UA=UXR[,active,drop=FALSE])
upsilon=by*0
if(pv<pv_thres|iter<5){
for(ii in 1:10){
Hinv=1/(Dvec+1/varinf)
upsilon=matrixVectorMultiply(Umat,outcome*Hinv)
for(jj in 1:5){
df=sum(Hinv)
varinf=min((sum(upsilon^2)+df)/idx,varinf_upper_boundary)
pv=ifelse(varinf<varinf_lower_boundary,runif(1,pv_thres,1),pv)
}
}
}
if(iter>5){
error=norm(beta-beta1,"2")/sqrt(length(beta))
}
iter=iter+1
}
var.upsilon=varinf
####################################### Getting the variables in 95% credible sets ####################################
causal.cs=group.pip.filter(pip.summary=summary(fit.causal)$var,pip.thres.cred=pip_thres_cred)
pip.alive=causal.cs$ind.keep
pip.remove=setdiff(1:ncol(XtX),pip.alive)
###################################### Preparing the results ################################
if(length(pip.alive)>0){
pip.remove=setdiff(1:ncol(XtX),pip.alive)
gammatheta=coef.susie(fit.causal)[-1]*sqrt(Noutcome)/sqrt(XtXadjust)
gammatheta[pip.remove]=0
gamma=rep(0,n)
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
theta=theta.pip=theta.se=theta.pratt=theta.cs=theta.cs.pip=rep(0,p)
gamma=gamma.pip=gamma.se=gamma.pratt=gamma.cs=gamma.cs.pip=rep(0,n)
names(theta)=names(theta.pip)=colnames(bXest)
names(gamma)=names(gamma.pip)=rownames(bXest)
}
A=list(theta=theta,gamma=gamma,theta.pip=theta.pip,gamma.pip=gamma.pip,theta.pratt=theta.pratt,gamma.pratt=gamma.pratt,theta.cs=theta.cs,gamma.cs=gamma.cs,theta.cs.pip=theta.cs.pip,gamma.cs.pip=gamma.cs.pip,upsilon=upsilon,var.upsilon=var.upsilon,pv.upsilon=pv,fit.causal=fit.causal,cs.summary=causal.cs,Bicvec=Bicvec)
Asummary=summary.tgvis(fit=A,bXest=bXest)
A$Summary=Asummary
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
}else{
A=tgvis_fast(by,bXest,LD,Lvec=L_vec,n_outcome=Noutcome,var_inf=var_inf,eigen_thres=eigen_thres,pip_thres_cred=pip_thres_cred,susie_iter=susie_iter,estimate_residual_variance=estimate_residual_variance,scaled_prior_variance=scaled_prior_variance,residual_variance=residual_variance,prior_weight_theta=prior_weight_theta,prior_weight_gamma=prior_weight_gamma,ebic_factor=ebic_beta,pleiotropy_rm=pleiotropy_rm)
return(A)
}
}
