#' tgvis: Function to estimate and select the optimal number of single effects using profile-likelihood and BIC
#'
#' This function estimates the number of single effects in a locus by combining profile-likelihood methods and Bayesian Information Criterion (BIC) to optimize the model. It includes resampling for estimating standard errors and performs score tests for infinitesimal effects.
#'
#' @param by A vector of Z-scores of the marginal effects from the outcome GWAS.
#' @param bXest A matrix of direct effect estimates based on the Z-scores of tissue-gene pairs.
#' @param LD The LD matrix of variants.
#' @param Noutcome The sample size of the outcome GWAS.
#' @param L.causal.vec A vector of candidate numbers of single effects used in BIC. Default is `c(1:8)`.
#' @param max.iter The maximum number of iterations for the profile-likelihood algorithm. Default is 50.
#' @param max.eps The convergence tolerance for the profile-likelihood algorithm. Default is 1e-3.
#' @param susie.iter The maximum number of iterations for `susie_rss` within the profile-likelihood algorithm. Default is 50.
#' @param pip.thres.cred The cumulative PIP threshold for variables in a credible set. Default is 0.95.
#' @param eigen.thres The threshold of eigenvalues for modelling the infinitesimal effect. Default is 1.
#' @param varinf.upper.boundary The upper boundary for the prior variance of infinitesimal effects, multiplied by var(y) to adapt to different locus variances. Default is 0.25.
#' @param varinf.lower.boundary The lower boundary for the prior variance of infinitesimal effects, not multiplied by var(y). Default is 0.001.
#' @param ebic.beta The extended BIC factor for causal effects of tissue-gene pairs and direct causal variants used in BIC computation. Default is 1.
#' @param ebic.upsilon The extended BIC factor for infinitesimal effects used in BIC computation. Default is 1.
#' @param pip.min The minimum PIP threshold for individual causal effects in the profile-likelihood. This is used to specify which tissue-gene pairs and direct causal variants to include in the score test of variance of infinitesimal effects. Default is 0.05.
#' @param pv.thres The p-value threshold for the score test. Default is 0.05.
#' @param pleiotropy.rm A vector of indices specifying which variants should not be considered as having direct causal effects.
#' @param prior.weight.theta A vector of prior weights of gene-tissue pairs, which will be used as input in SuSiE. Default is \code{NULL}.
#' @param prior.weight.gamma A vector of prior weights of direct causal variants, which will be used as input in SuSiE. Default is \code{NULL}.

#'
#' @return A list containing:
#' \item{theta}{The estimated effects for tissue-gene pairs, scaled by the outcome GWAS sample size.}
#' \item{gamma}{The estimated effects for direct causal variants, scaled by the outcome GWAS sample size.}
#' \item{theta.se}{Standard errors of the tissue-gene pair effects.}
#' \item{gamma.se}{Standard errors of the direct causal variant effects.}
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
tgvis=function(by,bXest,LD,Noutcome,L.causal.vec=c(1:8),max.iter=50,max.eps=1e-3,susie.iter=500,pip.thres.cred=0.95,eigen.thres=1,varinf.upper.boundary=0.25,varinf.lower.boundary=0.001,ebic.beta=1,ebic.upsilon=1,pip.min=0.05,pv.thres=0.05,pleiotropy.rm=NULL,prior.weight.theta=NULL,prior.weight.gamma=NULL){
############################## Preparing the data ##############################
n=length(by);p=dim(bXest)[2]
if(is.null(pleiotropy.rm)==T){
pleiotropy.rm=findUniqueNonZeroRows(bXest)
}
pleiotropy.keep=setdiff(1:n,pleiotropy.rm)
Theta=matrixInverse(LD)
varinf.upper.boundary=varinf.upper.boundary*sum(by*(Theta%*%by))/n
XR=cbind(matrixMultiply(LD,bXest),LD[,pleiotropy.keep])
XR1=cbind(bXest,diag(n)[,pleiotropy.keep])
XtX=matrixMultiply(t(XR),XR1)
dXtX <- diag(XtX)
dXtX[is.na(dXtX)] <- 1
dXtX[dXtX == 0] <- 1
XtX[is.na(XtX)] <- 0
diag(XtX) <- dXtX
XtXadjust=diag(XtX)
XtX <- cov2cor(XtX)
XtX=(t(XtX)+XtX)/2
fiteigen=matrixEigen(LD)
Umat=fiteigen$vectors
Dvec=fiteigen$values
Kthres=eigen_cumsum(Dvec,eigen.thres)
Umat=Umat[,1:Kthres]
Dvec=Dvec[1:Kthres]
LD2=matrixMultiply(Umat,t(Umat)*(Dvec^2))
if(is.null(prior.weight.theta)==T){
prior.weight.theta=rep(1/p,p)
}
if(is.null(prior.weight.gamma)==T){
prior.weight.gamma=rep(1/length(pleiotropy.keep),length(pleiotropy.keep))
}else{
prior.weight.gamma=prior.weight.gamma[pleiotropy.keep]
}
prior_weights=c(prior.weight.theta,prior.weight.gamma)
########################### Selecting the number of single effects ##############################
Bicvec=L.causal.vec
for(i in 1:length(L.causal.vec)){
upsilon=0*by
varinf=1
iter=0
error=1
beta=XR[1,]
fit.causal=NULL
while(error>max.eps&iter<max.iter){
beta1=beta
res.beta=by-matrixVectorMultiply(LD,upsilon)
Xty=c(t(bXest)%*%res.beta,res.beta[pleiotropy.keep])
XtyZ=Xty/sqrt(XtXadjust)
fit.causal=susie_rss(z=XtyZ,R=XtX,n=Noutcome,L=max(1,L.causal.vec[i]),estimate_prior_method="EM",max_iter=susie.iter,intercept=F,standardize=F,prior_weights=prior_weights,s_init=fit.causal)
beta=coef.susie(fit.causal)[-1]*sqrt(Noutcome)/sqrt(XtXadjust)
############# Score test needs to determine the fixed effect #######################
############# We remove the variants in the 95% credible sets with small PIP #######################
causal.cs=group.pip.filter(pip.summary=summary(fit.causal)$var,pip.thres.cred=pip.min)
pip.alive=causal.cs$ind.keep
beta[-pip.alive]=0
res.upsilon=by-matrixVectorMultiply(XR,beta)
outcome=matrixVectorMultiply(t(Umat),res.upsilon)
########### Performing Score test to pre-remoing infinitesimal effect in each step ###################
pv=inf.test(res.inf=res.upsilon,LD=LD,LD2=LD2,Theta=Theta,A=XR[,which(fit.causal$pip>pip.min)])
######################## Performing REML ###############################
upsilon=by*0
if(pv<pv.thres|iter<5){
for(ii in 1:5){
Hinv=1/(Dvec+1/varinf)
upsilon=matrixVectorMultiply(Umat,outcome*Hinv)
for(jj in 1:5){
df=sum(Hinv)
varinf=min((sum(upsilon^2)+df)/n,varinf.upper.boundary)
pv=ifelse(varinf<varinf.lower.boundary,0.5,pv)
}
}
}
if(iter>5){
error=norm(beta-beta1,"2")/sqrt(length(beta))
}
iter=iter+1
}
df=sum(Dvec*Hinv)*ifelse(sum(abs(upsilon))==0,0,1)
res=by-matrixVectorMultiply(XR,beta)-matrixVectorMultiply(LD,upsilon)
rss=sum(res*matrixVectorMultiply(Theta,res))
Bicvec[i]=log(rss)+(log(n)+ebic.beta*log(dim(XtX)[1]))/n*L.causal.vec[i]+(ebic.upsilon*log(n)+log(n))/n*df
}
################### Reestimating using optimal number of single effect #####################
istar=which.min(Bicvec)
upsilon=0*by
varinf=1
iter=0
error=1
fit.causal=NULL
while(error>max.eps&iter<max.iter){
beta1=beta
res.beta=by-matrixVectorMultiply(LD,upsilon)
Xty=c(t(bXest)%*%res.beta,res.beta[pleiotropy.keep])
XtyZ=Xty/sqrt(XtXadjust)
fit.causal=susie_rss(z=XtyZ,R=XtX,n=Noutcome,L=max(1,L.causal.vec[istar]),estimate_prior_method="EM",max_iter=susie.iter,intercept=F,standardize=F,prior_weights=prior_weights,s_init=fit.causal)
beta=coef.susie(fit.causal)[-1]*sqrt(Noutcome)/sqrt(XtXadjust)
causal.cs=group.pip.filter(pip.summary=summary(fit.causal)$var,pip.thres.cred=pip.min)
pip.alive=causal.cs$ind.keep
beta[-pip.alive]=0
res.upsilon=by-matrixVectorMultiply(XR,beta)
outcome=matrixVectorMultiply(t(Umat),res.upsilon)
pv=inf.test(res.inf=res.upsilon,LD=LD,LD2=LD2,Theta=Theta,A=XR[,which(fit.causal$pip>pip.min)])
upsilon=by*0
if(pv<pv.thres|iter<5){
for(ii in 1:10){
Hinv=1/(Dvec+1/varinf)
upsilon=matrixVectorMultiply(Umat,outcome*Hinv)
for(jj in 1:5){
df=sum(Hinv)
varinf=min((sum(upsilon^2)+df)/n,varinf.upper.boundary)
pv=ifelse(varinf<varinf.lower.boundary,runif(1,pv.thres,1),pv)
}
}
}
if(iter>5){
error=norm(beta-beta1,"2")/sqrt(length(beta))
}
iter=iter+1
}
var.upsilon=varinf
################################# Otaining SE using resampling #######################
fit.causal.sampling=susie.resampling(alpha=fit.causal$alpha,mu=fit.causal$mu,mu2=fit.causal$mu2)
fit.causal$beta.se=fit.causal.sampling$sd
####################################### Getting the variables in 95% credible sets ####################################
causal.cs=group.pip.filter(pip.summary=summary(fit.causal)$var,pip.thres.cred=pip.thres.cred)
pip.alive=causal.cs$ind.keep
pip.remove=setdiff(1:ncol(XtX),pip.alive)
###################################### Preparing the results ################################
if(length(pip.alive)>0){
pip.remove=setdiff(1:ncol(XtX),pip.alive)
gammatheta=coef.susie(fit.causal)[-1]
gammatheta[pip.remove]=0
gamma=LD[,1]*0
gamma[pleiotropy.keep]=gammatheta[-c(1:p)]
theta=gammatheta[1:p]
gammatheta.se=fit.causal$beta.se
gammatheta.se[pip.remove]=0
gamma.se=gamma
gamma.se[pleiotropy.keep]=gammatheta.se[-c(1:p)]
theta.se=gammatheta.se[1:p]
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
gamma.pratt=prattestimation(by=by,bXest=diag(length(by)),LD=LD,Theta=Theta,theta=gamma*sqrt(Noutcome))
theta.pratt=prattestimation(by=by,bXest=bXest,LD=LD,Theta=Theta,theta=theta*sqrt(Noutcome)/sqrt(XtXadjust[1:p]))
names(theta)=names(theta.pip)=colnames(bXest)
names(gamma)=names(gamma.pip)=rownames(bXest)
}else{
theta=theta.pip=theta.se=theta.pratt=theta.cs=theta.cs.pip=rep(0,p)
gamma=gamma.pip=gamma.se=gamma.pratt=gamma.cs=gamma.cs.pip=rep(0,n)
names(theta)=names(theta.pip)=colnames(bXest)
names(gamma)=names(gamma.pip)=rownames(bXest)
}
A=list(theta=theta*sqrt(Noutcome)/sqrt(XtXadjust[1:p]),gamma=gamma*sqrt(Noutcome),theta.se=theta.se*sqrt(Noutcome)/sqrt(XtXadjust[1:p]),gamma.se=gamma.se*sqrt(Noutcome),theta.pip=theta.pip,gamma.pip=gamma.pip,theta.pratt=theta.pratt,gamma.pratt=gamma.pratt,theta.cs=theta.cs,gamma.cs=gamma.cs,theta.cs.pip=theta.cs.pip,gamma.cs.pip=gamma.cs.pip,upsilon=upsilon,var.upsilon=var.upsilon,pv.upsilon=pv,fit.causal=fit.causal,cs.summary=causal.cs,Bicvec=Bicvec)
return(A)
}
