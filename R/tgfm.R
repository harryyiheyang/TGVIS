#' tgfm: Function to perform eQTL fine-mapping and causal tissue-gene pair estimation
#'
#' This function performs fine-mapping of eQTLs and estimates causal tissue-gene pairs using summary statistics, linkage disequilibrium (LD) matrix, and resampling methods. It returns a list of estimates and associated statistics, including credible sets, PIP resampling, and Pratt estimations.
#'
#' @param by A vector of Z-scores of the marginal effects from the outcome GWAS.
#' @param bX A matrix of Z-scores of marginal eQTL effect estimates for tissue-gene pairs.
#' @param LD The LD matrix of variants.
#' @param Nvec A vector representing the sample sizes, with the first element for the outcome GWAS and the subsequent elements for eQTL studies.
#' @param L.eqtl The number of single effects used in the eQTL fine-mapping step. Default is 5.
#' @param L.causal The number of single effects used for tissue-gene pairs and direct causal variants. Default is 10.
#' @param pip.thres.cred The cumulative PIP threshold for variables in a credible set. Default is 0.95.
#' @param eqtl.sampling.time The number of resampling iterations for eQTL effect estimation. Default is 100.
#' @param causal.sampling.time The number of resampling iterations for causal effect estimation. Default is 100.
#' @param eqtl.thres A threshold for the minimum individual PIPs in each credible set during eQTL fine-mapping. This threshold is used to pre-remove variables with very low PIPs to avoid overly large credible sets. Default is 0.05.
#' @param susie.iter The maximum number of iterations for `susie_rss` within the profile-likelihood algorithm. Default is 50.
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
#' \item{fit.causal}{The SuSiE fit object for the causal analysis.}
#' \item{cs.summary}{A summary of the credible sets obtained from the analysis.}
#' \item{pip.causal.samping}{The PIP resampling results from causal effect estimation.}
#' \item{estimate.causal.sampling}{The effect size resampling results from causal effect estimation.}
#' \item{bXest}{The matrix of resampled eQTL effects used in causal effect estimation.}
#'
#' @importFrom CppMatrix matrixInverse matrixMultiply matrixVectorMultiply matrixEigen matrixListProduct matrixGeneralizedInverse
#' @importFrom susieR susie_rss susie_get_cs coef.susie
#' @importFrom Matrix Matrix solve
#' @export
#'
tgfm=function(by,bX,LD,Nvec,L.eqtl=5,L.causal=10,pip.thres.cred=0.5,eqtl.sampling.time=100,causal.sampling.time=100,eqtl.thres=0.05,susie.iter=1000){
n=length(by);p=dim(bX)[2]
Theta=matrixInverse(LD)
eX=bX*0
eQTLList=list()
###################### Estimating eQTL Effect Step ######################
for(i in 1:p){
errorindicator <- FALSE
indx <- which(bX[,i] != 0)
eQTLList[[i]]=list(alpha=matrix(0.5,1,p),mu=matrix(0,1,p),mu2=matrix(1,1,p),index.causal=1,indx=indx)
tryCatch({
a <- LD[indx, indx]
fit <- susie_rss(z = bX[indx, i], R = a, n = Nvec[i + 1], L = L.eqtl,  estimate_prior_method="EM")
####################### We don't consider the credible set including too many variables ################
index.causal = intersect(unique(susie_get_cs_index(fit)),which(fit$pip>eqtl.thres))
eQTLList[[i]]=list(alpha=fit$alpha,mu=fit$mu,mu2=fit$mu2,index.causal=index.causal,indx=indx)
}, error = function(e){
cat("Error in iteration", i, ": ", e$message, "\n")
errorindicator <- TRUE
})
if(errorindicator) next
}
######################### Sampling Step ########################
ZArray= array(0,c(n,p,causal.sampling.time))
ZAll= matrix(0,n,p)
colnames(ZAll)=colnames(bX)
rownames(ZAll)=rownames(bX)
for (j in 1:p){
z = bX[eQTLList[[j]]$indx, 1] * 0
indj = unique(eQTLList[[j]]$index.causal)
z = tgfm.resampling(alpha = eQTLList[[j]]$alpha, mu = eQTLList[[j]]$mu, mu2 = eQTLList[[j]]$mu2, sampling = causal.sampling.time*eqtl.sampling.time) * sqrt(Nvec[j + 1])
zall = colMeans(z)
zall[abs(zall)<eqtl.thres]=0
if(length(indj)>0){
zall[-indj]=0
}
ZAll[eQTLList[[j]]$indx,j]=zall
for(jj in 1:causal.sampling.time){
indjj=c((eqtl.sampling.time*(jj-1)+1):(eqtl.sampling.time*jj))
zjj=z[indjj,]
if(eqtl.sampling.time>1){
zjj=colMeans(zjj)
}
zjj[abs(zjj)<eqtl.thres]=0
if(length(indj)>0){
zjj[-indj]=0
}
ZArray[eQTLList[[j]]$indx,j,jj]=zjj
}
}
#################################### Estimation Causal Effect Step ############################
eXi = ZAll
Xty = c(t(eXi) %*% by, by)
XtX = Matrix::bdiag(matrixMultiply(t(eXi) %*% LD, eXi), LD)
XtX = as.matrix(XtX)
XtX[1:p, -c(1:p)] = matrixMultiply(t(eXi), LD)
XtX[-c(1:p), c(1:p)] = t(XtX[1:p, -c(1:p)])
dXtX = diag(XtX); dXtX[is.na(dXtX)] = 1; dXtX[dXtX == 0] = 1;
XtX[is.na(XtX)] = 0; diag(XtX) = dXtX
Xadjust=diag(XtX)
XtX=cov2cor(XtX)
XtX=(t(XtX)+XtX)/2
XtyZ=Xty/sqrt(Xadjust)
prior.weight.theta=rep(1/p,p)
prior.weight.gamma=rep(1/n,n)
prior_weights=c(prior.weight.theta,prior.weight.gamma)
fit.causal = susie_rss(z=XtyZ,R=XtX,n=Nvec[1], L = L.causal, residual_variance = 1, estimate_prior_method="EM", prior_weights=prior_weights, intercept=F,max_iter=susie.iter)
############################# The second resmaping step ##################################
AA = AB = matrix(0, causal.sampling.time,n+p)
for (i in 1:causal.sampling.time) {
eXi = ZArray[,,i]
Xty = c(t(eXi) %*% by, by)
XtX = Matrix::bdiag(matrixMultiply(t(eXi) %*% LD, eXi), LD)
XtX = as.matrix(XtX)
XtX[1:p, -c(1:p)] = matrixMultiply(t(eXi), LD)
XtX[-c(1:p), c(1:p)] = t(XtX[1:p, -c(1:p)])
dXtX = diag(XtX); dXtX[is.na(dXtX)] = 1; dXtX[dXtX == 0] = 1;
XtX[is.na(XtX)] = 0; diag(XtX) = dXtX
Xadjusti=diag(XtX)
XtX=cov2cor(XtX)
XtX=(t(XtX)+XtX)/2
XtyZ=Xty/sqrt(Xadjusti)
fit.causali = susie_rss(z=XtyZ,R=XtX,n=Nvec[1],L=L.causal,estimate_prior_method="EM",s_init=fit.causal,prior_weights=prior_weights,intercept=F,max_iter=10)
AA[i,] = fit.causali$pip
AB[i,] = coef.susie(fit.causali)[-1]
}
################################# Preparing the result #####################################
fit.causal$pip.resampling=colMeans(AA)
fit.causal.sampling=susie.resampling(alpha=fit.causal$alpha,mu=fit.causal$mu,mu2=fit.causal$mu2)
fit.causal$beta.se=apply(AB,2,sd)
causal.cs=group.pip.filter(pip.summary=summary(fit.causal)$var,pip.thres.cred=pip.thres.cred)
pip.alive=causal.cs$ind.keep
pip.remove=setdiff(1:ncol(XtX),pip.alive)
thetagamma=coef.susie(fit.causal)[-1]
thetagamma[pip.remove]=0
thetagamma.se=fit.causal$beta.se
thetagamma.se[pip.remove]=0
thetagamma.pip=colMeans(AA)
thetagamma.pip[pip.remove]=0
thetagamma.cs=causal.cs$cs
thetagamma.cs.pip=causal.cs$cs.pip
thetagamma.cs.pip[pip.remove]=0
thetagamma.cs[pip.remove]=0
theta=thetagamma[1:p]
gamma=thetagamma[-c(1:p)]
theta.se=thetagamma.se[1:p]
gamma.se=thetagamma.se[-c(1:p)]
theta.pip=thetagamma.pip[1:p]
gamma.pip=thetagamma.pip[-c(1:p)]
theta.cs=thetagamma.cs[1:p]
gamma.cs=thetagamma.cs[-c(1:p)]
theta.cs.pip=thetagamma.cs.pip[1:p]
gamma.cs.pip=thetagamma.cs.pip[-c(1:p)]
gamma.pratt=prattestimation(by=by,bXest=diag(length(by)),LD=LD,Theta=Theta,theta=gamma*sqrt(Nvec[1]))
theta.pratt=prattestimation(by=by,bXest=ZAll,LD=LD,Theta=Theta,theta=theta*sqrt(Nvec[1])/sqrt(Xadjust[1:p]))
names(theta)=names(theta.pip)=colnames(bX)
names(gamma)=names(gamma.pip)=rownames(bX)
A=list(theta=theta*sqrt(Nvec[1])/sqrt(Xadjust[1:p]),gamma=gamma*sqrt(Nvec[1]),theta.se=theta.se*sqrt(Nvec[1])/sqrt(Xadjust[1:p]),gamma.se=gamma.se*sqrt(Nvec[1]),theta.pip=theta.pip,gamma.pip=gamma.pip,theta.pratt=theta.pratt,gamma.pratt=gamma.pratt,theta.cs=theta.cs,gamma.cs=gamma.cs,theta.cs.pip=theta.cs.pip,gamma.cs.pip=gamma.cs.pip,fit.causal=fit.causal,cs.summary=causal.cs,pip.causal.samping=AA,estimate.causal.sampling=AB,bXest=ZAll)
return(A)
}
