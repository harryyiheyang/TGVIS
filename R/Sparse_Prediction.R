#' @title Sparse Prediction of xQTL Effects using SuSiE and CARMA
#' @description
#' This function performs informative xQTL selection and sparse prediction using SuSiE and CARMA.
#'
#' @param bX A matrix of effect estimates from the exposure GWAS.
#' @param bXse A matrix of standard errors of effect estimates from the exposure GWAS.
#' @param LD The linkage disequilibrium (LD) matrix.
#' @param xQTL.method The method used in purifying the xQTLs. SuSiE or CARMA can be used here, where the latter can be more accurate but much most computationally costly. Defaults is SuSiE.
#' @param xQTL.selection.rule The method for purifying informative xQTLs within each credible set. Options include "minimum_pip", which selects all variables with PIPs exceeding a specified threshold, and "top_K", which ensures at least K variables are selected based on their PIP ranking. Defaults to "top_K".
#' @param top_K The maximum number of variables selected in each credible sets. Defaults to 1.
#' @param xQTL.pip.min The minimum empirical PIP used in purifying variables in each credible set. Defaults to \code{0.2}.
#' @param xQTL.pip.thres The threshold of individual PIP when selecting xQTL. Defaults to \code{0.5}.
#' @param xQTL.max.L When choosing \code{"SuSiE"}, the maximum number of L in estimating the xQTL effects. Defaults to 10.
#' @param xQTL.cred.thres When choosing \code{"SuSiE"}, the minimum empirical posterior inclusion probability (PIP) used in getting credible sets of xQTL selection. Defaults to \code{0.95}.
#' @param xQTL.Nvec When choosing \code{"SuSiE"}, the vector of sample sizes of exposures.
#' @param xQTL.weight When choosing \code{"SuSiE"}, the vector of weights used in specifying the prior weights of SuSiE. Defaults to \code{NULL}.
#' @param outlier.switch When choosing \code{"CARMA"}, an indicator of whether turning on outlier detection. Defaults to \code{F}.
#' @param Annotation When choosing \code{"CARMA"}, the annotation matrix of SNP. Default is NULL.
#' @param output.labels When choosing \code{"CARMA"}, output directory where output will be written while CARMA is running. Defaults to \code{NULL}, meaning that a temporary folder will be created and automatically deleted upon completion of the computation.
#' @param carma.iter When choosing \code{"CARMA"}, the maximum iterations for EM algorithm to run. Defaults to 5.
#' @param carma.inner.iter When choosing \code{"CARMA"}, the maximum iterations for Shotgun algorithm to run per iteration within EM algorithm. Defaults to 5.
#' @param xQTL.max.num When choosing \code{"CARMA"}, the maximum number of causal variants assumed per locus, which is similar to the number of single effects in SuSiE. Defaults to 10.
#' @param carma.epsilon.threshold When choosing \code{"CARMA"}, the convergence threshold measured by average of Bayes factors. Defaults to \code{1e-3}.
#' @param ridge.diff A ridge.parameter on the differences of causal effect estimate in one credible set. Defaults to \code{100}.
#' @importFrom MASS rlm ginv
#' @importFrom CppMatrix matrixInverse matrixMultiply matrixVectorMultiply matrixEigen matrixListProduct
#' @importFrom Matrix Matrix solve chol bdiag
#' @importFrom susieR susie_suff_stat coef.susie susie susie_rss susie_get_cs
#'
#' @return A list containing the results of the MRBEEX analysis using different methods:
#' \describe{
#'   \item{\code{bXest}}{A matrix where each column represents \code{R * beta_PLS} for each exposure.}
#'   \item{\code{bXestse}}{A matrix where each column contains the standard errors of \code{R * beta_PLS} for each exposure.}
#'   \item{\code{bXest0}}{A matrix where each column represents \code{beta_PLS} for each exposure.}
#'   \item{\code{bXestse0}}{A matrix where each column contains the standard errors of \code{beta_PLS} for each exposure.}
#'   \item{\code{xQTLfitList}}{A list containing the xQTL selection results, which can be used in \code{CisMRBEEX} by setting \code{xQTLfitList = xQTLfitList}.}
#' }
#' @export

Sparse_Prediction=function(bX,bXse,LD,xQTL.Nvec,ridge.diff=100,
                           xQTL.method="SuSiE",xQTL.selection.rule="top_K",
                           top_K=1,xQTL.pip.min=0.2,
                           xQTL.max.L=10,xQTL.cred.thres=0.95,
                           xQTL.pip.thres=0.5,xQTL.weight=NULL,
                           outlier.switch=T,Annotation=NULL,output.labels=NULL,
                           carma.iter=5,carma.inner.iter=5,xQTL.max.num=10,
                           carma.epsilon.threshold=1e-3){
p=ncol(bX)
m=nrow(bX)
############################ Selection Part ##################################
if(xQTL.method=="SuSiE"){
xQTLfitList=list()
if(is.null(xQTL.weight)==T){
xQTL.weight=rep(1,m)
}
pb <- txtProgressBar(min = 0, max = p, style = 3)
for(i in 1:p){
fit=susie_rss(z=bX[,i]/bXse[,i],R=LD,n=xQTL.Nvec[i],L=xQTL.max.L,max_iter=1000,prior_weights=xQTL.weight)
cs_found <- susie_get_cs(fit, coverage = xQTL.cred.thres)$cs
if (length(cs_found) > 0) {
fit=susie_rss(z=bX[,i]/bXse[,i],R=LD,n=xQTL.Nvec[i],L=length(susie_get_cs(fit,coverage=xQTL.cred.thres)$cs),max_iter=1000,prior_weights=xQTL.weight)
}
xQTLfitList[[i]]=fit
setTxtProgressBar(pb, i)
}
close(pb)
}

if(xQTL.method=="CARMA"){
if (!requireNamespace("CARMA", quietly = TRUE)) {
stop("Package 'CARMA' is required for xQTL.method='CARMA', but is not installed. ",
"Please install it with install.packages('CARMA').")
}
z.list=ld.list=w.list=lambda.list=list()
for(i in 1:p){
z.list[[i]]=bX[,i]/bXse[,i]
ld.list[[i]]=LD
w.list[[i]]=Annotation
lambda.list[[i]]=1
}
is.detect=F
if(is.null(output.labels)==T){
output.labels=tempfile("my_tmpdir_")
dir.create(output.labels)
is.delect=T
}
if(is.null(Annotation)==F){
fitxQTL=CARMA::CARMA(z.list,ld.list,w.list,lambda.list,outlier.switch=outlier.switch,num.causal=xQTL.max.num,printing.log=F,all.iter=carma.iter,all.inner.iter=carma.inner.iter,epsilon.threshold=carma.epsilon.threshold,output.labels=output.labels)
}else{
fitxQTL=CARMA::CARMA(z.list,ld.list,lambda.list = lambda.list,outlier.switch=outlier.switch,num.causal=xQTL.max.num,printing.log=F,all.iter=carma.iter,all.inner.iter=carma.inner.iter,epsilon.threshold=carma.epsilon.threshold,output.labels=output.labels)
}
if(is.delect==T){
unlink(output.labels, recursive = TRUE, force = TRUE)
}
xQTLfitList=fitxQTL
}

############################### Prediction Part ################################
bXest=bX
bXest0=bXestse0=bX*0
bXestse=bXestse0=matrix(1000,m,p)
if(xQTL.method=="SuSiE"){
for(i in 1:p){
fit=xQTLfitList[[i]]
if(xQTL.selection.rule=="top_K"){
indj=top_K_pip(summary(fit)$vars,top_K=top_K,pip.min.thres=xQTL.pip.min,xQTL.pip.thres=xQTL.pip.thres)
}else{
causal.cs=group.pip.filter.xQTL(pip.summary = summary(fit)$var,xQTL.cred.thres=xQTL.cred.thres,xQTL.pip.thres=xQTL.pip.min)
indj=union(causal.cs$ind.keep,which(fit$pip>xQTL.pip.thres))
}
if(length(indj)>0){
betaj=coef.susie(fit)[-1]
betaj[-indj]=0
Diff=generate_block_matrix(summary(fit)$vars,rep(1,m),betaj)
LDj=LD[indj,indj]
Thetaj=solve(LDj+Diff[indj,indj]*1e3)
bXest0[indj,i]=as.vector(Thetaj%*%(bX[indj,i]/bXse[indj,i]))*bXse[indj,i]
Thetajj=LD*0
Thetajj[indj,indj]=Thetaj
bXestse0[,i]=sqrt(diag(Thetajj))*bXse[,i]
LDjj=LD%*%Thetajj%*%LD
bXest[,i]=as.vector(LD%*%(bXest0[,i]/bXse[,i]))*bXse[,i]
bXestse[,i]=sqrt(diag(LDjj))*bXse[,i]
}else{
bXest[,i]=bX[,i]
bXestse[,i]=bXse[,i]
}
}
}

if(xQTL.method=="CARMA"){
fitxQTL=xQTLfitList
for(i in 1:p){
sumstat.result = data.frame(variable=c(1:nrow(bX)),variable_prob = fitxQTL[[i]]$PIPs, cs = rep(0,nrow(bX)))
if(length(fitxQTL[[i]]$`Credible set`[[2]])!=0){
for(l in 1:length(fitxQTL[[i]]$`Credible set`[[2]])){
sumstat.result$cs[fitxQTL[[i]]$`Credible set`[[2]][[l]]]=l
}
}
if(xQTL.selection.rule=="top_K"){
indj=top_K_pip(sumstat.result,top_K=top_K,pip.min.thres=xQTL.pip.min,xQTL.pip.thres=xQTL.pip.thres)
}else{
indj=which(sumstat.result$cs!=0&sumstat.result$pip>xQTL.pip.min)
}
if(length(indj)>0){
fitjj=susie_rss(z=bX[,i]/bXse[,i],R=LD,n=xQTL.Nvec[i],L=1+sum(sumstat.result$cs>0),max_iter=100)
betaj=coef.susie(fitjj)[-1]
Diff=generate_block_matrix_CARMA(sumstat.result,rep(1,m),betaj)
Diff=Diff[indj,indj]
LDj=LD[indj,indj]
Thetaj=solve(LDj+Diff*1e3)
bXest0[indj,i]=as.vector(Thetaj%*%(bX[indj,i]/bXse[indj,i]))*bXse[indj,i]
Thetajj=LD*0
Thetajj[indj,indj]=Thetaj
bXestse0[,i]=sqrt(diag(Thetajj))*bXse[,i]
LDjj=LD%*%Thetajj%*%LD
bXest[,i]=as.vector(LD%*%(bXest0[,i]/bXse[,i]))*bXse[,i]
bXestse[,i]=sqrt(diag(LDjj))*bXse[,i]
}else{
bXest[,i]=bX[,i]
bXestse[,i]=bXse[,i]
}
}
}

return(list(bXest0=bXest0,bXest=bXest,bXestse0=bXestse0,bXestse=bXestse,xQTLfitList=xQTLfitList))
}
