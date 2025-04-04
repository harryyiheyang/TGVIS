#' modified_predixcan: Function to estimate causal effects using a modified PrediXcan approach with BIC-based tuning parameter selection
#'
#' This function implements a modified version of the PrediXcan method, using an alternating direction method of multipliers (ADMM) algorithm and BIC for selecting the optimal tuning parameter. It estimates the causal effects of tissue-gene pairs while accounting for direct causal variants.
#'
#' @param by A vector of Z-scores of marginal effects from the outcome GWAS (same as in ctwas).
#' @param bxest A vector of direct xQTL effect estimates for a tissue-gene pair.
#' @param LD The LD matrix of variants (same as in ctwas).
#' @param tauvec A vector of tuning parameters used in penalizing the direct causal effect. Default is `seq(3,10,by=1)`.
#' @param rho.gamma A parameter set in the ADMM algorithm. Default is 1.5.
#' @param max.iter The maximum number of iterations for the ADMM algorithm. Default is 15.
#' @param max.eps The convergence tolerance for the ADMM algorithm. Default is 0.005.
#' @param ebic.factor The extended BIC factor for model selection. Default is 2.
#' @param pleiotropy.rm A vector of indices specifying which variants should not be considered as having direct causal effects.
#'
#' @return A list containing:
#' \item{theta}{The estimated effect size of the tissue-gene pair.}
#' \item{gamma}{The estimated effect sizes of the direct causal variants.}
#' \item{covtheta}{The covariance of the estimated effect size `theta`.}
#' \item{Bic}{The BIC values for each tuning parameter.}
#' \item{Btheta}{The estimated `theta` values for each tuning parameter.}
#' \item{Bgamma}{The estimated `gamma` values for each tuning parameter.}
#' \item{Eta}{The estimated linear predictor for the tissue-gene pair.}
#'
#' @importFrom CppMatrix matrixInverse matrixMultiply matrixVectorMultiply matrixEigen matrixListProduct
#' @export
#'
modified_predixcan=function(by,bxest,LD,tauvec=seq(3,10,by=1),rho.gamma=1.5,max.iter=15,max.eps=0.005,ebic.factor=2,normmax=2,pleiotropy.rm=NULL){
n=length(by)
pleiotropy.keep=setdiff(c(1:n),pleiotropy.rm)
dx=matrixVectorMultiply(LD,bxest)
Theta=matrixInverse(LD)
yinv=c(matrixVectorMultiply(Theta,by))
xtx=sum(bxest*dx)
xty=sum(by*bxest)
theta.ini=xty/xtx
Thetarho=matrixInverse(LD[pleiotropy.keep,pleiotropy.keep]+rho.gamma*diag(length(pleiotropy.keep)))
tauvec=sort(tauvec,decreasing=F)
w=length(tauvec)
Btheta=c(1:w)
Bgamma=matrix(0,n,w)
Bbic=c(1:w)
for(sss in c(w:1)){
theta=theta.ini
gamma=by*0
gamma1=gamma
u=rho.gamma*(gamma-gamma1)
theta1=theta*0
error=1
iter=1
while(error>max.eps&iter<max.iter){
theta1=theta
theta=(xty-sum(dx*gamma))/xtx
res=c(by-dx*theta-u+rho.gamma*gamma1)
gamma[pleiotropy.keep]=c(matrixVectorMultiply(Thetarho,res[pleiotropy.keep]))
gamma1=mcp(gamma+u/rho.gamma,tauvec[sss],ga=3)
u=u+rho.gamma*(gamma-gamma1)
u[pleiotropy.rm]=0
iter=iter+1
if(iter>3){
error=abs(theta-theta1)
}
Btheta[sss]=theta
Bgamma[,sss]=gamma1
r=vec(by-dx*theta-matrixVectorMultiply(LD,gamma1))
df=sum(gamma1!=0)
rss=sum(r*(matrixVectorMultiply(Theta,r)))
Bbic[sss]=n*log(rss)+log(n)*(1+ebic.factor)*df
}
}

star=which.min(Bbic)
theta=Btheta[star]
gamma=Bgamma[,star]
eta=bxest*theta

indgamma=which(gamma!=0)
effn=n-length(indgamma)
if(sum(indgamma)>0){
Z=cbind(bxest,diag(n)[,indgamma])
Hinv=matrixMultiply(t(Z),matrixMultiply(LD,Z))
Hinv=MASS::ginv(Hinv)
r=vec(yinv-bxest*theta-gamma)
varr=mean(r*matrixVectorMultiply(LD,r))*n/(effn-1)
covg=varr*Hinv
covg=covg[1,1]
covtheta=covg
}
if(sum(indgamma)==0){
r=vec(yinv-bxest*theta)
varr=mean(r*matrixVectorMultiply(LD,r))
covg=varr/xtx*n/(n-1)
covtheta=covg
}

A=list()
A$theta=theta
A$gamma=gamma
A$covtheta=as.numeric(covtheta)
A$Bic=Bbic
A$Btheta=Btheta
A$Bgamma=Bgamma
A$Eta=eta
return(A)
}
