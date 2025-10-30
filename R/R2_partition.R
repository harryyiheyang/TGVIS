#' R2_partition: Function to partition the Pratt index contributions of tissue-gene pairs, direct causal variants, and infinitesimal effects
#'
#' This function partitions the total Pratt index of all effects into contributions from tissue-gene pairs, direct causal variants, and infinitesimal effects, using  vectors of linear predictors and the outcome GWAS Z-scores.
#'
#' @param y A vector of Z-scores of marginal effects from the outcome GWAS.
#' @param eta.theta A vector of the linear predictor of tissue-gene causal effects.
#' @param eta.gamma A vector of the linear predictor of direct causal effects.
#' @param eta.upsilon A vector of the linear predictor of infinitesimal effects. Default is `NULL`.
#' @param LD A LD matrix.
#'
#' @return A data frame containing:
#' \item{r1}{Total Pratt index of all effects.}
#' \item{r2}{Pratt index contribution of tissue-gene pairs.}
#' \item{r3}{Pratt index contribution of direct causal variants.}
#' \item{r4}{Pratt index contribution of infinitesimal effects.}
#'
#' @export
#'
R2_partition <- function(y,eta.theta,eta.gamma,eta.upsilon=NULL,LD){
fitsqrt=matrixsqrt(LD)
TC=fitsqrt$wi
RC=fitsqrt$w
y=matrixVectorMultiply(TC,y)
eta.theta=matrixVectorMultiply(TC,eta.theta)
eta.gamma=matrixVectorMultiply(TC,eta.gamma)
if(is.null(eta.upsilon[1])==1){
A=R2est.adjust(pleiotropy=eta.gamma,y1=y,eta=eta.theta)
return(A)
}else{
eta.upsilon=matrixVectorMultiply(TC,eta.upsilon)
eta=eta.theta
y1=y
pleiotropy=eta.gamma
infeffect=eta.upsilon
Rsquare=data.frame(r1=0,r1se=0,r2=0,r2se=0,r3=0,r3se=0,r4=0,r4se=0)
if(sum(abs(infeffect))>0){
pleiotropy=vec(pleiotropy);eta=vec(eta);y1=vec(y1);infeffect=c(infeffect)
n=length(y1)
all=pleiotropy+eta+infeffect
theta.all=1
if(sum(abs(all))!=0){
theta.all=sd(all)/sd(y1)
all=all/sd(all)
}
sd.y=sd(y1)
y1=y1/sd.y
theta.eta=1
if(sum(abs(eta))!=0){
theta.eta=sd(eta)/sd.y
eta=eta/sd(eta)
}
theta.pleiotropy=1
if(sum(abs(pleiotropy))!=0){
theta.pleiotropy=sd(pleiotropy)/sd.y
pleiotropy=pleiotropy/sd(pleiotropy)
}
theta.infeffect=1
theta.infeffect=sd(infeffect)/sd.y
infeffect=infeffect/sd(infeffect)
beta.eta=cor(eta,y1);beta.eta.p=cor.test(eta,y1)$p.value;beta.eta.se=abs(beta.eta)/sqrt(qchisq(beta.eta.p,1,lower.tail=F))
beta.pleiotropy=cor(pleiotropy,y1);beta.pleiotropy.p=cor.test(pleiotropy,y1)$p.value;beta.pleiotropy.se=abs(beta.pleiotropy)/sqrt(qchisq(beta.pleiotropy.p,1,lower.tail=F))
beta.infeffect=cor(infeffect,y1);beta.infeffect.p=cor.test(infeffect,y1)$p.value;beta.infeffect.se=abs(beta.infeffect)/sqrt(qchisq(beta.infeffect.p,1,lower.tail=F))
beta.all=cor(all,y1);beta.all.p=cor.test(all,y1)$p.value;beta.all.se=abs(beta.all)/sqrt(qchisq(beta.all.p,1,lower.tail=F))
r1=theta.all*beta.all
r2=theta.eta*beta.eta
r3=theta.pleiotropy*beta.pleiotropy
r4=theta.infeffect*beta.infeffect
if(is.na(r1)==1) r1=0
if(is.na(r2)==1) r2=0
if(is.na(r3)==1) r3=0
if(is.na(r4)==1) r4=0
if(is.na(beta.eta)==1){
beta.eta=beta.eta.se=0
}
if(is.na(beta.pleiotropy)==1){
beta.pleiotropy=beta.pleiotropy.se=0
}
if(is.na(beta.infeffect)==1){
beta.infeffect=beta.infeffect.se=0
}
if(is.na(beta.all)==1){
beta.all=beta.all.se=0
}
r1se=sqrt(theta.all^2*beta.all.se^2+2*(1-r1)/length(y1)*r1)
r2se=sqrt(theta.eta^2*beta.eta.se^2+2*(1-r1)/length(y1)*r2)
r3se=sqrt(theta.pleiotropy^2*beta.pleiotropy.se^2+2*(1-r1)/length(y1)*r3)
r4se=sqrt(theta.infeffect^2*beta.infeffect.se^2+2*(1-r1)/length(y1)*r4)
if(r3==0){
r3se=0
}
Rsquare=data.frame(r1=r1,r2=r2,r3=r3,r4=r4)
}else{
Rsquare[1,1:4]=as.matrix(R2est.adjust(pleiotropy,y1,eta))
}
return(Rsquare)
}
}
