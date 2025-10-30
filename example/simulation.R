rm(list=ls())
library(MendelianRandomization)
library(MRBEEX)
devtools::document()
ARcov=function(p,rho){
  s=c(1:p)
  for(i in 1:p){
    s[i]=rho^(i-1)
  }
  return(toeplitz(s))
}
CScov=function(p,rho){
  return(matrix(rho,p,p)+(1-rho)*diag(p))
}
m=200 # number of IVs
p=10 # number of exposure
n1=3000 # exposure sample size
n0=5e5 # outcome sample size
pca.thres=0.999 # PCGMM threshold
pip.thres=0.3 # PIP threshold of SuSiE in cisMRBEE
ridge=0
Rbb=ARcov(p,-0.5) # exposure covariance
Ruv=ARcov(p+1,-0.3) # estimation error covariance
Nxy=c(rep(n1,p),n0) # sample size vector
Hxy=c(rep(.3,p),.001) # H2 vector for expsoures and outcome
Rnn=matrix(1,p+1,p+1) # sample overlap matrix
#We assume the cohorts of exposures are the subset of the cohort of outcome
theta0=c(1,0,-0.5,rep(0,p-3)) # true causal effect vector
UHP.var=1 # var(exposure)/var(UHP)=UHP.var
UHP.frac=0.01 # generate m x UHP.frac number of UHP
xQTL.pip.thres=0.5 # the PIP threshold of xQTL effect estimate
LD0=readRDS("example/LD.rds")%>%do.call(Matrix::bdiag,.)%>%as.matrix(.)
ldindex=sample(floor(2083/m),1,replace=F)
ldindex=c(((ldindex-1)*m+1):(ldindex*m))
LD1=LD0;diag(LD1)=0;
print(c(max(LD1),min(LD1)))
# The maximum correlation is 0.9
LD=LD0[ldindex,ldindex]*0.75+0.25*diag(m)
Theta=solve(LD)
# As the correlation may be very high, we consider a linear shrinkage
# of the true LD matrix
A=MRBEEX::summary_generation(theta=theta0,m=m,Rbb=Rbb,Ruv=Ruv,Rnn=Rnn,LD=LD,
                             Nxy=Nxy,non.zero.frac=rep(0.015,p),UHP.frac=UHP.frac,
                             CHP.frac=0,UHP.var=UHP.var,Hxy=Hxy,UHP.dis="normal")
bX=A$bX
by=A$by
var_inf=1e-6
upsilon=as.vector(LD%*%rnorm(n=m,0,sd=sqrt(var_inf)))
by=A$by+upsilon
if(sum(is.na(by))>0){
  next
}
bXse=A$bXse
byse=A$byse
Rxy=Rnn*Ruv
multiplefactor=colMeans(bXse)/mean(byse)
fit.cismrbee=MRBEEX::CisMRBEEX(by=by,bX=bX,byse=byse,bXse=bXse,xQTL.max.L=5,
                               xQTL.pip.thres=xQTL.pip.thres,xQTL.Nvec=rep(n1,p),LD=LD,
                               Rxy=Rnn*Ruv,reliability.thres=0.9,tauvec=seq(2,15,1),
                               ebic.gamma=2,ebic.theta=0,ridge=ridge,
                               causal.pip.thres=pip.thres)

fittgvis=tgvis(by=by/byse,bXest=fit.cismrbee$bXest0/bXse,LD=LD,Noutcome=n0,
                      L_vec=c(1:4),eigen_thres=1,pv_thres=5e-3,
                      ebic_beta=0,susie_iter=500)

fittgvis_new=tgvis(by=by/byse,bXest=fit.cismrbee$bXest0/bXse,LD=LD,L_vec=c(1:6),Noutcome=n0,var_inf=var_inf,estimate_inf=F)
fittgvis$theta/multiplefactor
fittgvis_new$theta/multiplefactor
fit.cismrbee$theta
