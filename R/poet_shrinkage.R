#' poet_shrinkage: Function to perform covariance matrix shrinkage using the POET method
#'
#' This function performs shrinkage of the covariance matrix using the POET (Principal Orthogonal complEment Thresholding) method, with an automatic selection of the number of latent factors and the optimal shrinkage parameter.
#'
#' @param LD A covariance matrix of linkage disequilibrium (LD).
#' @param KMax The maximum number of latent factors used in POET. Default is `min(15, round(nrow(LD)/2))`.
#' @param lamvec A vector of candidate shrinkage parameters. Default is `seq(0.025, 0.25, by = 0.025)`.
#' @param minvalue The threshold for the minimum eigenvalues used in determining the optimal shrinkage parameter. Default is `1e-3`.
#'
#' @return A covariance matrix that has been shrunk using the POET method.
#'
#' @importFrom CppMatrix matrixInverse matrixMultiply matrixVectorMultiply matrixEigen matrixListProduct
#' @export
#'
poet_shrinkage=function(LD,KMax=min(15,round(nrow(LD)/2)),lamvec=seq(0.025,0.25,by=0.025),minvalue=1e-3){
LD[is.na(LD)]=0;
diag(LD)=1
LD[abs(LD)<0.0001]=0
eig=matrixEigen(LD)
U=eig$vectors
d=eig$values
z=c()
for(j in 2:KMax){
z[j-1]=(d[j-1]-d[j])/(d[j]-d[j+1])
}
pck=which.max(z)+1
Uk=U[,1:pck];dk=d[1:pck]
hatc=(sum(diag(LD))-sum(dk))/(ncol(LD)-pck)
P=matrixMultiply(Uk,t(Uk)*dk)
E=LD-P;e=diag(E);e[e<0]=max(hatc,0.01);diag(E)=e
eigenvec=lamvec
for(i in 1:length(lamvec)){
E1=E*(1-lamvec[i])+diag(diag(E))*lamvec[i]
eigenvec=min(matrixEigen(E1)$values)
if(eigenvec>minvalue) break
}
hatLD=P+E1
return(cov2cor(hatLD))
}
