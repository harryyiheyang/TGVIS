vec=function(a){
as.vector(a)
}

#' @export
get_nonzero_rows <- function(mat, genes) {
if (!all(genes %in% colnames(mat))) {
stop("Some genes are not present in the matrix columns.")
}
selected_cols <- mat[, genes, drop = FALSE]
row_names <- apply(selected_cols, 2, function(x) {
nonzero_rows <- rownames(selected_cols)[x != 0]
if (length(nonzero_rows) > 0) {
paste(nonzero_rows, collapse = "~")
} else {
""
}
})
output_df <- data.frame(Genes = genes, NonzeroRows = row_names, stringsAsFactors = FALSE)
return(output_df)
}

colSD=function(A){
B=A
p=ncol(A)
for(i in 1:p){
B[,i]=sd(A[,i])
}
return(B)
}

soft=function(a,b){
c=abs(a)-b
c[c<0]=0
c=c*sign(a)
return(c)
}

positiveadj=function(A,min.eps=0.001){
a=matrixEigen(A)
d=c(a$values)
d[d<min.eps]=min.eps
B=matrixMultiply(a$vectors,t(a$vectors)*d)
return(B)
}

mcp=function(a,lam,ga=3.7){
b=abs(a)
z=soft(a,lam)/(1-1/ga)
z[which(b>(ga*lam))]=a[which(b>(ga*lam))]
return(z)
}

trace=function(A){
a=sum(diag(A))
return(a)
}

bimin=function(mat){
min_element <- min(mat)
min_indices <- which(mat == min_element, arr.ind = TRUE)
if (nrow(min_indices) > 1) {
min_indices <- min_indices[nrow(min_indices), ]
}
return(min_indices)
}

matrixsqrt=function(A){
fit=matrixEigen(A)
d=c(fit$value)
d1=d*0
d1[d>0]=1/d[d>0]
d=sqrt(d)
d1=sqrt(d1)
A=matrixMultiply(fit$vector,t(fit$vector)*d)
B=matrixMultiply(fit$vector,t(fit$vector)*d1)
C=list(w=A,wi=B,eigenfit=fit)
return(C)
}

R2est.adjust <- function(pleiotropy,y1,eta){
pleiotropy=vec(pleiotropy);eta=vec(eta);y1=vec(y1)
n=length(y1)
all=pleiotropy+eta
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
beta.eta=cor(eta,y1);beta.eta.p=cor.test(eta,y1)$p.value;beta.eta.se=abs(beta.eta)/sqrt(qchisq(beta.eta.p,1,lower.tail=F))
beta.pleiotropy=cor(pleiotropy,y1);beta.pleiotropy.p=cor.test(pleiotropy,y1)$p.value;beta.pleiotropy.se=abs(beta.pleiotropy)/sqrt(qchisq(beta.pleiotropy.p,1,lower.tail=F))
beta.all=cor(all,y1);beta.all.p=cor.test(all,y1)$p.value;beta.all.se=abs(beta.all)/sqrt(qchisq(beta.all.p,1,lower.tail=F))
r1=theta.all*beta.all
r2=theta.eta*beta.eta
r3=theta.pleiotropy*beta.pleiotropy
if(is.na(r1)==1) r1=0
if(is.na(r2)==1) r2=0
if(is.na(r3)==1) r3=0
if(is.na(beta.eta)==1){
beta.eta=beta.eta.se=0
}
if(is.na(beta.pleiotropy)==1){
beta.pleiotropy=beta.pleiotropy.se=0
}
if(is.na(beta.all)==1){
beta.all=beta.all.se=0
}
r1se=sqrt(theta.all^2*beta.all.se^2+2*(1-r1)/length(y1)*r1)
r2se=sqrt(theta.eta^2*beta.eta.se^2+2*(1-r1)/length(y1)*r2)
r3se=sqrt(theta.pleiotropy^2*beta.pleiotropy.se^2+2*(1-r1)/length(y1)*r3)
if(r3==0){
r3se=0
}
A=data.frame(r1=r1,r2=r2,r3=r3)
return(A)
}

susie.resampling=function(alpha,mu,mu2,sampling=500){
L=dim(mu)[1]
n=dim(mu)[2]
a=matrix(0,sampling,n)
for(i in 1:sampling){
N=mu*0
for(j in 1:L){
phi=c(rmultinom(1,1,alpha[j,]))
varr=c(mu2[j,])-c(mu[j,])^2
varr[varr<0]=0
varr=sqrt(varr)
N[j,]=(rnorm(n=n,mean=0,sd=1)*varr+c(mu[j,]))*phi
}
a[i,]=colSums(N)
}
return(list(mean=colMeans(a),sd=colSD(a),cov=cov(a)))
}

tgfm.resampling=function(alpha,mu,mu2,sampling=500){
L=dim(mu)[1]
n=dim(mu)[2]
a=c(mu[1,]*0)
M=matrix(0,sampling,n)
for(j in 1:L){
phi=t(rmultinom(sampling,1,alpha[j,]))
varr=c(mu2[j,])-c(mu[j,])^2
varr[varr<0]=0
N=MASS::mvrnorm(n=sampling,mu=mu[j,],Sigma=diag(varr))
N=N*phi
M=M+N
}
return(M)
}

susie_get_cs_index=function(fit){
cs=summary(fit)
s=cs$vars$variable[which(cs$vars$cs>0)]
return(s)
}

get_active_indices <- function(fit) {
cs = tryCatch(summary(fit), error = function(e) NULL)
if (!is.null(cs) && length(cs$cs) > 0) {
active_idx = unique(unlist(cs$cs$cs))
}else{
active_idx=NULL
}
return(active_idx)
}

group.pip.filter=function(pip.summary,pip.resamping=NULL,pip.thres.cred=0.5){
if(is.null(pip.resamping[1])==0){
pip.summary$variable_prob=pip.resamping[pip.summary$variable]
}
ind=which(pip.summary$cs>0)
if(length(ind)>0){
J=max(pip.summary$cs[ind])
pip.summary$cs.pip=pip.summary$variable_prob
for(i in 1:J){
indi=which(pip.summary$cs==i)
summaryi=pip.summary[indi,]
pip.cred=sum(summaryi$variable_prob)
pip.summary$cs.pip[indi]=pip.cred
}
ind.keep=which(pip.summary$cs.pip>=pip.thres.cred)
cs=pip.summary$cs
cs.pip=pip.summary$cs.pip
cs->cs[pip.summary$variable]
cs.pip->cs.pip[pip.summary$variable]
cs[which(cs==-1)]=0
}else{
ind.keep=NULL
cs=pip.summary$cs.pip*0
cs.pip=pip.summary$cs.pip*0
}
return(list(ind.keep=pip.summary$variable[ind.keep],cs=cs,cs.pip=cs.pip,result=pip.summary))
}

inf.test=function(res.inf,LD,LD2,Theta,A,var.res=1){
if(is.matrix(A)==T){
Varinf=LD
Var2=LD2
V=Theta/var.res
AT=matrixMultiply(t(A),V)
P=V-matrixMultiply(t(AT),matrixMultiply(matrixGeneralizedInverse(matrixMultiply(AT,A)),AT))
u=sum(res.inf^2)/2
PVar2=matrixMultiply(P,Var2)
e=sum(diag(PVar2))/2
h=sum(diag(matrixMultiply(PVar2,PVar2)))/2
kappa=h/2/e
v=2*e^2/h
}
if(is.vector(A)==T){
A=c(A)
Var2=LD2
V=Theta/var.res
delta=1/sum(A*matrixVectorMultiply(V,A))
AAT=matrixMultiply(t(t(A)),t(A))/delta
P=V-matrixListProduct(list(V,AAT,V))
u=sum(res.inf^2)/2/var.res
PVar2=matrixMultiply(P,Var2)
e=sum(diag(PVar2))/2
h=sum(diag(matrixMultiply(PVar2,PVar2)))/2
kappa=h/2/e
v=2*e^2/h
}
if(is.null(A)==1){
Varinf=LD
Var2=LD2
V=Theta/var.res
u=sum(res.inf^2)/2/var.res
VVar2=LD
e=sum(diag(VVar2))/2
h=sum(diag(LD2))/2
kappa=h/2/e
v=2*e^2/h
}
pv=pchisq(u/kappa,v,lower.tail=F)
return(pv)
}

prattestimation=function(by,bXest,LD,Theta,theta){
Thetay=matrixVectorMultiply(Theta,by)
vary=sum(by*Thetay)
varX=bXest[1,]
for(j in 1:length(varX)){
varX[j]=sum(bXest[,j]*matrixVectorMultiply(LD,bXest[,j]))
}
Xty=matrixVectorMultiply(t(bXest),by)
corXy=Xty/sqrt(varX)/sqrt(vary)
corXy[is.na(corXy)]=0
theta1=theta*sqrt(varX/vary)
pratt=corXy*theta1
return(pratt)
}

colSD=function(A){
p=ncol(A)
a=c(1:p)
for(i in 1:p){
a[i]=sd(A[,i])
}
return(a)
}

eigen_cumsum=function(D,thres=0.99){
d=cumsum(D)/sum(D)
ind=which(d>=thres)
if(length(ind)>1){
K=ind[1]
}
if(length(ind)==1){
K=length(D)
}
return(K)
}

#' @export
findUniqueNonZeroRows=function(bXest) {
M=bXest
nonZeroCounts <- colSums(M != 0)
uniqueCols <- which(nonZeroCounts == 1)
if(length(uniqueCols) == 0) {
uniqueRows=NULL
}
if(length(uniqueCols)==1){
uniqueRows=which(M[,uniqueCols]!=0)
}
if(length(uniqueCols)>1){
nonZeroCounts=rowSums(M[,uniqueCols]!=0)
uniqueRows <- unique(which(nonZeroCounts>0))
}
return(uniqueRows)
}

generate_block_matrix <- function(vars_df, s, theta) {
ind=which(theta==0)
vars_df$cs[which(vars_df$variable%in%ind)]=-1
concerned_vars <- vars_df[vars_df$cs != -1, ]
cs_values <- unique(concerned_vars$cs)
max_var_index <- max(vars_df$variable)
final_matrix <- matrix(0, nrow = max_var_index, ncol = max_var_index)
for (cs_val in cs_values) {
group_vars <- concerned_vars$variable[concerned_vars$cs == cs_val]
if (length(group_vars) > 1) {
group_s <- s[group_vars]
D <- generate_D_matrix(group_s,sign(theta[group_vars]))
final_matrix[group_vars, group_vars]=D
}
}
return(final_matrix)
}

generate_block_matrix_CARMA <- function(sumstat.result, s, theta) {
ind=which(theta==0)
concerned_vars <- sumstat.result[sumstat.result$cs >0, ]
cs_values <- unique(concerned_vars$cs)
max_var_index <- nrow(sumstat.result)
final_matrix <- matrix(0, nrow = max_var_index, ncol = max_var_index)
for (cs_val in cs_values) {
group_vars <- concerned_vars$variable[concerned_vars$cs == cs_val]
if (length(group_vars) > 1) {
group_s <- s[group_vars]
D <- generate_D_matrix(group_s,sign(theta[group_vars]))
final_matrix[group_vars, group_vars]=D
}
}
return(final_matrix)
}

top_K_pip=function(susie_summary,top_K=1,pip.min.thres=0.01,xQTL.pip.thres=0.5){
ind=which(susie_summary$cs>0&susie_summary$variable_prob>=pip.min.thres)
if(length(ind)>0){
susie_summary=susie_summary[ind,]
J=max(susie_summary$cs)
index=c()
for(j in 1:J){
indj=which(susie_summary$cs==j)
g=susie_summary[indj,]
if(length(indj)<=top_K){
index=c(index,g$variable)
}
if(length(indj)>top_K){
index=c(index,g$variable[top_K_indices(g$variable_prob,k=top_K)])
}
}
}
if(length(ind)==0){
index=which(susie_summary$variable_prob>=xQTL.pip.thres)
}
return(index)
}

group.pip.filter.xQTL=function(pip.summary,xQTL.cred.thres=0.95,xQTL.pip.thres=0.1){
ind=which(pip.summary$cs>0)
if(length(ind)>0){
J=max(pip.summary$cs[ind])
pip.summary$cs.pip=pip.summary$variable_prob
for(i in 1:J){
indi=which(pip.summary$cs==i)
summaryi=pip.summary[indi,]
pip.cred=sum(summaryi$variable_prob)
pip.summary$cs.pip[indi]=pip.cred
}
ind.keep=which(pip.summary$cs.pip>=xQTL.cred.thres&pip.summary$variable_prob>=xQTL.pip.thres)
cs=pip.summary$cs
cs.pip=pip.summary$cs.pip
cs->cs[pip.summary$variable]
cs.pip->cs.pip[pip.summary$variable]
cs[which(cs==-1)]=0
}else{
ind.keep=NULL
cs=pip.summary$cs.pip*0
cs.pip=pip.summary$cs.pip*0
}
return(list(ind.keep=pip.summary$variable[ind.keep],cs=cs,cs.pip=cs.pip,result=pip.summary))
}

generate_D_matrix <- function(s, sign_vec) {
p <- length(s)
if (length(sign_vec) != p) {
stop("Length of sign_vec must match length of s.")
}

if (p == 1) {
D <- 0
} else {
num_pairs <- p*(p-1)/2
D_all <- matrix(0, nrow = num_pairs, ncol = p)

row_idx <- 1
for (i in 1:(p-1)) {
for (j in (i+1):p) {
D_all[row_idx, i] <-  sign_vec[i] / s[i]
D_all[row_idx, j] <- -sign_vec[j] / s[j]
row_idx <- row_idx + 1
}
}

D <- t(D_all)%*%D_all
}

return(D)
}

top_K_indices <- function(vec, k=1) {
return(order(vec, decreasing = TRUE)[1:k])
}

find_cumvar_index <- function(eigenvalues, thres = 1) {
cumvar <- cumsum(eigenvalues) / sum(eigenvalues)
idx <- which(cumvar >= thres)[1]
if(is.na(idx)) {
warning("No index found where cumulative variance exceeds threshold")
return(length(eigenvalues))
}
return(idx)
}

prattestimation_theta=function(by,bXest,LD,Theta,theta){
bXest=as.matrix(bXest)
Thetay=matrixVectorMultiply(Theta,by)
vary=sum(by*Thetay)
varX=bXest[1,]
for(j in 1:length(varX)){
varX[j]=sum(bXest[,j]*matrixVectorMultiply(LD,bXest[,j]))
}
Xty=matrixVectorMultiply(t(bXest),by)
corXy=Xty/sqrt(varX)/sqrt(vary)
corXy[is.na(corXy)]=0
theta1=theta*sqrt(varX/vary)
pratt=corXy*theta1
return(pratt)
}

prattestimation_gamma=function(by,LD,Theta,gamma){
m=nrow(LD)
Thetay=matrixVectorMultiply(Theta,by)
vary=sum(by*Thetay)
varX=rep(1,m)
Xty=by
corXy=Xty/sqrt(varX)/sqrt(vary)
corXy[is.na(corXy)]=0
theta1=gamma*sqrt(varX/vary)
pratt=corXy*theta1
return(pratt)
}

#' @export
summary.tgvis <- function(fit, bXest = NULL, quiet = FALSE) {
nz_idx <- function(x) {
if (is.null(x)) return(integer(0))
x2 <- x
x2[is.na(x2)] <- 0
which(x2 > 0)
}
safe_get <- function(x, idx, prefix) {
if (is.null(x) || length(idx) == 0) return(setNames(numeric(0), NULL))
y <- x[idx]
nm <- names(x)
if (is.null(nm)) {
names(y) <- paste0(prefix, "_", idx)
} else {
names(y) <- nm[idx]
}
y
}
safe_vec <- function(x, len, fill = NA_real_) {
if (length(x) == 0) return(rep(fill, len))
if (length(x) == len) return(x)
if (length(x) == 1) return(rep(x, len))
y <- rep(fill, len)
y[seq_len(min(len, length(x)))] <- x[seq_len(min(len, length(x)))]
y
}
order_na_last <- function(...) {
mats <- list(...)
keys <- Map(function(v) list(is.na(v), v), mats)
do.call(order, c(unlist(keys, recursive = FALSE), list(na.last = TRUE)))
}

idx_theta <- nz_idx(fit$theta.cs)
theta    <- safe_get(fit$theta,       idx_theta, "Gene")
th_pip   <- safe_get(fit$theta.pip,   idx_theta, "Gene")
th_pratt <- safe_get(fit$theta.pratt, idx_theta, "Gene")
th_cs    <- safe_get(fit$theta.cs,    idx_theta, "Gene")
th_cspip <- safe_get(fit$theta.cs.pip,idx_theta, "Gene")

df_theta <- if (length(theta)) {
data.frame(
Variable = names(theta),
estimate = as.numeric(theta),
pip      = safe_vec(as.numeric(th_pip),   length(theta)),
pratt    = safe_vec(as.numeric(th_pratt), length(theta)),
cs       = safe_vec(as.numeric(th_cs),    length(theta)),
cs.pip   = safe_vec(as.numeric(th_cspip), length(theta)),
Type     = rep("Gene", length(theta)),
stringsAsFactors = FALSE
)
} else {
NULL
}

idx_gamma <- nz_idx(fit$gamma.cs)
gamma    <- safe_get(fit$gamma,       idx_gamma, "SNP")
ga_pip   <- safe_get(fit$gamma.pip,   idx_gamma, "SNP")
ga_pratt <- safe_get(fit$gamma.pratt, idx_gamma, "SNP")
ga_cs    <- safe_get(fit$gamma.cs,    idx_gamma, "SNP")
ga_cspip <- safe_get(fit$gamma.cs.pip,idx_gamma, "SNP")
df_gamma <- if (length(gamma)) {
data.frame(
Variable = names(gamma),
estimate = as.numeric(gamma),
pip      = safe_vec(as.numeric(ga_pip),   length(gamma)),
pratt    = safe_vec(as.numeric(ga_pratt), length(gamma)),
cs       = safe_vec(as.numeric(ga_cs),    length(gamma)),
cs.pip   = safe_vec(as.numeric(ga_cspip), length(gamma)),
Type     = rep("SNP", length(gamma)),
stringsAsFactors = FALSE
)
} else {
NULL
}
A <- if (!is.null(df_theta) || !is.null(df_gamma)) {
A <- rbind(df_theta, df_gamma)
rownames(A) <- NULL
A
} else {
if (!quiet) message("No direct causal variant or gene–tissue pair identified.")
return(data.frame(
Variable = character(0),
cs       = numeric(0),
cs.pip   = numeric(0),
cs.pratt = numeric(0),
xQTL     = character(0),
Type     = character(0),
estimate = numeric(0),
pip      = numeric(0),
pratt    = numeric(0),
stringsAsFactors = FALSE
))
}
if (!("cs" %in% names(A)))   A$cs   <- NA_real_
if (!("pratt" %in% names(A))) A$pratt <- NA_real_
cs_display <- A$cs
cs_group <- A$cs
cs_group[is.na(cs_group)] <- -Inf

A$cs.pratt <- ave(A$pratt, cs_group, FUN = function(v) sum(v, na.rm = TRUE))
A$cs <- cs_display
A$xQTL <- A$Variable
if (any(A$Type == "Gene", na.rm = TRUE) &&
!is.null(bXest) &&
exists("get_nonzero_rows", mode = "function")) {

exp_idx <- which(A$Type == "Gene")
mapped  <- tryCatch(
get_nonzero_rows(bXest, A$Variable[exp_idx]),
error = function(e) NULL
)
if (!is.null(mapped) && "NonzeroRows" %in% names(mapped)) {
A$xQTL[exp_idx] <- mapped$NonzeroRows
}
}
need_cols <- c("Variable","cs","cs.pip","cs.pratt","xQTL",
           "Type","estimate","pip","pratt")
for (cc in need_cols) if (!cc %in% names(A)) A[[cc]] <- NA
A <- A[, need_cols, drop = FALSE]
ord <- order_na_last(A$cs, A$cs.pratt, A$Type, A$Variable)
A <- A[ord, , drop = FALSE]
rownames(A) <- NULL
colnames(A)=tolower(colnames(A))
return(A)
}
