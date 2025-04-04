#' eQTLmapping_susie: Function to perform eQTL fine-mapping using SuSiE and optional resampling
#'
#' This function performs fine-mapping of eQTLs using the SuSiE (Sum of Single Effects) algorithm. It allows for optional resampling of eQTL effects and returns both the estimated effects and resampled effects.
#'
#' @param bX A matrix of Z-scores of marginal eQTL effect estimates for tissue-gene pairs.
#' @param LD The LD matrix of variants.
#' @param Nvec A vector representing the sample sizes of tissue-gene pair eQTL studies.
#' @param pip.thres A threshold for individual PIP when no credible set is found. Default is 0.2.
#' @param pip.min The minimum individual PIP in each 95\% credible set. Used to remove variables with low PIPs within credible sets. Default is 0.05.
#' @param L The number of single effects to be used in the SuSiE model. Default is 5.
#' @param resample A logical value indicating whether to resample eQTL effects. Default is `FALSE`.
#' @param sampling.time The number of resampling iterations to perform when `resample` is `TRUE`. Default is 100.
#'
#' @return A list containing:
#' \item{Estimate}{A matrix of estimated eQTL effects for tissue-gene pairs.}
#' \item{Sampling}{A matrix of resampled eQTL effects, if resampling is enabled. Otherwise, this will contain zeros.}
#'
#' @importFrom CppMatrix matrixInverse matrixMultiply matrixVectorMultiply matrixEigen matrixListProduct
#' @importFrom susieR susie_rss susie_get_cs
#' @importFrom Matrix Matrix solve
#' @export
#'
eQTLmapping_susie <- function(bX,LD,Nvec,pip.thres=0.5,pip.min=0.2,L=5,resample=F,sampling.time=100) {
p <- ncol(bX)
B <- bX * 0
colnames(B)=colnames(bX)
rownames(B)=rownames(bX)
C=B
pb <- txtProgressBar(min = 0, max = p, style = 3)
for (i in 1:p) {
indx <- which(bX[, i] != 0)
a <- LD[indx, indx]
y <- bX[indx, i]
errorindicate <- 0
tryCatch({
if (length(indx) > 5) {
fit <- susie_rss(z = y, R = a, n = Nvec[i], L = L, verbose = FALSE)
index.causal = intersect(susie_get_cs_index(fit),which(fit$pip>pip.min))
z <- coef(fit)[-1] * sqrt(Nvec[i])
if(length(index.causal)>0){
z[-index.causal]=0
}else{
z=z*(fit$pip>pip.thres)
}
if(sum(abs(z))>0){
if(length(index.causal)>0&resample==T){
z1 <- tgfm.resampling(alpha = fit$alpha, mu = fit$mu, mu2 = fit$mu2, sampling = sampling.time) * sqrt(Nvec[i])
z1 = colMeans(z1)
z1[-index.causal]=0
}else{
z1=z*0
}
}else{
z1=z
}
} else {
z <- matrixVectorMultiply(matrixInverse(a), y)
z1 = z
}
B[indx, i] <- z
C[indx, i] <- z1
}, error = function(e) {
message(paste("Error in iteration", i))
errorindicate <- 1
})
if(errorindicate == 1){
next
}
setTxtProgressBar(pb, i)
}
close(pb)
rownames(B) = rownames(C) = rownames(bX)
colnames(B) = colnames(C) = colnames(bX)
return(list(Estimate=B,Sampling=C))
}
