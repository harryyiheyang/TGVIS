#' eQTLmapping_susie: Function to perform eQTL fine-mapping using SuSiE and optional resampling
#'
#' This function performs fine-mapping of eQTLs using the SuSiE (Sum of Single Effects) algorithm. It allows for optional resampling of eQTL effects and returns both the estimated effects and resampled effects.
#'
#' @param bX A matrix of Z-scores of marginal eQTL effect estimates for tissue-gene pairs.
#' @param LD The LD matrix of variants.
#' @param Nvec A vector representing the sample sizes of tissue-gene pair eQTL studies.
#' @param pip_thres A threshold for individual PIP when no credible set is found. Default is 0.2.
#' @param pip.min The minimum individual PIP in each 95\% credible set. Used to remove variables with low PIPs within credible sets. Default is 0.05.
#' @param L The number of single effects to be used in the SuSiE model. Default is 8.
#' @param max_iter The maximum iterations in the SuSiE model. Default is 300.
#' @param coverage The coverage of credible set to be used in SuSiE. Default is 0.95.
#'
#' @return A matrix of estimated eQTL effects for tissue-gene pairs.
#'
#' @importFrom CppMatrix matrixInverse matrixMultiply matrixVectorMultiply matrixEigen matrixListProduct
#' @importFrom susieR susie_rss susie_get_cs
#' @importFrom Matrix Matrix solve
#' @export
#'
eQTLmapping_susie <- function(bX,LD,Nvec,pip_thres=0.5,pip.min=0.2,L=8,coverage=0.95,max_iter=300,...) {
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
if (length(indx) > 3) {
fit <- susie_rss(z = y, R = a, n = Nvec[i], L = L, verbose = FALSE, coverage=coverage, max_iter=max_iter,...)
fit <- susie_rss(z = y, R = a, n = Nvec[i], L = max(1,length(get_active_indices(fit))), max_iter=max_iter, verbose = FALSE, coverage=coverage,...)
index.causal = intersect(susie_get_cs_index(fit),which(fit$pip>pip.min))
z <- coef(fit)[-1] * sqrt(Nvec[i])
if(length(index.causal)>0){
z[-index.causal]=0
}else{
z=z*(fit$pip>pip_thres)
}
} else {
z <- matrixVectorMultiply(matrixInverse(a), y)
}
B[indx, i] <- z
}, error = function(e) {
message(paste0("Error in iteration ", i, ": ", conditionMessage(e)))
errorindicate <- 1
})
if(errorindicate == 1){
next
}
setTxtProgressBar(pb, i)
}
close(pb)
rownames(B) =  rownames(bX)
colnames(B) =  colnames(bX)
return(B)
}
