#' ctwas: Function to perform SuSiE-based causal tissue-gene pair and direct causal variant estimation
#'
#' This function performs the SuSiE (Sum of Single Effects) algorithm to estimate the causal tissue-gene pairs and direct causal variants using a provided set of summary statistics and linkage disequilibrium (LD) matrix. It returns a list of estimates and their associated statistics, including credible sets and Pratt estimations.
#'
#' @param by A vector of Z-scores of the marginal effects from the outcome GWAS.
#' @param bXest A matrix of direct effect estimates based on the Z-scores of tissue-gene pairs.
#' @param LD The LD matrix of variants.
#' @param Noutcome The sample size of the outcome GWAS.
#' @param L.causal The number of single effects used for tissue-gene pairs and direct causal variants. Default is 10.
#' @param pip.thres.cred The cumulative PIP threshold for variables in a credible set. Default is 0.95.
#' @param susie.iter The maximum number of iterations for `susie_rss` within the profile-likelihood algorithm. Default is 50.
#'
#' @return A list containing:
#' \item{theta}{The estimated effects for tissue-gene pairs, scaled by the outcome GWAS sample size.}
#' \item{gamma}{The estimated effects for direct causal variants, scaled by the outcome GWAS sample size.}
#' \item{theta.se}{Standard errors of the tissue-gene pair effects.}
#' \item{gamma.se}{Standard errors of the direct causal variant effects.}
#' \item{theta.pip}{Posterior inclusion probabilities (PIP) for tissue-gene pairs.}
#' \item{gamma.pip}{Posterior inclusion probabilities (PIP) for direct causal variants.}
#' \item{theta.pratt}{Pratt index for tissue-gene pairs.}
#' \item{gamma.pratt}{Pratt index for direct causal variants.}
#' \item{theta.cs}{Credible set indicators for tissue-gene pairs.}
#' \item{gamma.cs}{Credible set indicators for direct causal variants.}
#' \item{theta.cs.pip}{Cumulative PIP within credible sets for tissue-gene pairs.}
#' \item{gamma.cs.pip}{Cumulative PIP within credible sets for direct causal variants.}
#' \item{fit.causal}{The SuSiE fit object for the causal analysis.}
#' \item{cs.summary}{A summary of the credible sets obtained from the analysis.}
#' @importFrom CppMatrix matrixInverse matrixMultiply matrixVectorMultiply matrixEigen matrixListProduct matrixGeneralizedInverse
#' @importFrom susieR susie_rss susie_get_cs coef.susie
#' @importFrom Matrix Matrix solve
#' @export
#'
ctwas <- function(by, bXest, LD, Noutcome, L.causal = 10, pip.thres.cred = 0.95, susie.iter = 500) {
####################################### Preparing the data #####################################
n <- length(by)
p <- dim(bXest)[2]
Theta <- matrixInverse(LD)
XR <- cbind(matrixMultiply(LD, bXest), LD)
Xty <- c(t(bXest) %*% by, by)
XtX <- Matrix::bdiag(matrixMultiply(t(bXest) %*% LD, bXest), LD)
XtX <- as.matrix(XtX)
XtX[1:p, -c(1:p)] <- matrixMultiply(t(bXest), LD)
XtX[-c(1:p), c(1:p)] <- t(XtX[1:p, -c(1:p)])
dXtX <- diag(XtX)
dXtX[is.na(dXtX)] <- 1
dXtX[dXtX == 0] <- 1
XtX[is.na(XtX)] <- 0
diag(XtX) <- dXtX
Xadjust <- diag(XtX)
XtX <- cov2cor(XtX)
XtX=(t(XtX)+XtX)/2
XtyZ <- Xty / sqrt(Xadjust)
prior.weight.theta <- rep(1 / p, p)
prior.weight.gamma <- rep(1 / n, n)
prior_weights <- c(prior.weight.theta, prior.weight.gamma)
####################################### Performing SuSiE ##############################################################
fit.causal <- susie_rss(z = XtyZ, R = XtX, n = Noutcome, L = L.causal, residual_variance = 1, estimate_prior_method = "EM", prior_weights = prior_weights, intercept = FALSE, max_iter = susie.iter)
####################################### Resampling the SE ##############################################################
fit.causal.sampling <- susie.resampling(alpha = fit.causal$alpha, mu = fit.causal$mu, mu2 = fit.causal$mu2)
fit.causal$beta.se <- fit.causal.sampling$sd
####################################### Getting the variables in 95% credible sets ####################################
causal.cs <- group.pip.filter(pip.summary = summary(fit.causal)$var, pip.thres.cred = pip.thres.cred)
pip.alive <- causal.cs$ind.keep
pip.remove <- setdiff(1:ncol(XtX), pip.alive)
###################################### Preparing the results ################################
thetagamma <- coef.susie(fit.causal)[-1]
thetagamma[pip.remove] <- 0
thetagamma.se <- fit.causal$beta.se
thetagamma.se[pip.remove] <- 0
thetagamma.pip <- fit.causal$pip
thetagamma.pip[pip.remove] <- 0
thetagamma.cs <- causal.cs$cs
thetagamma.cs.pip <- causal.cs$cs.pip
thetagamma.cs.pip[pip.remove] <- 0
thetagamma.cs[pip.remove] <- 0
theta <- thetagamma[1:p]
gamma <- thetagamma[-c(1:p)]
theta.se <- thetagamma.se[1:p]
gamma.se <- thetagamma.se[-c(1:p)]
theta.pip <- thetagamma.pip[1:p]
gamma.pip <- thetagamma.pip[-c(1:p)]
theta.cs <- thetagamma.cs[1:p]
gamma.cs <- thetagamma.cs[-c(1:p)]
theta.cs.pip <- thetagamma.cs.pip[1:p]
gamma.cs.pip <- thetagamma.cs.pip[-c(1:p)]
gamma.pratt <- prattestimation(by = by, bXest = diag(length(by)), LD = LD, Theta = Theta, theta = gamma * sqrt(Noutcome))
theta.pratt <- prattestimation(by = by, bXest = bXest, LD = LD, Theta = Theta, theta = theta * sqrt(Noutcome) / sqrt(Xadjust[1:p]))
names(theta) <- names(theta.pip) <- colnames(bXest)
names(gamma) <- names(gamma.pip) <- rownames(bXest)
A <- list(
theta = theta * sqrt(Noutcome) / sqrt(Xadjust[1:p]),
gamma = gamma * sqrt(Noutcome),
theta.se = theta.se * sqrt(Noutcome) / sqrt(Xadjust[1:p]),
gamma.se = gamma.se * sqrt(Noutcome),
theta.pip = theta.pip,
gamma.pip = gamma.pip,
theta.pratt = theta.pratt,
gamma.pratt = gamma.pratt,
theta.cs = theta.cs,
gamma.cs = gamma.cs,
theta.cs.pip = theta.cs.pip,
gamma.cs.pip = gamma.cs.pip,
fit.causal = fit.causal,
cs.summary = causal.cs
)
return(A)
}
