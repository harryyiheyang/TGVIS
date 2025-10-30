#' Build Blockwise Linkage Disequilibrium (LD) Matrix
#'
#' This function constructs a genome-wide LD matrix by combining block-level
#' eigen-decomposed LD components stored in a specified directory. Each LD block
#' is reconstructed using eigenvectors (`U`) and eigenvalues (`lambda`) from
#' `SBayesRC::readEig`, and re-ordered to match the SNP order in `GWAS_Locus`.
#'
#' @param GWAS_Locus A data frame containing at least two columns:
#'   \describe{
#'     \item{Block}{Character or integer ID indicating LD block assignment for each SNP.}
#'     \item{SNP}{SNP identifier, must be unique.}
#'   }
#'
#' @param ldDir A character string specifying the directory path that contains
#'   LD block eigen-decomposition files (compatible with `SBayesRC::readEig`).
#'
#' @param snpinfo A data frame with columns:
#'   \describe{
#'     \item{SNP}{SNP identifier.}
#'     \item{Block}{Block ID corresponding to each SNP.}
#'   }
#'   Used to match and order SNPs correctly within each block.
#'
#' @return A symmetric numeric matrix of dimension \eqn{n_{SNP} \times n_{SNP}},
#'   where each block on the diagonal corresponds to an LD submatrix for that block.
#' @importFrom SBayesRC readEig
#' @import dplyr
#' @importFrom CppMatrix matrixMultiply
#' @export
#'
build_LD_matrix <- function(GWAS_Locus, ldDir, snpinfo) {
  if(!all(c("Block", "SNP") %in% colnames(GWAS_Locus))) {
    stop("GWAS_Locus must contain 'Block' and 'SNP' columns")
  }
  if(!all(c("SNP", "Block") %in% colnames(snpinfo))) {
    stop("snpinfo must contain 'SNP' and 'Block' columns")
  }
  blocks <- unique(GWAS_Locus$Block)
  n_snps <- nrow(GWAS_Locus)
  R <- matrix(0, n_snps, n_snps)
  snp_names <- character(n_snps)
  start_idx <- 1
  for(i in seq_along(blocks)) {
    block <- blocks[i]
    block_snps <- GWAS_Locus$SNP[GWAS_Locus$Block == block]
    n_block <- length(block_snps)
    A <- SBayesRC::readEig(ldDir = ldDir, block = block)
    Ri <- CppMatrix::matrixMultiply(A$U, t(A$U) * A$lambda)
    rownames(Ri) <- colnames(Ri) <- snpinfo$SNP[snpinfo$Block == block]
    Ri <- Ri[block_snps, block_snps]
    end_idx <- start_idx + n_block - 1
    R[start_idx:end_idx, start_idx:end_idx] <- Ri
    snp_names[start_idx:end_idx] <- block_snps
    start_idx <- end_idx + 1
  }
  rownames(R) <- colnames(R) <- snp_names
  if(!all(rownames(R) == GWAS_Locus$SNP)) {
    warning("Row names of R do not match the order of SNPs in GWAS_Locus")
  }

  return(R)
}
