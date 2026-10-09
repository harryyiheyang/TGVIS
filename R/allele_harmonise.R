#' allele_harmonise: Function to harmonize allele coding between a reference panel and GWAS/eQTL summary data
#'
#' This function harmonizes the allele coding between an LD reference panel and GWAS/eQTL summary data to ensure consistency in the effect allele (A1) and the other allele (A2). The Z-scores in the GWAS/eQTL summary data are adjusted accordingly if the alleles are flipped.
#'
#' Matching uses a hashed \code{match()} on SNP IDs, so neither input is sorted
#' or modified: the caller's \code{ref_panel} and \code{gwas_data} keep their
#' row order and keys. Each reference SNP is matched to its first row in
#' \code{gwas_data}.
#'
#' @param ref_panel A data.frame or data.table representing an LD reference panel. It must contain columns "SNP", "A1" (effect allele), and "A2" (other allele).
#' @param gwas_data A data.frame or data.table representing GWAS/eQTL summary data. It must contain columns "SNP", "A1", "A2", and "Zscore".
#'
#' @return A data.table with one row per reference SNP found in \code{gwas_data}
#'   with matching or swapped alleles, in the row order of \code{ref_panel}.
#'   It holds the GWAS/eQTL columns with \code{Zscore} flipped for swapped
#'   alleles, then \code{A1}/\code{A2} from the reference panel, then any other
#'   reference-panel columns (prefixed \code{i.} if the name is also a
#'   GWAS/eQTL column).
#' @importFrom data.table as.data.table
#' @export
#'
allele_harmonise <- function(ref_panel, gwas_data) {
stopifnot(all(c("SNP", "A1", "A2") %in% colnames(ref_panel)))
stopifnot(all(c("SNP", "A1", "A2", "Zscore") %in% colnames(gwas_data)))

ref_panel <- as.data.frame(ref_panel)
gwas_data <- as.data.frame(gwas_data)

idx <- match(as.character(ref_panel$SNP), as.character(gwas_data$SNP))
ga1 <- as.character(gwas_data$A1)[idx]
ga2 <- as.character(gwas_data$A2)[idx]
ra1 <- as.character(ref_panel$A1)
ra2 <- as.character(ref_panel$A2)

same <- !is.na(idx) & ga1 == ra1 & ga2 == ra2
swap <- !is.na(idx) & ga1 == ra2 & ga2 == ra1
keep <- which(same | swap)

out <- gwas_data[idx[keep], setdiff(colnames(gwas_data), c("A1", "A2")), drop = FALSE]
out$Zscore <- ifelse(swap[keep], -out$Zscore, out$Zscore)

out$A1 <- ra1[keep]
out$A2 <- ra2[keep]
extra <- setdiff(colnames(ref_panel), c("SNP", "A1", "A2"))
if (length(extra)) {
  add <- ref_panel[keep, extra, drop = FALSE]
  clash <- colnames(add) %in% colnames(out)
  colnames(add)[clash] <- paste0("i.", colnames(add)[clash])
  out <- cbind(out, add)
}
rownames(out) <- NULL

as.data.table(out)
}
