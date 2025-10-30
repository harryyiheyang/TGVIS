#' allele_harmonise: Function to harmonize allele coding between a reference panel and GWAS/eQTL summary data
#'
#' This function harmonizes the allele coding between an LD reference panel and GWAS/eQTL summary data to ensure consistency in the effect allele (A1) and the other allele (A2). The Z-scores in the GWAS/eQTL summary data are adjusted accordingly if the alleles are flipped.
#'
#' @param ref_panel A data.frame or data.table representing an LD reference panel. It must contain columns "SNP", "A1" (effect allele), and "A2" (other allele).
#' @param gwas_data A data.frame or data.table representing GWAS/eQTL summary data. It must contain columns "SNP", "A1", "A2", and "Zscore".
#'
#' @return A data.table with harmonized alleles and adjusted Z-scores. The output contains the harmonized alleles from the reference panel and the original alleles from the GWAS/eQTL data as separate columns.
#' @importFrom data.table setDT setkey setnames
#' @export
#'
allele_harmonise <- function(ref_panel, gwas_data) {
stopifnot(all(c("SNP", "A1", "A2") %in% colnames(ref_panel)))
stopifnot(all(c("SNP", "A1", "A2", "Zscore") %in% colnames(gwas_data)))

setDT(ref_panel)
setDT(gwas_data)
setkey(ref_panel, SNP)
setkey(gwas_data, SNP)

merged_data <- gwas_data[ref_panel, nomatch = 0]
setnames(merged_data, "A1", "A1.x")
setnames(merged_data, "A2", "A2.x")
setnames(merged_data, "i.A1", "A1.y")
setnames(merged_data, "i.A2", "A2.y")

ind <- (merged_data$A1.x == merged_data$A1.y & merged_data$A2.x == merged_data$A2.y) |
  (merged_data$A1.x == merged_data$A2.y & merged_data$A2.x == merged_data$A1.y)

merged_data <- merged_data[ind]

merged_data[, Zscore := ifelse(A1.x == A2.y & A2.x == A1.y, -Zscore, Zscore)]

setnames(merged_data, "A1.y", "A1")
setnames(merged_data, "A2.y", "A2")
setnames(merged_data, "A1.x", "A11")
setnames(merged_data, "A2.x", "A21")
merged_data$A11=merged_data$A21=NULL

return(merged_data)
}
