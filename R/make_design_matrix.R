#' make_design_matrix: Constructs a Design Matrix of Z-scores for Gene-Tissue Pairs
#'
#' This function constructs a design matrix from a `data.table` where each row represents a SNP, each column represents a tissue-gene pair, and the values are z-scores. If a tissue-gene pair does not have a z-score for a SNP, the corresponding entry will be `NA`.
#'
#' @param df A `data.table` with three columns: `SNP`, `Variable` (representing tissue-gene pairs), and `Zscore`.
#'
#' @importFrom tidyr spread
#'
#' @return A matrix where rows correspond to SNPs, columns correspond to tissue-gene pairs, and values are z-scores. Missing values (`NA`) indicate that the tissue-gene pair does not have a z-score for that SNP.
#'
#' @export
#'
make_design_matrix <- function(df) {
  names(df) = c("SNP", "Variable", "Zscore")
  df = df[!duplicated(df[, c("SNP", "Variable")]), ]
  Z = tidyr::spread(df, Variable, Zscore)
  rownames(Z) = Z[, 1]
  Z = Z[, -1]
  return(Z)
}
