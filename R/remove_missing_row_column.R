#' remove_missing_row_column: Removes Rows and Columns with Excessive Missing Values from a Matrix.
#'
#' This function processes a matrix, typically the output from `make_design_matrix`, by removing rows and columns that have a high percentage of missing values (`NA`). The function allows for flexibility in the order of operations, either removing rows first or columns first.
#'
#' @param M A matrix, typically the output from `make_design_matrix`, where rows represent SNPs, columns represent tissue-gene pairs, and values are z-scores.
#' @param rowthres A numeric threshold (between 0 and 1). Rows with a proportion of missing values greater than this threshold will be removed.
#' @param colthres A numeric threshold (between 0 and 1). Columns with a proportion of missing values greater than this threshold will be removed.
#' @param rowfirst A logical value. If `TRUE`, rows are processed before columns. If `FALSE`, columns are processed before rows.
#'
#' @return A matrix with rows and columns containing excessive missing values removed based on the specified thresholds.
#'
#' @export
#'
remove_missing_row_column <- function(M, rowthres = 0.95, colthres = 0.95, rowfirst = TRUE) {
  M1 = is.na(M)
  if (rowfirst == TRUE) {
    M2 = rowMeans(M1)
    ind1 = which(M2 <= rowthres)
    M2 = colMeans(M1[ind1, ])
    ind2 = which(M2 <= colthres)
  } else {
    M2 = colMeans(M1)
    ind2 = which(M2 <= colthres)
    M2 = rowMeans(M1[, ind2])
    ind1 = which(M2 <= rowthres)
  }
  M3 = M[ind1, ind2]
  return(M3)
}
