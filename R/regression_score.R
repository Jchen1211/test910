#' Calculate perturbation effect
#' @description
#' Calculate perturbation effect by ridge regression in single-cell CRISPR screening data.
#'
#' @param count A Seurat object. Noted that the object has to be normalized and scaled by standard procedure.
#' @param sgrna A perturbation-cell data frame. Noted that the data frame contain perturbation and cell parameters.
#' @param negctr A parameter for non-perturbed cells which is recommended as NT.
#' @param selected_gene Interested features list. Default: NULL.
#' @param permut The number of permutations for p value. Default: NULL.
#' @param paramt The Regularization coefficient. Default: 0.01.
#' @param filt A parameter for filtering low expression genes. Default: 0.01.
#' @param slot A parameter for extracting data from Seurat object. Default: 'scale.data'.
#'
#' @return A list for regression result.
#' @export
#'
#' @examples
#' reg_score <- regression_score(count,sgrna,negctr,selected_gene = NULL,
#'                                permut = 1000,paramt = 0.01,filt = 0.01,slot = 'scale.data')
regression_score <- function(count,sgrna,negctr,selected_gene = NULL,permut = NULL,paramt = 0.01,filt = 0.01,slot = 'scale.data'){
  if (!is.null(permut)) {
    permutation = as.integer(permut)
  } else {
    permutation = 1000
  }

  sgr_count<-sgrna
  if (sum(colnames(sgr_count) %in% c("cell", "perturbation")) != 2) {
    stop("cell, perturbation column names not found in sgrna file.")
  }

  guide = table(sgr_count[,which(colnames(sgr_count)=="cell")])
  ncnt = table(table(sgr_count[,which(colnames(sgr_count)=="cell")]))
  tar_count = count

  mat = sum(sgr_count[, which(colnames(sgr_count)=="cell")] %in% colnames(tar_count))
  if (mat == 0) {
    stop("cell names in expression matrix and sgrna matrix do not match")
  }

  ind <- indmatrix(tar_count, sgr_count)#

  mat_reg = gene_regression(count = tar_count, indmatrix = ind, selected = selected_gene, neg = negctr, gene_frac = filt, slot = slot)

  Xmat = mat_reg[[1]]
  Ymat = mat_reg[[2]]

  perscore = matrix_permutation(Xmat, Ymat, lamb = paramt, permut = permutation)
  scr = perscore[[1]]
  scr_pval = perscore[[2]]

  return(list(score=data.frame(Perturbedgene = rownames(scr), scr),
              pvalue=data.frame(Perturbedgene = rownames(scr), scr_pval)))

}
