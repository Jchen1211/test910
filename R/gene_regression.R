#' Construct data for gene regression
#'
#' @param count A Seurat object.Noted that the object has to be normalized and scaled by standard procedure.
#' @param indmatrix A perturbation-cell matrix calculated by indmatrix function.
#' @param selected Interested features list. Default: NULL.
#' @param neg A parameter for control cell. Non-perturbed cells are recommended to assign NT.
#' @param gene_frac A parameter for filtering low expression genes. Default: 0.01.
#' @param slot A parameter for extracting data from Seurat object. Default: 'scale.data'.
#'
#' @return A list for gene regression. Xmat represents cell identity matrix that contains cell and perturbation parameters. Ymat represents gene count matrix.
#' @export
#' @importFrom stats quantile
#' @examples
#' mat_reg <- gene_regression(count, ind, selected = NULL,neg,
#'                             gene_frac = 0.01, slot = 'scale.data')
gene_regression<-function(count, indmatrix, selected = NULL,neg, gene_frac = 0.01, slot = 'scale.data'){
  threshd = 0.95
  raw_count = GetAssayData(count, slot = 'counts')
  sc_count = GetAssayData(count, slot = slot)

  sel_genes = rownames(sc_count)
  sel_genes = sel_genes [!is.na(sel_genes)]
  if (is.null(selected) == FALSE) {
    sel_genes = sel_genes[sel_genes %in% selected]
    if (length(sel_genes) == 0) {
      stop("No genes left after selected")
    }
  }

  sel_genes = sel_genes[sel_genes %in% rownames(raw_count) ]
  raw_count2 = as.matrix(raw_count[sel_genes,])
  sel_genes = rownames(raw_count2)[which(rowSums(raw_count2 != 0) >= ncol(raw_count2) * gene_frac)]

  if (length(sel_genes) == 0) {
    stop("No genes left after filtered")
  }

  sel_cells = rownames(indmatrix)
  sel_cells = sel_cells[sel_cells %in% colnames(sc_count)]

  sel_cells = sel_cells[!is.na(sel_cells) & sel_cells %in% colnames(sc_count)]

  Ymat = as.matrix(sc_count[sel_genes, sel_cells])
  Ymat = t(Ymat)

  tgf = colnames(indmatrix)
  tgf[tgf %in% neg] = "NT"
  tgp = as.factor(tgf)

  Xmat = matrix(rep(0, length(sel_cells) * length(unique(tgp))), nrow = length(sel_cells))
  rownames(Xmat) = sel_cells
  colnames(Xmat) = levels(tgp)
  for (cnl in colnames(indmatrix)) {
    cellns = which(indmatrix[, cnl] == TRUE)
    if (cnl %in% neg) {
      Xmat[cellns, "NT"] = 1
    } else {
      Xmat[cellns, cnl] = 1
    }
  }
  Xmat[, "NT"] = 1

  Ymat_out = apply(Ymat, 2, function(X) {
    return(quantile(X, probs = threshd))
  })
  out_mat = t(matrix(rep(Ymat_out, nrow(Ymat)), ncol = nrow(Ymat)))
  Ymat_cor = ifelse(Ymat > out_mat, out_mat, Ymat)
  Ymat = Ymat_cor

  return(list(Xmat, Ymat))
}
