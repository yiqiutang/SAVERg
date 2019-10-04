#' Normalization
#'
#' Performs preprocessing and log-normalization on the raw count dataset.
#'
#' @param x The expression matrix that needs preprocessing or log-normalization.
#' @param percent Genes that are expressed in less than 100*percent\% of the cells are filtered out.
#'
#' @return The preprocessed and log-normalized expression matrix.
#'


#' @rdname normalization
#' @export
preprocessing <- function(x, percent){
  n <- dim(x)[2]
  gene.exprs.count <- rowSums(x != 0)
  x = x[gene.exprs.count > n * percent, ]
  return(x)
}

#' @rdname normalization
#' @export
log_normalization = function(x){
  sf = colSums(x)/median(colSums(x))
  return(log(sweep(x, 2, sf, '/')+1))
}
