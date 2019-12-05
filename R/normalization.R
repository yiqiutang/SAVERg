# normalization.R

#' Normalization
#'
#' Performs preprocessing on the raw count dataset.
#'
#' @references
#' * Standardize data columns in R.
#' Retrieved from https://stackoverflow.com/questions/15215457/standardize-data-columns-in-r.
#'
#' @param x The expression matrix that needs preprocessing or log-normalization.
#' @param percent Genes that are expressed in less than 100*percent\% of the cells are filtered out.
#' @return The preprocessed and log-normalized expression matrix.
#'
#' @examples
#' \dontrun{
#' preprocessing(ipsc_saver, percent = 0.1)
#' }
preprocessing <- function(x, percent){
  n <- dim(x)[2]
  gene.exprs.count <- rowSums(x != 0)
  x <- x[gene.exprs.count > n * percent, ]
  return(x)
}

#' Normalization
#'
#' Performs log-normalization on the raw count dataset.
#'
#' @param x The expression matrix that needs log-normalization.
#' @return The log-normalized expression matrix.
#'
#' @examples
#' \dontrun{
#' log_normalization(ipsc_saver)
#' }
log_normalization = function(x){
  sf <- colSums(x)/stats::median(colSums(x))
  return(log(sweep(x, 2, sf, '/')+1))
}

# [END]
