# trajectory_analysis.R

#' Pseudotime trajectory analysis
#'
#' Carries out cell trajectory analysis to assess the imputation accuracy of SAVER.
#'
#' @param count_data The count expression matrix. The rows correspond to genes and
#' the columns correspond to cells. Can be sparse.
#' @param cellLabels The cell labels of the dataset, which is the gold standard of the trajectory analysis.
#' @param percent Genes that are expressed in less than 100*percent\% of the cells are filtered out. Default is 0.1.
#' @param ncores Number of cores to use. Default is 4.
#' @param need.imputation Whether the input matrix needs imputation. Default is FALSE.
#' @param imputed.data Whether the input matrix has been imputed. Default is TRUE.
#' @return a list of two metrics: Pseudo-temporal Ordering Score (POS) and Kendall's rank correlation score,
#' which are used to assess the accuracy of the inferred cell trajectory and the plot of the cell trajectory.
#'
#' @author Yiqiu Tang
#'
#' @examples
#' library(SAVERg)
#' deng_cellLabels <- factor(colnames(deng_saver),
#'                           levels=c('zygote', 'early 2-cell',
#'                                    'mid 2-cell', 'late 2-cell',
#'                                    '4-cell', '8-cell', '16-cell', 'early blastocyst',
#'                                    'mid blastocyst', 'late blastocyst'))
#' trajectory_analysis(deng_saver, deng_cellLabels)
#'
#' @references
#' Ashenberg, O., Silverbush, D., & Gosik, K. (2019). ANALYSIS OF SINGLE CELL RNA-SEQ DATA.
#' Retrieved from https://broadinstitute.github.io/2019_scWorkshop/.
#'
#' Cole Trapnell Lab (2018). Monocle. Retrieved
#' from http://cole-trapnell-lab.github.io/monocle-release/docs/#constructing-single-cell-trajectories.
#' Huang, M., Zhang, N., & Li, M. (2019, November 13). Single-Cell RNA-Seq Gene Expression Recovery
#'  [R package SAVER version 1.1.2]. Retrieved from https://cran.r-project.org/web/packages/SAVER/index.html.
#'
#'  Huang, M., Zhang, N., & Li, M. (2019, November 13). SAVER Tutorial.
#'   https://cran.r-project.org/web/packages/SAVER/vignettes/saver-tutorial.html.
#'
#' Trajectory Analysis. Retrieved
#' from https://www.rdocumentation.org/packages/geomorph/versions/3.0.7/topics/trajectory.analysis.
#'
#' @export
trajectory_analysis <- function(count_data, cellLabels, percent=0.1, ncores=4,
                                need.imputation=FALSE, imputed.data = TRUE) {
  if (!imputed.data) {
    message("Starting preprocessing ...")
    count_data <- preprocessing(count_data, percent = percent)
    message("Done!")

    if (need.imputation) {
      message("Starting log-normalization ...")
      count_data <- log_normalization(count_data)
      message("Done!")

      message("Starting imputation using SAVER ...")
      count_data <- SAVER::saver(count_data, ncores)$estimate
    }
  }

  message("Starting cell trajectory inference ...")
  colnames(count_data) = c(1:ncol(count_data))
  procdata <- TSCAN::preprocess(count_data)
  lpsmclust <- TSCAN::exprmclust(procdata)
  lpsorder <- TSCAN::TSCANorder(lpsmclust, orderonly=F)
  Pseudotime <- lpsorder$Pseudotime[match(colnames(count_data),lpsorder$sample_name)]
  cor.kendall <- stats::cor(Pseudotime, as.numeric(cellLabels), method = "kendall", use = "complete.obs")
  subpopulation <- data.frame(cell = colnames(count_data), sub = as.numeric(cellLabels)-1)
  POS <- TSCAN::orderscore(subpopulation, lpsorder)[1]
  message("Done!")

  names(POS) <- NULL
  out <- list(cor.kendall=round(cor.kendall, digits = 4),
              POS = round(POS, digits = 4))
  message("Output POS and Kendall's rank correlation score:")
  print(out)
  message("Plot the inferred cell trajectory ...")
  message("Done!")

  g <- TSCAN::plotmclust(lpsmclust, show_cell_names = F)
  print(g)

  # returns
  return(out)
}

# [END]
