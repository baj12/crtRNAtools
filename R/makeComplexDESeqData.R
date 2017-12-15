#' Make a simulated DESeqDataSet with a potentially complex design
#'
#' Constructs a simulated dataset of Negative Binomial data from
#' two conditions. By default, there are no fold changes between
#' the two conditions, but this can be adjusted with the \code{betaSD} argument.
#'
#' @param n number of rows (should be a multiple of 4)
#' @param m number of columns
#' @param betaSD the standard deviation for non-intercept betas, i.e. beta ~ N(0,betaSD)
#' @param interceptMean the mean of the intercept betas (log2 scale)
#' @param interceptSD the standard deviation of the intercept betas (log2 scale)
#' @param dispMeanRel a function specifying the relationship of the dispersions on
#' \code{2^trueIntercept}
#' @param sizeFactors multiplicative factors for each sample
#' @param designStr design string one of "~ condition1 + condition2" ,
#'                  "~ condition1 * condition2", "~ condition1"
#'
#' @return a \code{\link{DESeqDataSet}} with true dispersion,
#' intercept and beta values in the metadata columns.  Note that the true
#' betas are provided on the log2 scale.
#'
#' @examples
#' library(DESeq2)
#' dds <- makeComplexDESeqData()
#' dds
#'
#' @export
makeComplexDESeqData <- function (n = 1000, m = 24, betaSD = 0,
                                  interceptMean = 4, interceptSD = 2,
                                  designStr = "~ condition1 + condition2",
          dispMeanRel = function(x) 4/x + 0.1, sizeFactors = rep(1,m))
{
  beta <- cbind(rnorm(n, interceptMean, interceptSD), rnorm(n,
                                                            0, betaSD))
  dispersion <- dispMeanRel(2^(beta[, 1]))
  colData <- S4Vectors::DataFrame(condition1 = factor(rep(c("A", "B"),
                                              times = c(ceiling(m/2),
                                                        floor(m/2)))),
                       condition2 = factor(rep(c("C", "D","E","F"),
                                               times = m/4)
                                           ))
  x <- if (m > 1) {
    stats::model.matrix.default(~ colData$condition1 )
  }
  else {
    cbind(rep(1, m), rep(0, m))
  }
  mu <- t(2^(x %*% t(beta)) * sizeFactors)
  countData <- matrix(rnbinom(m * n, mu = mu, size = 1/dispersion),
                      ncol = m)
  mode(countData) <- "integer"
  colnames(countData) <- paste("sample", 1:m, sep = "")
  rowRanges <- GenomicRanges::GRanges("1", IRanges::IRanges(start = (1:n - 1) * 100 +
                                      1, width = 100))
  names(rowRanges) <- paste0("gene", 1:n)
  design <- if (m > 1) {
    as.formula(designStr, env = .GlobalEnv)
  }
  else {
    as.formula("~ 1", env = .GlobalEnv)
  }
  object <- DESeq2::DESeqDataSetFromMatrix(countData = countData, colData = colData,
                                   design = design, rowRanges = rowRanges)
  trueVals <- S4Vectors::DataFrame(trueIntercept = beta[, 1], trueBeta = beta[,
                                                                   2], trueDisp = dispersion)
  SummarizedExperiment::mcols(trueVals) <- S4Vectors::DataFrame(type = rep("input", ncol(trueVals)),
                               description = c("simulated intercept values", "simulated beta values",
                                               "simulated dispersion values"))
  SummarizedExperiment::mcols(object) <- cbind(SummarizedExperiment::mcols(object), trueVals)
  return(object)
}
