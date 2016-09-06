#' Rearrange data and select all genes with complete observation in order to analyze it using GOpro package.
#' Prepare an input dataset or use \code{RTCGA.PANCAN12} data.
#'
#' @param RTCGA logical, if TRUE the data is taken from the RTCGA.PANCAN12 package, otherwise provide data using parameter /textit{data} as a vector of data frames. Default value is FALSE.
#' @param cohorts a vector of names of cohorts to be derived from the RTCGA.PANCAN12 package if RTCGA is TRUE, at least two names of: 'colon', 'ovarian', bladder', 'leukemia', 'glioblastoma', 'lung adenocarcinoma', 'BRCA',
#' lung squamous carcinoma', 'kidney', 'endometrioid', 'rectal', 'head and neck', or a vector of names of cohorts provided by data parameter.
#' @param data a list of data frames each with rows named with genes' aliases. The number of a column in each data frame should correspond to the same case.
#'
#' @return A list of data frames of genes expressions or other gene activity related measure with rows named as genes' aliases.
#'
#' @examples
#' prepareData(RTCGA = TRUE, cohorts = c('BRCA', 'ovarian', 'colon'))
#'
#' @export
#'
#' @import  RTCGA.PANCAN12

prepareData <- function(RTCGA = FALSE, cohorts = c('colon', 'bladder', 'leukemia'), data = NULL)
{
  stopifnot(is.logical(RTCGA), is.vector(cohorts), is.character(cohorts))
  if((!is.null(data)) & (RTCGA == FALSE))
  {
    stopifnot(is.list(data), all(sapply(data, ncol) == ncol(data[[1]])), all(sapply(data, nrow) == nrow(data[[1]])), length(data) != length(cohorts))
  }

  if (RTCGA == TRUE)
  {
    library(RTCGA.PANCAN12)
    expr.all <- rbind(expression.cb1, expression.cb2)
    rownames(expr.all) <- NULL
    clinical <- clinical.cb
    TCGAcohortsDict <- c('TCGA Colon Cancer', 'TCGA Bladder Cancer', 'TCGA Acute Myeloid Leukemia', 'TCGA Glioblastoma',
                         'TCGA Ovarian Cancer', 'TCGA Lung Adenocarcinoma', 'TCGA Breast Cancer', 'TCGA Lung Squamous Cell Carcinoma', 'TCGA Kidney Clear Carcinoma',
                         'TCGA Endometrioid Cancer', 'TCGA Rectal Cancer', 'TCGA Head and Neck Cancer')
    names(TCGAcohortsDict) <- c('colon', 'bladder', 'leukemia', 'glioblastoma', 'ovarian', 'lung adenocarcinoma', 'BRCA',
                                'lung squamous carcinoma', 'kidney', 'endometrioid', 'rectal', 'head and neck')
    if (!all(cohorts %in% names(TCGAcohortsDict)))
    {
      stop(paste('The cohorts names specified are not valid cohorts names. Use one of:', paste(names(TCGAcohortsDict), collapse = ', ' )))
    }
    TCGAcohorts <- TCGAcohortsDict[cohorts]
    data <- lapply(TCGAcohorts, function (x) expr.all[ , c(1, na.omit(match(gsub('-', '.', clinical$sampleID[clinical$X_cohort == x]), names(expr.all))))])
  }

  # find and remove rows with a missing observation
  n <- rowSums(sapply(data, function (x) rowSums(apply(x, 2, (is.na)))))
  data <- lapply(data, function(x) x[n == 0, ])

  if(RTCGA == TRUE)
  {
    data <- lapply(data, function(x) {row.names(x) <- as.character(x[ , 1]); x})
    data <- lapply(data, function(x) x[, -1])
  }

  names(data) <- cohorts

  data
}
