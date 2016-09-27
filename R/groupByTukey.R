#' Determine group of genes by all pairwise comparisons (Tukey's method).
#' @noRd
#' @param expr A named list of data frames of significantly
#'  different gene expressions
#'  with genes in rows and samples in columns.
#'  Rows should be named with genes' aliases.
#'  The order of genes has to be the same for each data frame.
#'  The names of list elements should refer to groups of different cases.
#' @param parallel A logical value indicating if a task
#'  should be run on more than one core.
#'  The default value is \code{FALSE}.
#' @param sig.level A numeric value, a significance level used in BH correction
#'  for multiple testing.
#'  The default value is 0.05.
#'
#' @return A data frame with groups of genes named by profiles.
#'  Profiles are determined by the differences between groups.
#'
#' @examples
#' groupByTukey(exanova, sig.level = 0.1)
#'

groupByTukey <- function(expr, parallel = FALSE, sig.level = 0.05) {
    if (!all(sapply(expr, function(x)
        sapply(expr, function(y) row.names(x) == row.names(y))))) {
        stop("Not all genes are named or not all genes
            are in the same order for all data frames,
            or some genes are missing.")
    }
    if (missing(expr) | (is.list(expr) & !is.null(dim(expr))) |
        (length(expr) < 2)) {
        stop("At least two groups are needed for this comparison.")
    }
    if (is.null(names(expr)) | any(names(expr) == "")) {
        names(expr) <- paste0("group", seq_along(expr))
        warning(paste0("Not all elements of the list were named."))
    }
    if (!is.numeric(unlist(expr))) {
        stop("expr can only contain numbers.")
    }
    if (!is.numeric(sig.level) | (sig.level < 0)) {
        sig.level <- 0.05
        warning("No correct significance level specified.
                Taking 0.05 as the sig.level value.")
    }
    tukey.results <- tukeyHSDTest(expr, parallel)
    groupByProfiles(expr, tukey.results, sig.level)
    }
