#' Core function for calculating p-values across all groups of genes
#' using Tukey's Honest Distance test.
#' @noRd
#' @param expr A list of data frames of gene expressions
#'  with genes in rows and samples in columns.
#'  the names of rows represent genes' aliases.
#'  The order of genes has to be the same for each data frame.
#' @param z A number of evaluation of the function needed for
#'  \code{foreach} package.
#' @param group A character name of a group of genes (a cohort).
#'
#' @return A named vector of concatenated: mean value of observations
#'  by groups and the p-value of the Tukey HSD test.

TukeyCore <- function(z, expr, group) {
    n <- length(expr)
    ex <- unlist(sapply(expr, function(x) x[z, ]))
    gene <- cbind.data.frame(ex, group)
    # needed information exactly about all groups
    if (any(sum(tapply(gene$ex, gene$group, length) -
                tapply(gene$ex, gene$group,
                       function(x) sum(is.na(x))) == 0))) {
        rep(NA, n + (n * n - n)/2)
        # a number of groups for which a mean value would be
        # computed plus a number of all possible comparisons in Tukey's test
    } else {
        gene <- na.omit(gene)
        testing <- TukeyHSD(aov(ex ~ group, data = gene))
        mean.val <- tapply(gene$ex, gene$group,
                           mean)
        p.val <- testing$group[seq_len((n * n - n)/2), "p adj"]
        result <- c(mean.val, p.val)
        if (n == 2) {
            names(result) <- "GROUP1-GROUP2"
        }
        result
    }
}

#' Call the calculation of p-values across all groups of genes
#' using Tukey's Honest Distance test optionally on more than
#' one core (actually \code{#(available cores) - 1}).
#' @noRd
#' @param expr A list of data frames of gene expressions with
#'  genes in rows and samples in columns.
#'  The names of rows represent genes' aliases.
#'  The order of genes has to be the same for each data frame.
#' @param parallel A logical value indicating if a task should be
#'  run on more than one core. The default is \code{FALSE}.
#'
#' @return A list of \code{c(mean.val, p.val)}.
#'
#' @import doParallel

tukeyHSDTest <- function(expr, parallel = FALSE) {
    group.names <- names(expr)
    nc <- sapply(expr, ncol)
    if (is.null(nc)) {
        stop("The data are insufficient.
             Provide observations for more than one gene.")
    }
    group <- rep(group.names, nc)
    z <- 0
    if (parallel == TRUE) {
        cores <- detectCores() - 1
        cl <- makeCluster(cores)
        registerDoParallel(cl)
        out <- foreach(z = seq_len(nrow(expr[[1]])),
                       .export = "TukeyCore") %dopar% TukeyCore(z, expr, group)
        stopCluster(cl)
        out
    } else {
        out <- foreach(z = seq_len(nrow(expr[[1]]))) %do%
            TukeyCore(z, expr, group)
    }
    out
}
